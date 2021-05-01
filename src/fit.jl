
"""
fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm,tbl::DataFrame,data::Array{T,3},times)
fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm,tbl::DataFrame,data::Array{T,2},basisfunction::BasisFunction)

Generates Designmatrix & fits model, either mass-univariate (one model per epoched-timepoint) or time-expanded (modeling linear overlap).


# Examples
Mass Univariate Linear
```julia-repl
julia> data,evts = loadtestdata("testCase1")
julia> data_r = reshape(data,(1,:))
julia> data_e,times = Unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1.9),sfreq=10) # cut the data into epochs. data_e is now ch x times x epoch

julia> f  = @formula 0~1+continuousA+continuousB # 1
julia> model,results_long = fit(UnfoldLinearModel,f,evts,data_e,times)
```
Timexpanded Univariate Linear
```julia-repl
julia> basisfunction = firbasis(τ=(-1,1),sfreq=10,name="A")
julia> model,results_long = fit(UnfoldLinearModel,f,evts,data_r,basisfunction)
```

"""
function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm,tbl::DataFrame,data::Array{T,3},times;kwargs...) where {T}
    @assert size(data,2) == length(times)

    to = TimerOutput()
    # Generate the designmatrices
    @timeit to "designmatrix" Xs = designmatrix(type,f,tbl;kwargs...)

    # Fit the model
    @timeit to "fit" ufModel = unfoldfit(type,Xs,data;kwargs...)

    @timeit to "condense" c = condense_long(ufModel,times)
    # Condense output and return
    @debug(to)
    return ufModel,c
end


function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm, tbl::DataFrame, data::Array{T,2},basisfunction::BasisFunction; kwargs...) where {T}
# backward compatible call, typically these models are run with multiple events
fit(type,Dict(Any=>(f,basisfunction)), tbl, data; kwargs...)
end

# Timeexpanded Model
function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::Dict, tbl::DataFrame, data::Array{T,2}; kwargs...) where {T}
    to = TimerOutput()

    # Generate the designmatrices
    @timeit to "unfoldDesignmat" Xs = designmatrix(type,f,tbl;kwargs...)

    # Fit the model
    @timeit to "unfoldfit" ufModel = unfoldfit(type,Xs,data;kwargs...)

    @timeit to "unfoldCondense" c = condense_long(ufModel)
    @debug(to)
    return ufModel,c
end

# helper function for 1 channel data
function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f, tbl::DataFrame, data::Array{T,1}, args...; kwargs...) where {T}
    @debug("data array is size (X,), reshaping to (1,X)")
    data = reshape(data,1,:)
    return fit(type,f,tbl,data,args...;kwargs...)
end



"""
unfoldfit(::Type{UnfoldLinearMixedModel},Xobj::DesignMatrix,data::Union{<:AbstractArray{T,2},<:AbstractArray{T,3}}) where {T<:Union{Missing, <:Number}}
unfoldfit(::Type{UnfoldLinearModel},Xobj::DesignMatrix,data::Union{<:AbstractArray{T,2},<:AbstractArray{T,3}}) where {T<:Union{Missing, <:Number}}

Fit a DesignMatrix against a 2D/3D Array data along its last dimension
Data is typically interpreted as channel x time (with basisfunctions) or channel x time x epoch (for mass univariate)

Returns an UnfoldModel object with fields .beta, .sigma (in case of mixed model) and .modelinfo with a summary of the modelfit

Note: Might be renamed/refactored to fit! at a later point

# Examples
```julia-repl
```

"""
function unfoldfit(::Type{UnfoldLinearMixedModel},Xobj::DesignMatrix,data::Union{<:AbstractArray{T,2},<:AbstractArray{T,3}};kwargs...) where {T<:Union{Missing, <:Number}}
    # function content partially taken from MixedModels.jl bootstrap.jl
    df = Array{NamedTuple,1}()
    dataDim = length(size(data)) # surely there is a nicer way to get this but I dont know it

    # If we have3 dimension, we have a massive univariate linear mixed model for each timepoint
    if dataDim == 3
        firstData = data[1,1,:]
        ntime = size(data,2)
    else
        # with only 2 dimension, we run a single time-expanded linear mixed model per channel/voxel
        firstData = data[1,:]
        ntime = 1
    end
    nchan = size(data,1)

    #_,data = zeropad(Xobj.Xs[1],data)
    # get a un-fitted mixed model object
    
    
    mm = LinearMixedModel_wrapper(Xobj.formulas,firstData,Xobj.Xs)

    # prepare some variables to be used
    βsc, θsc= similar(coef(mm)), similar(mm.θ) # pre allocate
    p,k = length(βsc), length(θsc)
    #β_names = (Symbol.(fixefnames(mm))..., )
    
    β_names = (Symbol.(vcat(fixefnames(mm)...))...,)
    β_names = (unique(β_names)...,)

    @assert(length(β_names) == length(βsc),"Beta-Names & coefficient length do not match. Did you provide two identical basis functions?")

    @debug println("beta_names $β_names")
    @debug println("uniquelength: $(length(unique(β_names))) / $(length(β_names))")
    # for each channel
    prog = Progress(nchan*ntime,.1)
    #@showprogress .1 
    for ch in range(1,stop=nchan)
        # for each time
        for t in range(1,stop=ntime)

            #@debug "ch:$ch/$nchan, t:$t/$ntime"
            @debug "data-size: $(size(data))"
            #@debug println("mixedModel: $(mm.feterms)")
            if ndims(data) == 3
                refit!(mm,data[ch,t,:])
            else
                refit!(mm,data[ch,:])
            end
            #@debug println(MixedModels.fixef!(βsc,mm))
            println(length(βsc))
            β = NamedTuple{β_names}(MixedModels.fixef!(βsc, mm))

            out = (
            objective = mm.objective,
            σ = mm.σ,
            β = NamedTuple{β_names}(MixedModels.fixef!(βsc, mm)),
            se = SVector{p,Float64}(MixedModels.stderror!(βsc, mm)), #SVector not necessary afaik, took over from MixedModels.jl
            θ = SVector{k,Float64}(MixedModels.getθ!(θsc, mm)),
            channel = ch,
            timeIX = ifelse(dataDim==2,NaN,t)
            )
            push!(df,out)
            ProgressMeter.next!(prog; showvalues = [(:channel,ch), (:time,t)])
        end
    end
    
    modelinfo = MixedModelBootstrap(
        df,
        deepcopy(mm.λ),
        getfield.(mm.reterms, :inds),
        copy(mm.optsum.lowerbd),
        NamedTuple{Symbol.(fnames(mm))}(map(t -> (t.cnames...,), mm.reterms)),
    )

    beta = [x.β for x in MixedModels.tidyβ(modelinfo)]

    sigma = [x.σ for x in MixedModels.tidyσs(modelinfo)]

    if ndims(data)==3
        beta = permutedims(reshape(beta,:,ntime,nchan),[3 2 1])
        sigma = permutedims(reshape(sigma,:,ntime,nchan),[3 2 1])
    else
        beta = reshape(beta,:,nchan)'
        sigma = reshape(sigma,:,nchan)'
    end



    return UnfoldLinearMixedModel(beta,sigma,modelinfo,Xobj)
end

 function unfoldfit(::Type{UnfoldLinearModel},Xobj::DesignMatrix,data;solver=(x,y)->solver_default(x,y),kwargs...)
     X = Xobj.Xs
    # mass univariate, data = ch x times x epochs
    X,data = zeropad(X,data)

    @debug "UnfoldLinearModel, datasize: $(size(data))"
    # mass univariate
    beta,modelinfo = solver(X,data)

    return UnfoldLinearModel(beta,modelinfo,Xobj)

end



"""
$(SIGNATURES)

Wrapper to generate a LinearMixedModel. Code taken from MixedModels.jl and slightly adapted.

"""
function LinearMixedModel_wrapper(form,data::Array{<:Union{TData},1},Xs;wts = []) where {TData<:Number}
    #    function LinearMixedModel_wrapper(form,data::Array{<:Union{Missing,TData},1},Xs;wts = []) where {TData<:Number}


    # XXX Push this to utilities zeropad
    # Make sure X & y are the same size
    m = size(Xs[1])[1]


    if m != size(data)[1]
        Xs = changeMatSize!(size(data)[1],Xs[1],Xs[2:end])
    end
    
    y = (reshape(float(data), (:, 1)))

    MixedModels.LinearMixedModel(y, Xs, form, wts)
 end

 function MixedModels.LinearMixedModel(y, Xs, form::Array, wts)
    

    form_combined = form[1]
    for f =form[2:end]
        
        form_combined = form_combined.lhs ~ MatrixTerm(form_combined.rhs[1] + f.rhs[1]) + form_combined.rhs[2:end] + f.rhs[2:end]
    end
    MixedModels.LinearMixedModel(y, Xs, form_combined, wts)
 end
