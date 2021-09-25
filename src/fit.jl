
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
# helper function encapsulating the call into an array
function StatsModels.fit(UnfoldModelType::Type{T},f::FormulaTerm,tbl::DataFrame,data::AbstractArray,basisOrTimes::Union{BasisFunction,AbstractArray};kwargs...) where{T<:Union{<:UnfoldModel}}
    fit(UnfoldModelType,Dict(Any=>(f,basisOrTimes)),tbl,data;kwargs...)
end


function StatsModels.fit(UnfoldModelType::Type{T},design::Dict,tbl::DataFrame,data::AbstractArray;kwargs...) where {T<:Union{<:UnfoldModel}}
    to = TimerOutput()
    if UnfoldModelType == UnfoldModel
        UnfoldModelType = designToModeltype(design)
    end
    uf = UnfoldModelType(design)

    @timeit to "designmatrix" designmatrix!(uf,tbl;kwargs...)
    @timeit to "fit" fit!(uf,data;kwargs...)

    return uf
end

function StatsModels.fit(UnfoldModelType::Type{T},X::DesignMatrix,data::AbstractArray;kwargs...) where {T<:Union{<:UnfoldModel}}
    if UnfoldModelType == UnfoldModel
        error("Can't infer model automatically, specify with e.g. fit(UnfoldLinearModel...) instead of fit(UnfoldModel...)")
    end
    uf = UnfoldModelType(Dict(),X)
    
    fit!(uf,data;kwargs...)

    return uf
end

isMixedModelFormula(f::ConstantTerm) = false
isMixedModelFormula(f::FormulaTerm) = isMixedModelFormula(f.rhs)

function isMixedModelFormula(f::Tuple)
    ix = [isa(t,FunctionTerm) for t in f]
    return any([isa(t.forig,typeof(|)) for t in f[ix]])
end
function designToModeltype(design)
       # autoDetect
       tmp = collect(values(design))[1]
       f = tmp[1] # formula
       t = tmp[2] # Vector or BasisFunction

        isMixedModel = isMixedModelFormula(f)
        
       if typeof(t) <: BasisFunction
        if isMixedModel
            UnfoldModelType = UnfoldLinearMixedModelContinuousTime
        else
            UnfoldModelType = UnfoldLinearModelContinuousTime
        end
    else
        if isMixedModel
            UnfoldModelType = UnfoldLinearMixedModel
        else
            UnfoldModelType = UnfoldLinearModel
        end
    end
    return UnfoldModelType
end

# helper function for 1 channel data
function StatsModels.fit(UnfoldModelType::Type{T}, design::Dict, tbl::DataFrame, data::AbstractVector,args...; kwargs...) where {T<:Union{<:UnfoldModel}}
    
    @debug("data array is size (X,), reshaping to (1,X)")
    data = reshape(data,1,:)
    return fit(UnfoldModelType,design,tbl,data,args...;kwargs...)
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


function StatsModels.fit!(uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},data::AbstractArray;kwargs...)
    # function content partially taken from MixedModels.jl bootstrap.jl
    df = Array{NamedTuple,1}()
    dataDim = length(size(data)) # surely there is a nicer way to get this but I dont know it

    Xs = modelmatrix(uf)
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
    
    
    mm = LinearMixedModel_wrapper(formula(uf),firstData,Xs)

    # prepare some variables to be used
    βsc, θsc= similar(MixedModels.coef(mm)), similar(mm.θ) # pre allocate
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
    
    uf.modelfit = UnfoldMixedModelFitCollection(
        df,
        deepcopy(mm.λ),
        getfield.(mm.reterms, :inds),
        copy(mm.optsum.lowerbd),
        NamedTuple{Symbol.(fnames(mm))}(map(t -> (t.cnames...,), mm.reterms)),
    )


    return uf.modelfit
end

function StatsModels.coef(uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})
    beta = [x.β for x in MixedModels.tidyβ(modelfit(uf))]
    return reshape_lmm(uf,beta)
end

function MixedModels.ranef(uf::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})
    sigma = [x.σ for x in MixedModels.tidyσs(modelfit(uf))]
    return reshape_lmm(uf,sigma)
end

function reshape_lmm(uf::UnfoldLinearMixedModel,est)
        ntime = length(collect(values(design(uf)))[1][2])
        nchan = modelfit(uf).fits[end].channel
        return permutedims(reshape(est,:,ntime,nchan),[3 2 1])
end
    function reshape_lmm(uf::UnfoldLinearMixedModelContinuousTime,est)
        nchan = modelfit(uf).fits[end].channel
        return reshape(est,:,nchan)'
   
    end





 function StatsModels.fit!(uf::Union{UnfoldLinearModelContinuousTime,UnfoldLinearModel},data;solver=(x,y)->solver_default(x,y),kwargs...)
    @assert ~isempty(designmatrix(uf))
    X = modelmatrix(uf)
    # mass univariate, data = ch x times x epochs
    X,data = zeropad(X,data)

    @debug "UnfoldLinearModel, datasize: $(size(data))"

    uf.modelfit = solver(X,data)

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
