
"""
fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm,tbl::DataFrame,data::Array{T,3},times)
fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm,tbl::DataFrame,data::Array{T,2},basisfunction::BasisFunction)

Generates Designmatrix & fits model, either mass-univariate (one model per epoched-timepoint) or time-expanded (modeling linear overlap).


# Examples
Mass Univariate Linear
```julia-repl
julia> data,evts = loadtestdata("testCase1") 
julia> data_r = reshape(data,(1,:))
julia> data_e,times = unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1.9),sfreq=10) # cut the data into epochs. data_e is now ch x times x epoch

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


# Timeexpanded Model
function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm, tbl::DataFrame, data::Array{T,2}, basisfunction::BasisFunction; kwargs...) where {T}
    to = TimerOutput()

    # Generate the designmatrices
    @timeit to "unfoldDesignmat" Xs = designmatrix(type,f,tbl,basisfunction;kwargs...)

    # Fit the model
    @timeit to "unfoldfit" ufModel = unfoldfit(type,Xs,data)

    @timeit to "unfoldCondense" c = condense_long(ufModel)
    @debug(to)
    return ufModel,c
end

# helper function for 1 channel data
function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm, tbl::DataFrame, data::Array{T,1}, basisfunction::BasisFunction; kwargs...) where {T}
    @debug("data array is size (X,), reshaping to (1,X)")
    data = reshape(data,1,:)
    return fit(type,f,tbl,data,basisfunction;kwargs...)
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

    _,data = zeropad(Xobj.Xs[1],data)
    print(size(data))
    print(size(Xobj.Xs[1]))
    # get a un-fitted mixed model object
    mm = LinearMixedModel_wrapper(Xobj.formulas,firstData,Xobj.Xs)

    # prepare some variables to be used
    βsc, θsc= similar(coef(mm)), similar(mm.θ) # pre allocate
    p,k = length(βsc), length(θsc)
    β_names = (Symbol.(fixefnames(mm))..., )

    # for each channel
    @showprogress .1 for ch in range(1,stop=nchan)
        # for each time
        for t in range(1,stop=ntime)

            @debug print("ch:$ch/$nchan, t:$t/$ntime")
            if ndims(data) == 3
                refit!(mm,data[ch,t,:])
            else
                refit!(mm,data[ch,:])
            end

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
    if size(Xs[1],1) > size(data,1)
        @warn("Timeexpanded designmat longer than data, adding zeros to data. Future versions will fix this")
        target = size(Xs[1],1)-size(data,1)
        push!(data,zeros(target)...)
        #data[end+1:size(Xs[1],1)] .= 0
    else
        data = data[1:size(Xs[1],1)]
    end

    # TODO what follows will be updated and potentially replaced when LinearMixedModel is refactored
    # It is copied 0.98:1 from MixedModels.jl (the missing part is new)

    data= (reshape(float(data), (:, 1)))
    T = eltype(TData)

    reterms = ReMat{T}[]
    feterms = FeMat{T}[]
    for (i, x) in enumerate(Xs)

        if isa(x, ReMat{T})
            push!(reterms, x)
        else
            cnames = coefnames(form.rhs[i])
            push!(feterms, MixedModels.FeMat(x, isa(cnames, String) ? [cnames] : collect(cnames)))
        end
    end

    push!(feterms, MixedModels.FeMat(data, [""]))

    # detect and combine RE terms with the same grouping var
    if length(reterms) > 1
        reterms = MixedModels.amalgamate(reterms)
    end

    sort!(reterms, by = MixedModels.nranef, rev = true)
    allterms = convert(Vector{Union{ReMat{T},FeMat{T}}}, vcat(reterms, feterms))

    A, L = MixedModels.createAL(allterms)
    lbd = foldl(vcat, MixedModels.lowerbd(c) for c in reterms)
    θ = foldl(vcat, MixedModels.getθ(c) for c in reterms)
    optsum = OptSummary(θ, lbd, :LN_BOBYQA, ftol_rel = T(1.0e-12), ftol_abs = T(1.0e-8))
    fill!(optsum.xtol_abs, 1.0e-10)
    X = first(feterms)
    LinearMixedModel(form, allterms, sqrt.(convert(Vector{T}, wts)), MixedModels.mkparmap(reterms),    (n = size(X, 1), p = size(X,2), nretrms = length(reterms)), A, L, optsum)
end


"""
$(SIGNATURES)

Type Piracy. can be removed once MixedModels fully supports sparse FeMats
https://github.com/JuliaStats/MixedModels.jl/pull/309

"""
## Piracy it is!
function MixedModels.FeMat(X::SparseMatrixCSC, cnms)
    #println("Pirated Sparse FeMat")
    FeMat{eltype(X),typeof(X)}(X, X, range(1,stop=size(X,2)), minimum(size(X)), cnms)
end
