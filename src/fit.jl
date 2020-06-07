# Mass Univariate Model
function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm,tbl::DataFrame,data::Array{T,3},times;kwargs...) where {T}
    @assert size(data,2) == length(times)

    to = TimerOutput()
    # Generate the designmatrices
    @timeit to "designmatrix" Xs = designmatrix(type,f,tbl;kwargs...)

    # Fit the model
    @timeit to "fit" ufModel = unfoldfit(type,Xs,data)

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


## actual fitting functions

function unfoldfit(::Type{UnfoldLinearMixedModel},Xobj::DesignMatrix,data)#AbstractArray{T,3} where {T<:Union{Missing, <:Number}}
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
    for ch in range(1,stop=nchan)
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

 function unfoldfit(::Type{UnfoldLinearModel},Xobj::DesignMatrix,data::AbstractArray{T,3};optimizer=undef) where {T<:Union{Missing, <:Number}}
     X = Xobj.Xs
    # mass univariate, data = ch x times x epochs
    X,data = zeropad(X,data)

    print(size(data))
    # mass univariate
    beta = Array{Union{Missing,Number}}(undef,size(data,1),size(data,2),size(X,2))
    for ch in 1:size(data,1)
        for t in 1:size(data,2)
            @debug("$(ndims(data,)),$t,$ch")
            ix = .!ismissing.(data[ch,t,:])
            beta[ch,t,:] = X[ix,:] \ data[ch,t,ix]
        end
    end

    modelinfo = [undef] # no history implemented (yet?)

    return UnfoldLinearModel(beta,modelinfo,Xobj)

end

"""
$(SIGNATURES)
Equates the length of data and designmatrix by cutting the shorter one

The reason we need this is because when generating the designmatrix, we do not know how long the data actually are. We only assume that event-latencies are synchronized with the data
"""
function zeropad(X,data::AbstractArray{T,2})where {T<:Union{Missing, <:Number}}
    @debug("2d zeropad")
    if size(X,1) > size(data,2)
        X = X[1:size(data,2),:]
    else
        data =data[:,1:size(X,1)]
    end
    return X,data
end
function zeropad(X,data::AbstractArray{T,3})where {T<:Union{Missing, <:Number}}
    @debug("3d zeropad")
    if size(X,1) > size(data,3)
        X = X[1:size(data,3),:]
    else
        data =data[:,:,1:size(X,1)]
    end
    return X,data
end
function unfoldfit(::Type{UnfoldLinearModel},Xobj::DesignMatrix,data::AbstractArray{T,2};optimizer=undef) where {T<:Union{Missing, <:Number}}
    # timeexpanded, data = ch x time
    # X is epochs x predictor

    X = Xobj.Xs # extract designmatrix
    # Cut X or y depending on which is larger
    X,data = zeropad(X,data)


    # likely much larger matrix, using lsqr
    modelinfo = []
    beta = Array{Float64}(undef,size(data,1),size(X,2))

    for ch in 1:size(data,1)
        ix = .!ismissing.(data[ch,:])
        beta[ch,:],h = lsqr(X[ix,:],data[ch,ix],log=true)
        push!(modelinfo,h)
    end

    return UnfoldLinearModel(beta,modelinfo,Xobj)
end


function LinearMixedModel_wrapper(form,data::Array{<:Union{TData},1},Xs;wts = []) where {TData<:Number}
    #    function LinearMixedModel_wrapper(form,data::Array{<:Union{Missing,TData},1},Xs;wts = []) where {TData<:Number}


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
    LinearMixedModel(form, allterms, sqrt.(convert(Vector{T}, wts)), A, L, optsum)
end


## Piracy it is!
function MixedModels.FeMat(X::SparseMatrixCSC, cnms)
    #println("Pirated Sparse FeMat")
    FeMat{eltype(X),typeof(X)}(X, X, range(1,stop=size(X,2)), minimum(size(X)), cnms)
end
