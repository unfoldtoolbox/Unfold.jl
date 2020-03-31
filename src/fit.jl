# Mass Univariate Model
function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm,tbl::DataFrame,data::Array{T,3},times;kwargs...) where {T}
    @assert size(data,2) == length(times)

    to = TimerOutput()
    # Generate the designmatrices
    @timeit to "unfoldDesignmatrix" Xs = unfoldDesignmatrix(type,f,tbl;kwargs...)

    # Fit the model
    @timeit to "fit" df = unfoldFit(type,Xs,dropdims(data,dims=1))

    @timeit to "condense" c = condense(df,tbl,times)
    # Condense output and return
    @debug(to)
    return c
end


# Timeexpanded Model
function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm, tbl::DataFrame, data::Array{T,2}, basisfunction::BasisFunction; kwargs...) where {T}
    global to = TimerOutput()

    # Generate the designmatrices
    @timeit to "unfoldDesignmat" Xs = unfoldDesignmatrix(type,f,tbl,basisfunction;kwargs...)

    # Fit the model

    @timeit to "unfoldFit" m = unfoldFit(type,Xs,dropdims(data,dims=2))

    @timeit to "unfoldCondense" c = condense(m,tbl)
    @debug(to)
    return c
end

# helper function for 1 channel data
function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm, tbl::DataFrame, data::Array{T,1}, basisfunction::BasisFunction; kwargs...) where {T}
    @debug("data array is size (X,), reshaping to (1,X)")
    data = reshape(data,:,1)
    fit(type,f,tbl,data,basisfunction;kwargs...)
end


## UnfoldFit functions
# Mass Univariate Linear MOdel
function unfoldFit(::Type{UnfoldLinearModel},X::UnfoldDesignmatrix,data)
    beta,optim = fit_lm(X.Xs, data)
    m = UnfoldLinearModel(beta,optim,X.formulas,X.Xs)
    return m
end


# Massive Univariate Mixed Model
function unfoldFit(::Type{UnfoldLinearMixedModel},X::UnfoldDesignmatrix,data::AbstractArray{T,2}) where {T<:Union{Missing, <:Number}}
    df = Array{LinearMixedModel,1}()
    for t in range(1,stop=size(data,1))
        #println("calculating t $t from $(size(data,1))")
        mm = LinearMixedModel_wrapper(X.formulas,data[t,:],X.Xs)
        fit!(mm)
        push!(df,mm)
    end
    return df
end
# Timeexpanded Mixed Model
function unfoldFit(::Type{UnfoldLinearMixedModel},X::UnfoldDesignmatrix,data::Array{T,1}) where{T}
    mm = LinearMixedModel_wrapper(X.formulas,data,X.Xs)
    fit!(mm)
    return mm
end

function fit_lm(X,data::AbstractArray{T,2}) where {T<:Union{Missing, <:Number}}
    # msas univariate, data = times x epochs
    if size(X,1) > size(data,2)
        X = X[1:size(data,2),:]
        #data[end+1:size(Xs[1],1)] .= 0
    else
        data =data[:,1:size(X,1)]
    end
    # mass univariate
    beta = Array{Union{Missing,Number}}(undef,size(X,2),size(data,1))

    for t in 1:size(data,1)
        #println("$t,$(sum()")
        ix = .!ismissing.(data[t,:])
        beta[:,t] = X[ix,:] \ data[t,ix]
    end

    #beta = X \ data'

    #println(dump(beta))
    history = [] # no history implemented (yet?)
    #beta = dropdims(beta,dims=1) #
    return(beta,[history])

end
function fit_lm(X,data::AbstractArray{T,1}) where {T<:Union{Missing, <:Number}}
    # timeexpanded, data = vector
    # X is epochs x predictor


    # Cut X or y depending on which is larger
    if size(X,1) > size(data,1)
        X = X[1:size(data,1),:]
        #data[end+1:size(Xs[1],1)] .= 0
    else
        data =data[1:size(X,1)]
    end
    ix = .!ismissing.(data)
    # likely much larger matrix, using lsqr
    println(typeof(X[ix,:]))
    println(typeof(data[ix,:]))
    println(eltype(data[ix,:]))
    println(typeof(ix))
    println(eltype(X[ix,:]))
    #println(X[ix,:])
    #println(data[ix])
    beta,history = lsqr(X[ix,:],data[ix],log=true)

    return(beta,history)
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
        reterms = amalgamate(reterms)
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
