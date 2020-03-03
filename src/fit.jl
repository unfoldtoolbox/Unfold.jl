# Mass Univariate Model
function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm,tbl::DataFrame,data::Array{T,3},times;kwargs...) where {T}
    @assert size(data,2) == length(times)

    to = TimerOutput()
    # Generate the designmatrices
    @timeit to "unfoldDesignmatrix" Xs = unfoldDesignmatrix(type,f,tbl;kwargs...)

    # Fit the model
    @timeit to "fit" df = unfoldFit(type,Xs,dropdims(data,dims=1))

    # Condense output and return
    @timeit to "condense" c = condense(df,tbl,times)
    println(to)
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
    println(to)
    return c
end

# helper function for 1 channel data
function fit(type::Type{<:Union{UnfoldLinearModel,UnfoldLinearMixedModel}},f::FormulaTerm, tbl::DataFrame, data::Array{T,1}, basisfunction::BasisFunction; kwargs...) where {T}
    println("data array is size (X,), reshaping to (1,X)")
    data = reshape(data,:,1)
    fit(type,f,tbl,data,basisfunction;kwargs...)
end


## UnfoldFit functions
# Mass Univariate Linear MOdel
function unfoldFit(::Type{UnfoldLinearModel},X::UnfoldDesignmatrix,data::Array{T,2}) where {T}
    beta,optim = fit_lm(X.Xs, data)
    m = UnfoldLinearModel(beta,optim,X.formulas,X.Xs)
    return m
end
# Timeexpanded Model
function unfoldFit(::Type{UnfoldLinearModel},X::UnfoldDesignmatrix,data::Array{T,1}) where {T}
    beta,optim = fit_lm(X.Xs,data)
    m = UnfoldLinearModel(beta,optim,X.formulas,X.Xs)
    return m
end

# Massive Univariate Mixed Model
function unfoldFit(::Type{UnfoldLinearMixedModel},X::UnfoldDesignmatrix,data::Array{T,2}) where {T}
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

function fit_lm(X,data::Array{T,2}) where {T}
    # msas univariate, data = times x epochs
    if size(X,1) > size(data,2)
        X = X[1:size(data,2),:]
        #data[end+1:size(Xs[1],1)] .= 0
    else
        data =data[:,1:size(X,1)]
    end
    println("mass univariate case")
    # mass univariate
    beta = X \ data'
    history = [] # no history implemented (yet?)
    #beta = dropdims(beta,dims=1) #
    return(beta,[history])

end
function fit_lm(X,data::Array{T,1}) where {T}
    # timeexpanded, data = vector
    # X is epochs x predictor
    println("fit_lm \n, data $(size(data))\n X $(size(X)) \n)")

    # Cut X or y depending on which is larger
    if size(X,1) > size(data,1)
        X = X[1:size(data,1),:]
        #data[end+1:size(Xs[1],1)] .= 0
    else
        data =data[1:size(X,1)]
    end
    println(size(data))

    println("time expanded case")
    # likely much larger matrix, using lsqr
    beta,history = lsqr(X,data,log=true)

    println("fitlm")
    println(size(beta))
    return(beta,history)
end


function LinearMixedModel_wrapper(form,data::Array{T2,1},Xs;wts = []) where {T2}
    # Make sure X & y are the same size
    if size(Xs[1],1) > size(data,1)
        println("Timeexpanded designmat longer than data, adding zeros to data. Future versions will fix this")
        target = size(Xs[1],1)-size(data,1)
        push!(data,zeros(target)...)
        #data[end+1:size(Xs[1],1)] .= 0
    else
        data = data[1:size(Xs[1],1)]
    end

    # TODO what follows will be updated and potentially replaced when LinearMixedModel is refactored
    # It is copied 1:1 from MixedModels.jl

    data= (reshape(float(data), (:, 1)))
    T = eltype(data)

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
