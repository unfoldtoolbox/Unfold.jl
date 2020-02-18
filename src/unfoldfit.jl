function fit(::Type{UnfoldLinearMixedModel},f::FormulaTerm,tbl::DataFrame,data::Array{T,3},times;
       contrasts = Dict{Symbol,Any}()) where {T}
    # yes random, no timeexpansion
    # TODO add contrasts


    tbl,y = dropMissingEpochs(tbl,data)
    println(typeof(y))
    println(typeof(Array{Float64}(y)))
    form  = apply_schema(f,schema(f,tbl,contrasts),LinearMixedModel)

    X = modelcols(form.rhs,tbl)
    #

    df = []
    for t in range(1,stop=size(data,2))

        mm = LinearMixedModel_wrapper(form,y[:,t,:],X)
        fit!(mm)
        push!(df,mm)
    end
    condense.(df,Ref(tbl),times)

end

function fit(::Type{UnfoldLinearMixedModel},f::FormulaTerm,tbl::DataFrame,data::Array{T,2},basisfunction::BasisFunction;    contrasts = Dict{Symbol,Any}()) where {T}
    # yes random, yes timeexpansion
    Xs,form = LinearMixedModel_formula(f,tbl,basisfunction,contrasts=contrasts)
    # TODO perform calculation per channel
    mm = LinearMixedModel_wrapper(form,data,Xs)
    fit!(mm)
    mm = condense(mm,tbl)
    mm
end

function fit(::Type{UnfoldLinearModel},f::FormulaTerm,tbl::DataFrame,y::Array{T,3},times) where {T}
        # no random, no timeexpansion
    form  = apply_schema(f,schema(f,tbl))
    X = modelcols(form.rhs,tbl)

    X,y = dropMissingEpochs(X,y)
    # TODO loop over channels => replaces the dropdims
    beta = X \ dropdims(y,dims=1)'

    m = UnfoldLinearModel(beta,[],form,X)
    return condense(m,tbl,times)

end

function fit(::Type{UnfoldLinearModel},f::FormulaTerm,tbl::DataFrame,y::Array{T,2},basisfunction::BasisFunction) where {T}
    # no random, yes timeexpansion
    form  = apply_schema(f,schema(f,tbl))#+(1|subject)
    formula = FormulaTerm(form.lhs,unfold.TimeExpandedTerm(form.rhs,basisfunction))
    X = modelcols(formula.rhs,tbl)
    y = dropdims(y,dims=1) # XXX change for channels

    if size(X,1)<=length(y)
        y = y[1:size(X,1)]
    elseif size(X,1)>length(y)
        X = X[1:size(y,1),:]
    end

    beta,optim = fit_lm(X,y)
    m = UnfoldLinearModel(beta,optim,formula,X)
    return condense(m,tbl)
end


function fit_lm(X,y)
    # due to our lazy evaluation of X, it could be that it is smaller than y.
    # In that case we fill it with nans

    if size(X,1)<=length(y)
        y = y[1:size(X,1)]
    end

    b,history = lsqr(X,y,log=true)
    return(b,history)
end


function LinearMixedModel_formula(f,tbl,basisfunction;contrasts)


    # TODO This will be updated once LinearMixedModel is refactored a bit
    # TODO: perform missing_omit() after apply_schema() when improved
    # missing support is in a StatsModels release
    #tbl, _ = StatsModels.missing_omit(tbl, f)
    form = apply_schema(f, schema(f, tbl, contrasts), LinearMixedModel)
    form = FormulaTerm(form.lhs, unfold.TimeExpandedTerm.(form.rhs,Ref(basisfunction)))
    _,Xs = modelcols(form, tbl)

    #println(size(data))
    #println(size(Xs[1]))
    return Xs,form
end

function LinearMixedModel_wrapper(
    form,
    y,
    Xs;
    wts = [],
)
# TODO This will be updated once LinearMixedModel is refactored a bit
    y = dropdims(y,dims=1)
    if size(Xs[1],1) > size(y,1)
        println("Timeexpanded designmat longer than data, adding zeros to data. Future versions will fix this")
        target = size(Xs[1],1)-size(y,1)
        push!(y,zeros(target)...)
        #data[end+1:size(Xs[1],1)] .= 0
    else
        y = y[1:size(Xs[1],1)]
    end
    y = (reshape(float(y), (:, 1)))
    T = eltype(y)

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

    push!(feterms, MixedModels.FeMat(y, [""]))

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
