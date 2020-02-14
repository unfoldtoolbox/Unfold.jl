function fit(::Type{UnfoldLinearMixedModel},f::FormulaTerm,tbl::DataFrame,data,basisfunction::BasisFunction;    contrasts = Dict{Symbol,Any}(),
)
    # TODO: perform missing_omit() after apply_schema() when improved
    # missing support is in a StatsModels release
    #tbl, _ = StatsModels.missing_omit(tbl, f)
    form = apply_schema(f, schema(f, tbl, contrasts), LinearMixedModel)
    form = FormulaTerm(form.lhs, unfold.TimeExpandedTerm.(form.rhs,Ref(basisfunction)))
    _,Xs = modelcols(form, tbl)
    data = data[1:size(Xs[1],1)]
    LinearMixedModel_wrapper(form,data,Xs)

end

function fit(::Type{UnfoldLinearModel},f::FormulaTerm,tbl::DataFrame,y::AbstractVector,basisfunction::BasisFunction)
    #TODO Allow data to be matrix instead of vector
    formula  = apply_schema(f,schema(f,tbl))#+(1|subject)
    formula = unfold.TimeExpandedTerm(formula.rhs,basisfunction)
    X = modelcols(formula,tbl)
    beta,optim = fit_lm(X,y)
    # "cut" betas & change to a dataFrame
    model = UnfoldLinearModel(formula,tbl,X,beta,optim)
    return(model)
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



function LinearMixedModel_wrapper(
    form,
    y,
    Xs;
    wts = [],
)


    y = SparseMatrixCSC(reshape(float(y), (:, 1)))
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
    println("pushing y",typeof(y))
    push!(feterms, MixedModels.FeMat(y, [""]))

    # detect and combine RE terms with the same grouping var
    if length(reterms) > 1
        reterms = amalgamate(reterms)
    end

    sort!(reterms, by = MixedModels.nranef, rev = true)
    println("reterms: $(size(reterms))")
    println("feterms: $(size(feterms))")
    allterms = convert(Vector{Union{ReMat{T},FeMat{T}}}, vcat(reterms, feterms))

    A, L = MixedModels.createAL(allterms)
    lbd = foldl(vcat, MixedModels.lowerbd(c) for c in reterms)
    θ = foldl(vcat, MixedModels.getθ(c) for c in reterms)
    optsum = OptSummary(θ, lbd, :LN_BOBYQA, ftol_rel = T(1.0e-12), ftol_abs = T(1.0e-8))
    fill!(optsum.xtol_abs, 1.0e-10)
    LinearMixedModel(form, allterms, sqrt.(convert(Vector{T}, wts)), A, L, optsum)
end


## Piracy it is!
function FeMat(X::SparseMatrixCSC, cnms)
        println("sparseFeMat")
        FeMat{eltype(X),typeof(X)}(X, X, [], maximum(size(X)), cnms)
end
