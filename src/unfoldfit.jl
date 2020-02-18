function fit(::Type{UnfoldLinearMixedModel},f::FormulaTerm,tbl::DataFrame,data,basisfunction::BasisFunction;    contrasts = Dict{Symbol,Any}())
    Xs,data,form = LinearMixedModel_formula(f,tbl,data,basisfunction,contrasts=contrasts)
    mm = LinearMixedModel_wrapper(form,data,Xs)
    fit!(mm)
    mm = condense(mm,tbl)
    mm
end

function LinearMixedModel_formula(f,tbl,data,basisfunction;contrasts)
    # TODO: perform missing_omit() after apply_schema() when improved
    # missing support is in a StatsModels release
    #tbl, _ = StatsModels.missing_omit(tbl, f)
    form = apply_schema(f, schema(f, tbl, contrasts), LinearMixedModel)
    form = FormulaTerm(form.lhs, unfold.TimeExpandedTerm.(form.rhs,Ref(basisfunction)))
    _,Xs = modelcols(form, tbl)
    data = dropdims(data,dims=1)
    if size(Xs[1],1) > size(data,1)
        println("Timeexpanded designmat longer than data, adding zeros to data. Future versions will fix this")
        target = size(Xs[1],1)-size(data,1)
        push!(data,zeros(target)...)
        #data[end+1:size(Xs[1],1)] .= 0
    else
        data = data[1:size(Xs[1],1)]
    end
    println(size(data))
    println(size(Xs[1]))
    return Xs,data,form
end

function fit(::Type{UnfoldLinearModel},f::FormulaTerm,tbl::DataFrame,y::AbstractVector,basisfunction::BasisFunction)
    #TODO Allow data to be matrix instead of vector
    form  = apply_schema(f,schema(f,tbl))#+(1|subject)
    formula = FormulaTerm(form.lhs,unfold.TimeExpandedTerm(form.rhs,basisfunction))
    X = modelcols(formula.rhs,tbl)


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



function LinearMixedModel_wrapper(
    form,
    y,
    Xs;
    wts = [],
)


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
        println("Priated Sparse FeMat")
        FeMat{eltype(X),typeof(X)}(X, X, range(1,stop=size(X,2)), minimum(size(X)), cnms)
end
