function fit(::Type{UnfoldLinearMixedModel},f::FormulaTerm,tbl::DataFrame,data,basisfunction::BasisFunction;    contrasts = Dict{Symbol,Any}())
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
        println(target)
        println(size(data))
        println(size(Xs[1]))
        push!(data,zeros(target)...)
        #data[end+1:size(Xs[1],1)] .= 0
    else
        data = data[1:size(Xs[1],1)]
    end
    println(size(data))
    println(size(Xs[1]))
    mm = LinearMixedModel_wrapper(form,data,Xs)
    fit!(mm)
    mm = condense(mm,tbl)
end

function condense(mm,tbl)

    # TODO loop over al timeexpanded basisfunctions, in case there are multiple ones in case of multiple events
    # TODO random correlations
    # TODO somehow get rid of the split of the coefficient names. What if there is a ":" in a coefficient name?
    cnames = [c[1] for c in split.(coefnames(mm.formula.rhs[1])," :")]
    times =  repeat(mm.formula.rhs[1].basisfunction.times,length(unique(cnames)))

    # fixefs are easy
    results = DataFrame(term=cnames,estimate=fixef(mm),stderror=stderror(mm),group="fixed",time=times)

    # ranefs more complex
    cnames,σvec,group = condense_ranef(mm)
    # combine
    results = vcat(results,DataFrame(term=cnames,estimate=σvec,stderror=NaN,group=group,time=times))

    #return
    UnfoldModel(mm,mm.formula,tbl,results)
end


function condense_ranef(mm)
    vc = VarCorr(mm)
    σρ = vc.σρ

    cnames = string.(foldl(vcat, [keys(sig)...] for sig in getproperty.(values(σρ), :σ)))
    cnames = [c[1] for c in split.(cnames," :")]

    σvec = vcat(collect.(values.(getproperty.(values(σρ), :σ)))...)


    nmvec = string.([keys(σρ)...])
    nvec = length.(keys.(getproperty.(values(σρ),:σ)))
    nmvec2 = []
    for n in zip(nmvec,nvec)
        append!(nmvec2,repeat([n[1]],n[2]))
    end
    return cnames,σvec,nmvec2
end


function fit(::Type{UnfoldLinearModel},f::FormulaTerm,tbl::DataFrame,y::AbstractVector,basisfunction::BasisFunction)
    #TODO Allow data to be matrix instead of vector
    formula  = apply_schema(f,schema(f,tbl))#+(1|subject)
    formula = unfold.TimeExpandedTerm(formula.rhs,basisfunction)
    X = modelcols(formula,tbl)
    beta,optim = fit_lm(X,y)
    # "cut" betas & change to a dataFrame
    #model = UnfoldModel(XXX,formula,tbl,results)
    #return(model)
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
function MixedModels.FeMat(X::SparseMatrixCSC, cnms)
        println("sparseFeMat")
        FeMat{eltype(X),typeof(X)}(X, X, range(1,stop=size(X,2)), minimum(size(X)), cnms)
end
