function condense(m,tbl,colnames_basis)
    # no random effects no Timeexpansion
    cnames = coefnames(m.formula.rhs)
    if typeof(cnames)<:String
        cnames = [cnames]
    end

    cnames_rep = repeat(cnames,length(colnames_basis))

    colnames_basis_rep = repeat(colnames_basis,1,length(cnames))
    colnames_basis_rep = dropdims(reshape(colnames_basis_rep',:,1),dims=2)


    betas = dropdims(reshape(m.beta,:,1),dims=2)

    results = DataFrame(term=cnames_rep,estimate=betas,stderror=Missing,group="mass univariate",colnames_basis=colnames_basis_rep)
    return UnfoldModel(m,m.formula,tbl,results)

end

function condense(mm_array::Array{LinearMixedModel,1},tbl,colnames_basis)
    # with random effects, no timeexpansion

    results = condense_fixef.(mm_array,colnames_basis)

    results = vcat(results,condense_ranef.(mm_array,colnames_basis))
    # XXX TODO Return an Array of Models, not only the first one!
    return UnfoldModel(mm_array[1],mm_array[1].formula,tbl,vcat(results...))

end



function condense(mm::UnfoldLinearModel,tbl)
    # no random effects, timeexpansion

    if typeof(mm.formula)<: Array
        # build new UnfoldLinearModel for each formula

        colnames_basis = [formula.rhs.basisfunction.colnames for formula in mm.formula]
        fromTo = [0 cumsum(length.(colnames_basis),dims=2)]

        results = DataFrame()
        for k = 1:length(mm.formula)
            mm_tmp = UnfoldLinearModel(mm.beta[fromTo[k]+1:fromTo[k+1]],mm.optim,mm.formula[k],mm.X[:,fromTo[k]+1:fromTo[k+1]])
            res = condense_fixef(mm_tmp,colnames_basis[k])
            results = vcat(results,res)
        end



    else
        colnames_basis = mm.formula.rhs.basisfunction.colnames

        results = condense_fixef(mm,colnames_basis)
    end

    return UnfoldModel(mm,mm.formula,tbl,results)
end


function condense(mm::LinearMixedModel,tbl)
    # with random effects, timeexpansion
    # TODO loop over al timeexpanded basisfunctions, in case there are multiple ones in case of multiple events
    # TODO random correlations => new function in MixedModels
    # TODO somehow get rid of the split of the coefficient names. What if there is a ":" in a coefficient name?
    colnames_basis = mm.formula.rhs[1].basisfunction.colnames
    results = condense_fixef(mm,colnames_basis)

    # ranefs more complex
    results = vcat(results,condense_ranef(mm,colnames_basis))

    #return
    return UnfoldModel(mm,mm.formula,tbl,results)
end


function condense_fixef(mm,colnames_basis)
    if typeof(mm.formula.rhs) <: Tuple
        fixefPart = mm.formula.rhs[1]
    else
        fixefPart = mm.formula.rhs
    end

    cnames = [c[1] for c in split.(coefnames(fixefPart)," :")]

    colnames_basis =  repeat(colnames_basis,length(unique(cnames)))


    return DataFrame(term=cnames,estimate=MixedModels.fixef(mm),stderror=MixedModels.stderror(mm),group="fixed",colnames_basis=colnames_basis)
end

function MixedModels.fixef(m::UnfoldLinearModel)
    # condense helper for the linear model which just returns a vector of fixef
     return m.beta
end
function MixedModels.stderror(m::UnfoldLinearModel)
    # for now we don't have an efficient way to calculate SEs for single subjects
    return fill(NaN,size(m.beta))
end
function condense_ranef(mm,colnames_basis)
    vc = VarCorr(mm)
    σρ = vc.σρ

    cnames = string.(foldl(vcat, [keys(sig)...] for sig in getproperty.(values(σρ), :σ)))
    cnames = [c[1] for c in split.(cnames," :")]

    σvec = vcat(collect.(values.(getproperty.(values(σρ), :σ)))...)


    nmvec = string.([keys(σρ)...])
    nvec = length.(keys.(getproperty.(values(σρ),:σ)))
    group = []
    for n in zip(nmvec,nvec)
        append!(group,repeat([n[1]],n[2]))
    end

    colnames_basis =  repeat(colnames_basis,length(unique(cnames)))
    # combine
    return DataFrame(term=cnames,estimate=σvec,stderror=NaN,group=group,time=colnames_basis)

end
