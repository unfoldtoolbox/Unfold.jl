
function condense(mm,tbl)

    # TODO loop over al timeexpanded basisfunctions, in case there are multiple ones in case of multiple events
    # TODO random correlations
    # TODO somehow get rid of the split of the coefficient names. What if there is a ":" in a coefficient name?
    if typeof(mm) == MixedModel
        fixefPart = mm.formula.rhs[1]
    else
        fixefPart = mm.formula.rhs
    end
    cnames = [c[1] for c in split.(coefnames(fixefPart)," :")]
    times =  repeat(fixefPart.basisfunction.times,length(unique(cnames)))

    # fixefs are easy
    results = DataFrame(term=cnames,estimate=fixef(mm),stderror=stderror(mm),group="fixed",time=times)

    # ranefs more complex
    if typeof(mm) == MixedModel
        cnames,σvec,group = condense_ranef(mm)
        times =  repeat(fixefPart.basisfunction.times,length(unique(cnames)))
        # combine
        results = vcat(results,DataFrame(term=cnames,estimate=σvec,stderror=NaN,group=group,time=times))
    end
    #return
    return UnfoldModel(mm,mm.formula,tbl,results)
end


function fixef(m::UnfoldLinearModel)
    # condense helper for the linear model which just returns a vector of fixef
     return m.beta
end
function stderror(m::UnfoldLinearModel)
    # for now we don't have an efficient way to calculate SEs for single subjects
    return fill(NaN,size(m.beta))
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
