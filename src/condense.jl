function condense(m,tbl,times)
    # no random effects no Timeexpansion
    cnames = coefnames(m.formula.rhs)
    cnames_rep = repeat(cnames,length(times))

    times_rep = repeat(times,1,length(cnames))
    times_rep = dropdims(reshape(times_rep',:,1),dims=2)

    println(size(m.beta))
    betas = dropdims(reshape(m.beta,:,1),dims=2)
    println(size(cnames))
    println(size(times))
    println(size(cnames_rep))
    println(size(times_rep))

    println(size(betas))
    results = DataFrame(term=cnames_rep,estimate=betas,stderror=Missing,group="mass univariate",time=times_rep)
    return UnfoldModel(m,m.formula,tbl,results)

end

function condense(mm_array::Array{LinearMixedModel,1},tbl,times)
    # with random effects, no timeexpansion
    println(MixedModels.sqrtpwrss(mm_array[1]))
    results = condense_fixef.(mm_array,times)

    results = vcat(results,condense_ranef.(mm_array,times))
    # XXX TODO Return an Array of Models, not only the first one!
    return UnfoldModel(mm_array[1],mm_array[1].formula,tbl,vcat(results...))

end



function condense(mm::UnfoldLinearModel,tbl)
    # no random effects, timeexpansion

    if typeof(mm.formula)<: Array
        # build new UnfoldLinearModel for each formula

        times = [formula.rhs.basisfunction.times for formula in mm.formula]
        fromTo = [0 cumsum(length.(times),dims=2)]

        results = DataFrame()
        for k = 1:length(mm.formula)
            mm_tmp = UnfoldLinearModel(mm.beta[fromTo[k]+1:fromTo[k+1]],mm.optim,mm.formula[k],mm.X[:,fromTo[k]+1:fromTo[k+1]])
            res = condense_fixef(mm_tmp,times[k])
            results = vcat(results,res)
        end



    else
        times = mm.formula.rhs.basisfunction.times
        results = condense_fixef(mm,times)
    end

    return UnfoldModel(mm,mm.formula,tbl,results)
end


function condense(mm::LinearMixedModel,tbl)
    # with random effects, timeexpansion
    # TODO loop over al timeexpanded basisfunctions, in case there are multiple ones in case of multiple events
    # TODO random correlations => new function in MixedModels
    # TODO somehow get rid of the split of the coefficient names. What if there is a ":" in a coefficient name?
    times = mm.formula.rhs[1].basisfunction.times
    results = condense_fixef(mm,times)

    # ranefs more complex
    results = vcat(results,condense_ranef(mm,times))

    #return
    return UnfoldModel(mm,mm.formula,tbl,results)
end


function condense_fixef(mm,times)
    if typeof(mm.formula.rhs) <: Tuple
        #println("I am an array")
        fixefPart = mm.formula.rhs[1]
    else
        fixefPart = mm.formula.rhs
    end
    if typeof(times) <: Number
        times = [times]
    end
    cnames = [c[1] for c in split.(coefnames(fixefPart)," :")]

    times =  repeat(times,length(unique(cnames)))
    println("condense_fixef")
    println(size(cnames))
    println(size(stderror(mm)))
    println(size(MixedModels.fixef(mm)))

    return DataFrame(term=cnames,estimate=MixedModels.fixef(mm),stderror=MixedModels.stderror(mm),group="fixed",time=times)
end

function MixedModels.fixef(m::UnfoldLinearModel)
    # condense helper for the linear model which just returns a vector of fixef
     return m.beta
end
function MixedModels.stderror(m::UnfoldLinearModel)
    # for now we don't have an efficient way to calculate SEs for single subjects
    return fill(NaN,size(m.beta))
end
function condense_ranef(mm,times)
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
    if typeof(times) <: Number
        times = [times]
    end
    times =  repeat(times,length(unique(cnames)))
    # combine
    return DataFrame(term=cnames,estimate=σvec,stderror=NaN,group=group,time=times)

end
