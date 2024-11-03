"""
$(SIGNATURES)

Returns a partial LMM model (non-functional due to lacking data) to be used in likelihoodratiotests.
`k` to selcet which of the modelfit's to fake
"""
function fake_lmm(data::AbstractArray{<:Number,3}, m::UnfoldLinearMixedModel, k::Int)
    mm = modelmatrix(m)
    @assert length(mm) == 1 "LRT is currently not implemented for fitting multiple events at the same time"
    fcoll = Unfold.modelfit(m)
    feterm = mm[1][1]
    #reterm = mm[2:end]
    #fakeY = zeros(size(feterm, 1))
    channel = fcoll.fits[k].channel
    time_ix = fcoll.fits[k].timeIX
    y = data[channel, time_ix, :]
    lmm = LinearMixedModel_wrapper(Unfold.formulas(m), y, mm[1])

    lmm.optsum.feval = 1
    lmm.optsum.fmin = 1
    lmm.optsum.sigma = fcoll.fits[k].σ
    lmm.optsum.optimizer = :unfold
    lmm.optsum.returnvalue = :success
    setθ!(lmm, fcoll.fits[k].θ)
    updateL!(lmm)
    #MixedModels.setβ!(lmm, fcoll.fits[k].β)
    return lmm
end

fake_lmm(data::AbstractArray, m::UnfoldLinearMixedModel) =
    fake_lmm.(data, m, 1:length(modelfit(m).fits))

"""
$(SIGNATURES)

Calculate likelihoodratiotest
"""
function MixedModels.likelihoodratiotest(data::AbstractArray, m::UnfoldLinearMixedModel...)
    #@info lrtest(fake_lmm.(m,1)...)
    n = length(Unfold.modelfit(m[1]).fits)
    lrt = Array{MixedModels.LikelihoodRatioTest}(undef, n)
    sizehint!(lrt, n)
    for k = 1:n

        ms = fake_lmm.(Ref(data), m, k)
        #@info objective.(ms)
        lrt[k] = MixedModels.likelihoodratiotest(ms...)
    end
    return lrt
end

"""
$(SIGNATURES)
Unfold-Method: return pvalues of likelihoodratiotests, typically calculated:

# Examples
julia> pvalues(likelihoodratiotest(m1,m2))

where m1/m2 are UnfoldLinearMixedModel's

Tipp: if you only compare two models you can easily get a vector of p-values:

julia> vcat(pvalues(likelihoodratiotest(m1,m2))...)


Multiple channels are returned linearized at the moment, as we do not have access to the amount of channels after the LRT, you can do:

julia> reshape(vcat(pvalues(likelihoodratiotest(m1,m2))...),ntimes,nchan)'

"""
function pvalues(lrtvec::Vector{MixedModels.LikelihoodRatioTest})
    [lrt.pvalues for lrt in lrtvec]
end



function MixedModels.rePCA(m::UnfoldLinearMixedModel)
    oneresult = MixedModels.rePCA(fake_lmm(m, 1))
    emptyPCA = (;)
    for s in keys(oneresult)
        res = hcat(
            [
                a[s] for
                a in MixedModels.rePCA.(fake_lmm.(Ref(m), 1:length(modelfit(m).fits)))
            ]...,
        )
        newPCA = NamedTuple{(s,)}([res])
        emptyPCA = merge(emptyPCA, newPCA)
    end
    return emptyPCA
end
