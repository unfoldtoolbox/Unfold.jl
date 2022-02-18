"""
$(SIGNATURES)

Returns a partial LMM model (non-functional due to lacking data) to be used in likelihoodratiotests.
`k` to selcet which of the modelfit's to fake
"""
function fake_lmm(m::UnfoldLinearMixedModel,k::Int)
    (feterm,reterm) = Unfold.designmatrix(m).Xs
    fakeY = zeros(size(feterm,1))
    lmm = Unfold.LinearMixedModel_wrapper(Unfold.formula(m),fakeY,designmatrix(m).Xs)

    fcoll = Unfold.modelfit(m)
    #lmm.objective .= fcoll.fits[1].objective
    #lmm.optsum.feval .= 1
    #lmm.optsum.fmin .= 1
    lmm.optsum.sigma = fcoll.fits[k].σ
    lmm.optsum.optimizer = :unfold
    lmm.optsum.returnvalue = :success
    setθ!(lmm,fcoll.fits[k].θ)
    return lmm
end
    fake_lmm(m::UnfoldLinearMixedModel) = fake_lmm.(m,1:length(modelfit(m).fits))
	
    """
    $(SIGNATURES)

    Calculate likelihoodratiotest
    """
	function MixedModels.likelihoodratiotest(m::UnfoldLinearMixedModel...)
		#@info lrtest(fake_lmm.(m,1)...)
		n = length(Unfold.modelfit(m[1]).fits)
		lrt = Array{MixedModels.LikelihoodRatioTest}(undef,n)
		sizehint!(lrt,n)
		for k = 1:n
			ms = fake_lmm.(m,k)
			@info objective.(ms)
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

    """
	function pvalues(lrtvec::Vector{MixedModels.LikelihoodRatioTest})
		[lrt.pvalues for lrt in lrtvec]
	end
