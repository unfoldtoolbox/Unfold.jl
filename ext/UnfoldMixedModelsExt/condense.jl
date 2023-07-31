
# extracts betas (and sigma's for mixed models) with string grouping indicator
# returns as a ch x beta, or ch x time x beta (for mass univariate)
function Unfold.make_estimate(
    m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    estimate = cat(coef(m), ranef(m), dims = ndims(coef(m)))
    if ndims(coef(m)) == 3
        group_f = repeat(
            [nothing],
            size(coef(m), 1),
            size(coef(m), 2),
            size(coef(m), ndims(coef(m))),
        )

        ranef_group = [x.group for x in MixedModels.tidyσs(modelfit(m))]
        
        # reshape to pred x time x chan and then invert to chan x time x pred
        ranef_group = permutedims(
            reshape(ranef_group, :, size(coef(m), 2), size(coef(m), 1)),
            [3 2 1],
        )


        group_s = ranef_group
        stderror_fixef = Unfold.stderror(m)
        stderror_ranef = fill(nothing, size(ranef(m)))
        stderror = cat(stderror_fixef, stderror_ranef, dims = 3)
    else
        group_f = repeat([nothing], size(coef(m), 1), size(coef(m), 2))
        group_s = repeat(["ranef"], size(coef(m), 1), size(ranef(m), 2))
        stderror = fill(nothing, size(estimate))
    end
    group = cat(group_f, group_s, dims = ndims(coef(m)))
    return Float64.(estimate), stderror, group
end

function Unfold.stderror(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})
    return permutedims(
        reshape(vcat([[b.se...] for b in modelfit(m).fits]...), reverse(size(coef(m)))),
        [3, 2, 1],
    )
end


function Unfold.get_coefnames(uf::UnfoldLinearMixedModelContinuousTime)
    # special case here, because we have to reorder the random effects to the end, else labels get messed up as we concat (coefs,ranefs)
 #   coefnames = Unfold.coefnames(formula(uf))
#    coefnames(formula(uf)[1].rhs[1])
    formulas = formula(uf)
    if !isa(formulas,AbstractArray) # in case we have only a single basisfunction
        formulas = [formulas]
    end
    fe_coefnames = vcat([coefnames(f.rhs[1]) for f in formulas]...)
    re_coefnames = vcat([coefnames(f.rhs[2:end]) for f in formulas]...)
    return vcat(fe_coefnames,re_coefnames)
end