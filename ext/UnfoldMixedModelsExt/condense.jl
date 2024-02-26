
function MixedModels.tidyσs(m::UnfoldLinearMixedModel)
    #    using MixedModels: AbstractReTerm
    t = MixedModels.tidyσs(modelfit(m))
    reorder_tidyσs(t, Unfold.formula(m))
end

MixedModels.tidyβ(m::UnfoldLinearMixedModel) = MixedModels.tidyβ(modelfit(m))

"""
    randomeffectgroupings(t::MixedModels.AbstractReTerm)
Returns the random effect grouping term (rhs), similar to coefnames, which returns the left hand sides
"""
randomeffectgroupings(t::AbstractTerm) = repeat([nothing], length(t.terms))
randomeffectgroupings(t::MixedModels.AbstractReTerm) =
    repeat([t.rhs.sym], length(t.lhs.terms))

"""
    reorder_tidyσs(t, f)
This function reorders a MixedModels.tidyσs output, according to the formula and not according to the largest RandomGrouping.

"""
function reorder_tidyσs(t, f)

    # get the order from the formula, this is the target
    f_order = randomeffectgroupings.(f.rhs) # formula order
    f_order = vcat(f_order...)

    # find the fixefs
    fixef_ix = isnothing.(f_order)



    f_order = string.(f_order[.!fixef_ix])
    f_name = coefnames(f.rhs)[.!fixef_ix]

    # get order from tidy object
    t_order = [string(i.group) for i in t if i.iter == 1]
    t_name = [string(i.column) for i in t if i.iter == 1]

    # combine for formula and tidy output the group + the coefname
    f_comb = f_order .* f_name
    t_comb = t_order .* t_name

    # find for each formula output, the fitting tidy permutation
    reorder_ix = Int[]
    for f_coef in f_comb
        ix = findall(t_comb .== f_coef)
        @assert length(ix) == 1 "error in reordering of MixedModels - please file a bugreport!"
        push!(reorder_ix, ix[1])
    end
    @assert length(reorder_ix) == length(t_comb)

    # repeat and build the index for all timepoints
    reorder_ix_all = repeat(reorder_ix, length(t) ÷ length(reorder_ix))
    for k = 1:length(reorder_ix):length(t)
        reorder_ix_all[k:k+length(reorder_ix)-1] .+= (k - 1)
    end

    return t[reorder_ix_all]


end

"""
    Unfold.make_estimate(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
extracts betas (and sigma's for mixed models) with string grouping indicator

returns as a ch x beta, or ch x time x beta (for mass univariate)
"""
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

        ranef_group = [x.group for x in MixedModels.tidyσs(m)]

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

function Unfold.stderror(
    m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime},
)
    return permutedims(
        reshape(vcat([[b.se...] for b in modelfit(m).fits]...), reverse(size(coef(m)))),
        [3, 2, 1],
    )
end


function Unfold.get_coefnames(uf::UnfoldLinearMixedModelContinuousTime)
    # special case here, because we have to reorder the random effects to the end, else labels get messed up as we concat (coefs,ranefs)
    #   coefnames = Unfold.coefnames(formula(uf))
    #    coefnames(formula(uf)[1].rhs[1])
    formulas = Unfold.formula(uf)
    if !isa(formulas, AbstractArray) # in case we have only a single basisfunction
        formulas = [formulas]
    end
    fe_coefnames = vcat([coefnames(f.rhs[1]) for f in formulas]...)
    re_coefnames = vcat([coefnames(f.rhs[2:end]) for f in formulas]...)
    return vcat(fe_coefnames, re_coefnames)
end
