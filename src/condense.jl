

function extract_coef_info(coefs, ix)
    # 1 = basisname, 2 = coefname, 3 = colname
    return [c[ix] for c in split.(coefs, " : ")]
end

function get_coefnames(uf::Union{UnfoldModel,DesignMatrix})
    coefnames = Unfold.coefnames(formula(uf))
    coefnames = vcat(coefnames...) # gets rid of an empty Any() XXX not sure where it comes from, only in MixedModel Timexpanded case
end

function get_coefnames(uf::UnfoldLinearMixedModelContinuousTime)
    # special case here, because we have to reorder the random effects to the end, else labels get messed up as we concat (coefs,ranefs)
 #   coefnames = Unfold.coefnames(formula(uf))
#    coefnames(formula(uf)[1].rhs[1])
    fe_coefnames = vcat([coefnames(f.rhs[1]) for f in formula(uf)]...)
    re_coefnames = vcat([coefnames(f.rhs[2:end]) for f in formula(uf)]...)
    return vcat(fe_coefnames,re_coefnames)
end
    
modelfit(uf::UnfoldModel) = uf.modelfit
StatsModels.coef(uf::UnfoldModel) = coef(modelfit(uf))
StatsModels.coef(mf::LinearModelFit) = mf.estimate


function StatsModels.coeftable(
    uf::Union{UnfoldLinearModelContinuousTime,UnfoldLinearMixedModelContinuousTime},
)
    coefsRaw = get_coefnames(uf)
    coefs = extract_coef_info(coefsRaw, 2)
    #colnames_basis_raw = get_colnames_basis(formula(uf))# this is unconverted basisfunction basis,
    colnames_basis = extract_coef_info(coefsRaw, 3) # this is converted to strings! 
    basisnames = extract_coef_info(coefsRaw, 1)
    @debug coefs
    @debug colnames_basis

    nchan = size(coef(uf), 1)

    coefs_rep = permutedims(repeat(coefs, 1, nchan), [2, 1])

    if length(colnames_basis) == 1
        colnames_basis = [colnames_basis]
    end


    #colnames_basis_rep = permutedims(repeat(colnames_basis_raw,Int(length(colnames_basis)/length(colnames_basis_raw)),nchan),[2,1])
    # XXX HOTFIX for above lines
    colnames_basis_rep = permutedims(repeat(colnames_basis, 1, nchan), [2, 1])
    try
        colnames_basis_rep = parse.(Float64, colnames_basis_rep)
    catch
        #do nothing
    end
    chan_rep = repeat(1:nchan, 1, size(colnames_basis_rep, 2))

    basisnames_rep = permutedims(repeat(basisnames, 1, nchan), [2, 1])


    return make_long_df(
        uf,
        coefs_rep,
        chan_rep,
        colnames_basis_rep,
        basisnames_rep,
        collabel(uf),
    )
end


function StatsModels.coeftable(uf::Union{UnfoldLinearModel,UnfoldLinearMixedModel})
    # Mass Univariate Case
    coefnames = get_coefnames(uf)

    colnames_basis = collect(values(design(uf)))[1][2]#times#get_colnames_basis(m.X.formulas,times)

    @debug "coefs: $(size(coefnames)),colnames_basis:$(size(colnames_basis)))"
    @debug "coefs: $coefnames,colnames_basis:$colnames_basis"
    nchan = size(coef(uf), 1)
    ncols = length(colnames_basis)
    ncoefs = length(coefnames)

    # replicate
    coefs_rep = permutedims(repeat(coefnames, outer = [1, nchan, ncols]), [2, 3, 1])
    colnames_basis_rep = permutedims(repeat(colnames_basis, 1, nchan, ncoefs), [2 1 3])
    chan_rep = repeat(1:nchan, 1, ncols, ncoefs)
    basisnames_rep = repeat(["mass-univariate"], nchan, ncols, ncoefs)
    #
    results =
        make_long_df(uf, coefs_rep, chan_rep, colnames_basis_rep, basisnames_rep, :time)

    return results
end


#---
# Returns a long df given the already matched
function make_long_df(m, coefs, chans, colnames, basisnames, collabel)
    @assert all(size(coefs) .== size(chans)) "coefs, chans and colnames need to have the same size at this point, $(size(coefs)),$(size(chans)),$(size(colnames)), should be $(size(coef(m))))"
    @assert all(size(coefs) .== size(colnames)) "coefs, chans and colnames need to have the same size at this point"
    estimate, stderror, group = make_estimate(m)


    return DataFrame(
        Dict(
            :coefname => String.(linearize(coefs)),
            :channel => linearize(chans),
            :basisname => linearize(basisnames),
            collabel => linearize(colnames),
            :estimate => linearize(estimate),
            :stderror => linearize(stderror),
            :group => linearize(group),
        ),
    )
end
#---------

# extracts betas (and sigma's for mixed models) with string grouping indicator
# returns as a ch x beta, or ch x time x beta (for mass univariate)
function make_estimate(
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

stderror(m::UnfoldModel) = stderror(modelfit(m))
function stderror(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})
    return permutedims(
        reshape(vcat([[b.se...] for b in modelfit(m).fits]...), reverse(size(coef(m)))),
        [3, 2, 1],
    )
end
function stderror(m::LinearModelFit)
    if isempty(m.standarderror)
        stderror = fill(nothing, size(coef(m)))
    else
        stderror = Float64.(m.standarderror)
    end
end

function make_estimate(uf::UnfoldModel)
    return Float64.(coef(uf)), stderror(uf), fill(nothing, size(coef(uf)))
end

# Return the column names of the basis functions.
function get_colnames_basis(formula::FormulaTerm)
    return get_colnames_basis(formula.rhs)
end

function get_colnames_basis(rhs::Tuple)
    return colnames(rhs[1].basisfunction)
end

function get_colnames_basis(rhs::AbstractTerm)
    return colnames(rhs.basisfunction)
end

function get_basis_name(m::UnfoldModel)
    return extract_coef_info(Unfold.get_coefnames(m), 1)
end
function get_basis_name(rhs::AbstractTerm)
    return name(rhs.basisfunction)
end

function get_colnames_basis(formulas::AbstractArray{<:FormulaTerm})
    # In case of multiple basisfunction we can have an array of formulas.
    # in that case we have to add an unique identifier
    colnames_all = []
    for formula in formulas
        append!(colnames_all, get_colnames_basis(formula.rhs))
    end
    #get_colnames_basis_name(formula.rhs) .*
    return colnames_all
end
