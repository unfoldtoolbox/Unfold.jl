

function extract_coef_info(coefs, ix)
    # 1 = eventname, 2 = coefname, 3 = colname
    return [c[ix] for c in split.(coefs, " : ")]
end

function get_coefnames(uf::Union{<:UnfoldModel,<:AbstractDesignMatrix})
    coefnames = Unfold.coefnames(formulas(uf))
    coefnames = vcat(coefnames...) # gets rid of an empty Any() XXX not sure where it comes from, only in MixedModel Timexpanded case
end



modelfit(uf::UnfoldModel) = uf.modelfit

StatsModels.coef(uf::UnfoldModel) = coef(modelfit(uf))
StatsModels.coef(mf::LinearModelFit) = mf.estimate





@traitfn function StatsModels.coeftable(
    uf::T,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}}
    coefsRaw = get_coefnames(uf)
    coefs = extract_coef_info(coefsRaw, 2)
    #colnames_basis_raw = get_colnames_basis(formulas(uf))# this is unconverted basisfunction basis,
    colnames_basis = extract_coef_info(coefsRaw, 3) # this is converted to strings! 
    #eventnames = extract_coef_info(coefsRaw, 1)
    @debug coefs
    @debug colnames_basis

    nchan = size(coef(uf), 1)

    coefs_rep = permutedims(repeat(coefs, 1, nchan), [2, 1])

    if length(colnames_basis) == 1
        colnames_basis = [colnames_basis]
    end


    colnames_basis_rep = permutedims(repeat(colnames_basis, 1, nchan), [2, 1])
    try
        colnames_basis_rep = parse.(Float64, colnames_basis_rep)
    catch
        #do nothing
    end
    chan_rep = repeat(1:nchan, 1, size(colnames_basis_rep, 2))



    designkeys = collect(first.(design(uf)))
    if length(designkeys) == 1
        # in case of 1 event, repeat it by ncoefs
        eventnames = repeat([designkeys[1]], length(chan_rep))
    else
        eventnames = String[]
        sizehint!(eventnames, length(chan_rep))
        for (ix, evt) in enumerate(designkeys)
            push!(
                eventnames,
                repeat([evt], size(modelmatrices(designmatrix(uf))[ix], 2))...,
            )
        end
    end
    eventnames_rep = permutedims(repeat(eventnames, 1, nchan), [2, 1])


    return make_long_df(
        uf,
        coefs_rep,
        chan_rep,
        colnames_basis_rep,
        eventnames_rep,
        collabel(uf),
    )
end


@traitfn function StatsModels.coeftable(
    uf::T,
) where {T <: UnfoldModel; !ContinuousTimeTrait{T}}
    # Mass Univariate Case
    coefnames = get_coefnames(uf)

    colnames_basis = collect(last.(last.(design(uf)[1])))

    @debug "coefs: $(size(coefnames)),colnames_basis:$(size(colnames_basis)))"
    @debug "coefs: $coefnames,colnames_basis:$colnames_basis"
    nchan = size(coef(uf), 1)
    ncols = length(colnames_basis)
    ncoefs = length(coefnames)

    # replicate
    coefs_rep = permutedims(repeat(coefnames, outer = [1, nchan, ncols]), [2, 3, 1])
    colnames_basis_rep = permutedims(repeat(colnames_basis, 1, nchan, ncoefs), [2 1 3])
    chan_rep = repeat(1:nchan, 1, ncols, ncoefs)

    designkeys = collect(first.(design(uf)))
    if length(designkeys) == 1
        # in case of 1 event, repeat it by ncoefs
        eventnames = repeat([designkeys[1]], ncoefs)
    else
        eventnames = String[]
        for (ix, evt) in enumerate(designkeys)
            push!(eventnames, repeat([evt], size(modelmatrix(uf)[ix], 2))...)
        end
    end

    eventnames_rep = permutedims(repeat(eventnames, 1, nchan, ncols), [2, 3, 1])
    #
    results =
        make_long_df(uf, coefs_rep, chan_rep, colnames_basis_rep, eventnames_rep, :time)

    return results
end


#---
# Returns a long df given the already matched
function make_long_df(m, coefs, chans, colnames, eventnames, collabel)
    @assert all(size(coefs) .== size(chans)) "coefs, chans and colnames need to have the same size at this point, $(size(coefs)),$(size(chans)),$(size(colnames)), should be $(size(coef(m))))"
    @assert all(size(coefs) .== size(colnames)) "coefs, chans and colnames need to have the same size at this point"
    estimate, stderror, group = make_estimate(m)


    return DataFrame(
        Dict(
            :coefname => String.(linearize(coefs)),
            :channel => linearize(chans),
            :eventname => linearize(eventnames),
            collabel => linearize(colnames),
            :estimate => linearize(estimate),
            :stderror => linearize(stderror),
            :group => linearize(group),
        ),
    )
end
#---------

stderror(m::UnfoldModel) = stderror(modelfit(m))

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