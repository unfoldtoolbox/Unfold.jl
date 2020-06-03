
function extract_term_info(terms,ix)
    # 1 = basisname, 2 = coefname, 3 = colname
        return [c[ix] for c in split.(terms," :")]
end

function get_nUniqueterms(m::UnfoldLinearModel)
    return length(unique(extract_term_info(get_terms(m),[1,2])))
end
function get_nUniqueterms(m::UnfoldLinearMixedModel)
    return sum(length.(unique.(extract_term_info.(coefnames.(m.X.formulas.rhs),Ref([1,2])))))
end
function get_terms(m)
    terms = coefnames(m.X.formulas)
    terms = vcat(terms...) # gets rid of an empty Any() XXX not sure where it comes from, only in MixedModel Timexpanded case
end


function condense_long(m)
    terms = get_terms(m)
    terms = extract_term_info(terms,2)
    colnames_basis = get_colnames_basis(m.X.formulas)#get_colnames_basis(m.X.formulas,times)
    @debug terms
    @debug colnames_basis

    nchan = size(m.beta,1)
    nUniqueterms = get_nUniqueterms(m)
    @debug nUniqueterms
    terms_rep = permutedims(repeat(terms,1,nchan),[2,1])

    if length(colnames_basis)==1
        colnames_basis = [colnames_basis]
    end
    colnames_basis_rep = permutedims(repeat(colnames_basis,nUniqueterms,nchan),[2,1])

    chan_rep = repeat(1:nchan,1,size(colnames_basis_rep,2))

    return make_long_df(m,terms_rep,chan_rep,colnames_basis_rep) #DataFrame(term=linearize(terms_rep),estimate=linearize(m.beta),stderror=linearize(MixedModels.stderror(m)),channel=linearize(chan_rep),group="fixed",colnames_basis=linearize(colnames_rep))
end


function condense_long(m,times::AbstractArray)
    # Mass Univariate Case
    terms = coefnames(m.X.formulas)
    terms = vcat(terms...)
    colnames_basis = times#get_colnames_basis(m.X.formulas,times)

    @debug "terms: $(size(terms)),colnames_basis:$(size(colnames_basis)))"
    @debug "terms: $terms,colnames_basis:$colnames_basis"
    nchan = size(m.beta,1)
    ncols = length(colnames_basis)
    nterms = length(terms)

    # replicate
    terms_rep = permutedims(repeat(terms,outer=[1,nchan,ncols]),[2,3,1])
    colnames_basis_rep = permutedims(repeat(colnames_basis,1,nchan,nterms),[2 1 3])
    chan_rep = repeat(1:nchan,1,ncols,nterms)

    #
    results = make_long_df(m,terms_rep,chan_rep,colnames_basis_rep)

    return results
end
#---
# Returns a long df given the already matched
function make_long_df(m,terms,chans,colnames)
    @assert all(size(terms) .== size(chans)) "terms, chans and colnames need to have the same size at this point, $(size(terms)),$(size(chans)),$(size(colnames)), should be $(size(m.beta)))"
    @assert all(size(terms) .== size(colnames)) "terms, chans and colnames need to have the same size at this point"

    estimate,group = make_estimate(m)
    results = DataFrame(term=linearize(terms),
        channel = linearize(chans),
        colnames_basis=linearize(colnames),
        estimate=linearize(estimate),
        stderror=Missing,
        group = linearize(group),
        )
end
#---------

# extracts betas (and sigma's for mixed models) with string grouping indicator
# returns as a ch x beta, or ch x time x beta (for mass univariate)
function make_estimate(m::UnfoldLinearMixedModel)
    estimate = cat(m.beta,m.sigma,dims=ndims(m.beta))
    group_f = repeat(["fixef"],1,size(m.beta,ndims(m.beta)))
    group_s = repeat(["ranef"],1,size(m.sigma,ndims(m.beta)))

    group = hcat(repeat(group_f,size(m.beta,1)),repeat(group_s,size(m.beta,1)))
    return estimate,group
end
function make_estimate(m::UnfoldLinearModel)
    return m.beta,"mass univariate"
end

# Return the column names of the basis functions.
function get_colnames_basis(formula::FormulaTerm)
    return get_colnames_basis(formula.rhs)
end

function get_colnames_basis(rhs::Tuple)
    return rhs[1].basisfunction.colnames
end

function get_colnames_basis(rhs::AbstractTerm)
    return rhs.basisfunction.colnames
end

function get_basis_name(rhs::AbstractTerm)
    return rhs.basisfunction.name
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
