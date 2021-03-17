
function extract_term_info(terms,ix)
    # 1 = basisname, 2 = coefname, 3 = colname
        return [c[ix] for c in split.(terms," : ")]
end

function get_terms(m)
    terms = coefnames(m.X.formulas)
    terms = vcat(terms...) # gets rid of an empty Any() XXX not sure where it comes from, only in MixedModel Timexpanded case
end


function condense_long(m)
    termsRaw = get_terms(m)
    terms = extract_term_info(termsRaw,2)
    colnames_basis_raw = get_colnames_basis(m.X.formulas)# this is unconverted basisfunction basis,
    colnames_basis = extract_term_info(termsRaw,3) # this is converted to strings! 
    basisnames = extract_term_info(termsRaw,1)
    @debug terms
    @debug colnames_basis

    nchan = size(m.beta,1)

    terms_rep = permutedims(repeat(terms,1,nchan),[2,1])

    if length(colnames_basis)==1
        colnames_basis = [colnames_basis]
    end
    
    #println(colnames_basis_raw[1:100])
    #println(colnames_basis[1:100])
    #println(size(terms_rep))
    #println("$(length(colnames_basis_raw)),$(length(colnames_basis)),$nchan")
    #colnames_basis_rep = permutedims(repeat(colnames_basis_raw,Int(length(colnames_basis)/length(colnames_basis_raw)),nchan),[2,1])
    # XXX HOTFIX for above lines
    colnames_basis_rep = permutedims(repeat(colnames_basis,1,nchan),[2,1])
    try
        colnames_basis_rep = parse.(Float64,colnames_basis_rep)
    catch
        #do nothing
    end
    chan_rep = repeat(1:nchan,1,size(colnames_basis_rep,2))

    basisnames_rep = permutedims(repeat(basisnames,1,nchan),[2,1])
    return make_long_df(m,terms_rep,chan_rep,colnames_basis_rep,basisnames_rep) #DataFrame(term=linearize(terms_rep),estimate=linearize(m.beta),stderror=linearize(MixedModels.stderror(m)),channel=linearize(chan_rep),group="fixed",colnames_basis=linearize(colnames_rep))
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
    basisnames_rep = repeat(["mass-univariate"],nchan,ncols,nterms)
    #
    results = make_long_df(m,terms_rep,chan_rep,colnames_basis_rep,basisnames_rep)

    return results
end
#---
# Returns a long df given the already matched
function make_long_df(m,terms,chans,colnames,basisnames)
    @assert all(size(terms) .== size(chans)) "terms, chans and colnames need to have the same size at this point, $(size(terms)),$(size(chans)),$(size(colnames)), should be $(size(m.beta)))"
    @assert all(size(terms) .== size(colnames)) "terms, chans and colnames need to have the same size at this point"
    estimate,stderror,group = make_estimate(m)
    
    
    results = DataFrame(term=String.(linearize(terms)),
        channel = linearize(chans),
        basisname = linearize(basisnames),
        colname_basis=linearize(colnames),
        estimate=linearize(estimate),
        stderror=linearize(stderror),
        group = linearize(group),
        )
end
#---------

# extracts betas (and sigma's for mixed models) with string grouping indicator
# returns as a ch x beta, or ch x time x beta (for mass univariate)
function make_estimate(m::UnfoldLinearMixedModel)
    estimate = cat(m.beta,m.sigma,dims=ndims(m.beta))
    if ndims(m.beta) == 3
        group_f = repeat([nothing],size(m.beta,1),size(m.beta,2),size(m.beta,ndims(m.beta)))

        ranef_group = [x.group for x in MixedModels.tidyσs(m.modelinfo)]

            # reshape to pred x time x chan and then invert to chan x time x pred
        ranef_group = permutedims(reshape(ranef_group,:,size(m.beta,2),size(m.beta,1)),[3 2 1])
            
        
        #    beta = reshape(beta,:,nchan)'
        #    sigma = reshape(sigma,:,nchan)'
        
        #group_s = repeat(["ranef"],size(m.beta,1),size(m.beta,2),size(m.sigma,ndims(m.beta)))
        group_s = ranef_group
        stderror_fixef = permutedims(reshape(vcat([[b.se...] for b in m.modelinfo.bstr]...),reverse(size(m.beta))),[3,2,1])
        stderror_ranef = fill(nothing,size(m.sigma))
        stderror = cat(stderror_fixef,stderror_ranef,dims=3)
    else
        group_f = repeat([nothing],size(m.beta,1),size(m.beta,ndims(m.beta)))
        group_s = repeat(["ranef"],size(m.beta,1),size(m.sigma,ndims(m.beta)))
        stderror = fill(nothing,size(estimate))
    end
    group = cat(group_f,group_s,dims=ndims(m.beta))
    return Float64.(estimate),stderror,group
end
function make_estimate(m::UnfoldLinearModel)
    return Float64.(m.beta),fill(Missing,size(m.beta)),nothing
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
