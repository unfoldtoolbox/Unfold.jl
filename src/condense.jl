
function extract_term_info(terms,ix)
    # 1 = basisname, 2 = coefname, 3 = colname
        return [c[ix] for c in split.(terms," : ")]
end

function get_terms(uf)
    terms = coefnames(formula(uf))
    terms = vcat(terms...) # gets rid of an empty Any() XXX not sure where it comes from, only in MixedModel Timexpanded case
end

modelfit(uf::UnfoldModel) = uf.modelfit
StatsModels.coef(uf::UnfoldModel) = coef(modelfit(uf))
StatsModels.coef(mf::LinearModelFit) = mf.estimate


function StatsModels.coeftable(uf::UnfoldModel)
    termsRaw = get_terms(uf)
    terms = extract_term_info(termsRaw,2)
    #colnames_basis_raw = get_colnames_basis(formula(uf))# this is unconverted basisfunction basis,
    colnames_basis = extract_term_info(termsRaw,3) # this is converted to strings! 
    basisnames = extract_term_info(termsRaw,1)
    @debug terms
    @debug colnames_basis

    nchan = size(coef(uf),1)

    terms_rep = permutedims(repeat(terms,1,nchan),[2,1])

    if length(colnames_basis)==1
        colnames_basis = [colnames_basis]
    end
    
    
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
    return make_long_df(uf,terms_rep,chan_rep,colnames_basis_rep,basisnames_rep) #DataFrame(term=linearize(terms_rep),estimate=linearize(coef(m)),stderror=linearize(MixedModels.stderror(m)),channel=linearize(chan_rep),group="fixed",colnames_basis=linearize(colnames_rep))
end


function StatsModels.coeftable(uf::Union{UnfoldLinearModel,UnfoldLinearMixedModel})
    # Mass Univariate Case
    terms = coefnames(formula(uf))
    terms = vcat(terms...)
    
    colnames_basis = collect(values(design(uf)))[1][2]#times#get_colnames_basis(m.X.formulas,times)

    @debug "terms: $(size(terms)),colnames_basis:$(size(colnames_basis)))"
    @debug "terms: $terms,colnames_basis:$colnames_basis"
    nchan = size(coef(uf),1)
    ncols = length(colnames_basis)
    nterms = length(terms)

    # replicate
    terms_rep = permutedims(repeat(terms,outer=[1,nchan,ncols]),[2,3,1])
    colnames_basis_rep = permutedims(repeat(colnames_basis,1,nchan,nterms),[2 1 3])
    chan_rep = repeat(1:nchan,1,ncols,nterms)
    basisnames_rep = repeat(["mass-univariate"],nchan,ncols,nterms)
    #
    results = make_long_df(uf,terms_rep,chan_rep,colnames_basis_rep,basisnames_rep)

    return results
end
#---
# Returns a long df given the already matched
function make_long_df(m,terms,chans,colnames,basisnames)
    @assert all(size(terms) .== size(chans)) "terms, chans and colnames need to have the same size at this point, $(size(terms)),$(size(chans)),$(size(colnames)), should be $(size(coef(m))))"
    @assert all(size(terms) .== size(colnames)) "terms, chans and colnames need to have the same size at this point"
    estimate,stderror,group = make_estimate(m)
    
    
    return  DataFrame(term=String.(linearize(terms)),
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
function make_estimate(m::Union{UnfoldLinearMixedModel,UnfoldLinearMixedModelContinuousTime})
    estimate = cat(coef(m),ranef(m),dims=ndims(coef(m)))
    if ndims(coef(m)) == 3
        group_f = repeat([nothing],size(coef(m),1),size(coef(m),2),size(coef(m),ndims(coef(m))))

        ranef_group = [x.group for x in MixedModels.tidyσs(modelfit(m))]

            # reshape to pred x time x chan and then invert to chan x time x pred
        ranef_group = permutedims(reshape(ranef_group,:,size(coef(m),2),size(coef(m),1)),[3 2 1])
            
       
        group_s = ranef_group
        stderror_fixef = Unfold.stderror(m)
        stderror_ranef = fill(nothing,size(ranef(m)))
        stderror = cat(stderror_fixef,stderror_ranef,dims=3)
    else
        group_f = repeat([nothing],size(coef(m),1),size(coef(m),ndims(coef(m))))
        group_s = repeat(["ranef"],size(coef(m),1),size(ranef(m),ndims(coef(m))))
        stderror = fill(nothing,size(estimate))
    end
    group = cat(group_f,group_s,dims=ndims(coef(m)))
    return Float64.(estimate),stderror,group
end

stderror(m::UnfoldModel) = stderror(modelfit(m))
function stderror(m::Union{UnfoldLinearMixedModel, UnfoldLinearMixedModelContinuousTime})
    return permutedims(reshape(vcat([[b.se...] for b in modelfit(m).fits]...),reverse(size(coef(m)))),[3,2,1])
end
function stderror(m::LinearModelFit)
    if isempty(m.standarderror)
        stderror = fill(nothing,size(coef(m)))
    else
        stderror = Float64.(m.standarderror)
    end
end

function make_estimate(uf::UnfoldModel)
    return Float64.(coef(uf)),stderror(uf),fill(nothing,size(coef(uf)))
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

function get_basis_name(m::UnfoldModel)
    return extract_term_info(Unfold.get_terms(m),1)
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
