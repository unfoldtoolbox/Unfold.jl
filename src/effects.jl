import Effects: effects
import Effects:_reference_grid
import Effects:typify
import Effects.typify
import Effects:_symequal
import StatsModels.collect_matrix_terms
import Base.getproperty
using Effects

"""
effects(design::AbstractDict, model::UnfoldModel;typical=mean)

Calculates marginal effects for all term-combinations in `design`.

 Implementation based on Effects Package; likely could repackage in UnfoldEffects; somebody wants to do it? This would make it easier to cross-maintain it to changes/bugfixes in the Effects.jl Package
 `design` is a Dictionary containing those predictors (als keys) with levels (as values), that you want to evaluate. The `typical` refers to the value, that other predictors not in the Dictionary should take on.

 # Example
 ```julia-repl
 julia> f = @formula 0 ~ categoricalA + continuousA + continuousB
 julia> uf = fit(UnfoldModel,(Any=>(f,times)),data,events)
 julia> d = Dict(:categoricalA=>["levelA","levelB"],:continuousB=>[-2,0,2])
 julia> effects(d,uf)
```
 will result in 6 predicted values: A/-2, A/0, A/2, B/-2, B/0, B/2.
""" 


function effects(design::AbstractDict, model::UnfoldModel;typical=mean)
    reference_grid = _reference_grid(design)
    form = Unfold.formula(model) # get formula

    # replace non-specified fields with "constants"
    m = modelmatrix(model,false) # get the modelmatrix without timeexpansion
    
    if isa(form,AbstractMatrix)

        form_typical = Array{Any}(undef,1, length(form))
        for f = 1:length(form)
            
            # strip of basisfunction and put it on afterwards again
            tmpf = deepcopy(form[f])
            
            # create a Formula without Timeexpansion
            tmpf = FormulaTerm(tmpf.lhs,tmpf.rhs.term)

            # typify that
            tmpf = typify(reference_grid, tmpf, m[f]; typical=typical) 
            
            # regenerate TimeExpansion
            tmpf = Unfold.TimeExpandedTerm(tmpf,form[f].rhs.basisfunction;eventfields=form[f].rhs.eventfields)
            form_typical[f] = tmpf
        end
    else
        form_typical = [typify(reference_grid, form, m; typical=typical)]
    end
    eff = yhat(model,form_typical,reference_grid)

    # because coefficients are 2D/3D arry, we have to cast it correctly to one big dataframe
    if isa(eff,Tuple) 
        # TimeContinuous Model, we also get back other things like times & fromWhereToWhere a BasisFunction goes
        if typeof(form) <: FormulaTerm
            form = [form]
        end
        # figure out the basisname
        bnames = [[form[ix].rhs.basisfunction.name] for ix in 1:length(form)]
        # repeat it as much as necessary
        bnames = repeat.(bnames,[e.stop+e.step-1 for e in eff[1]])

        result = DataFrame(cast_referenceGrid(reference_grid,eff[3],eff[2] ;basisname=vcat(bnames...)))
        select!(result,Not(:latency)) # remove the latency column if it was added
    else
        # normal mass univariate model
        result = DataFrame(cast_referenceGrid(reference_grid,eff,times(model)[1] ))
    end
    
return result   
end
 
 Effects.typify(reference_grid,form::Matrix,X;kwargs...) = typify.(Ref(reference_grid),form,Ref(X);kwargs...)

 
function cast_referenceGrid(r,eff,times;basisname=nothing)
    nchan = size(eff, 2) # correct
    neff = size(r,1) # how many effects requested
    neffCol = size(r,2) # how many predictors
    ncols = size(eff,1) ÷ neff # typically ntimes

    
   
    # replicate
    # for each predictor in r (reference grid), we need this at the bottom

    if isnothing(basisname)
        nbases = 1
    else
        nbases = length(unique(basisname))
    end
    coefs_rep = Array{Array}(undef,nbases,neffCol)
    
    
    for k = 1:neffCol
        # in case we have only a single basis (e.g. mass univariate), we can directly fill in all values
        ixList = []
        if isnothing(basisname)
            ix = ones(ncols) .== 1.
            append!(ixList,[ix])
        else
            #in case of multiple bases, we have to do it iteratively, because the bases can be different length
            for b = unique(basisname)
                ix = basisname[1:neff:end].==b
                append!(ixList,[ix])
            end
        end
        for i_ix = 1:length(ixList)
            
            coefs_rep[i_ix,k] = linearize(permutedims(repeat(r[:,k], outer = [1, nchan, sum(ixList[i_ix])]), [2,3,1]))
        end
    end
    
    # often the "times" vector
    if length(times)  == neff*ncols
        # in case we have timeexpanded, times is already in long format and doesnt need to be repeated for each coefficient
        colnames_basis_rep = permutedims(repeat(times, 1, nchan, 1), [2 1 3])
    else
        colnames_basis_rep = permutedims(repeat(times, 1, nchan, neff), [2 1 3])
    end

    # for multiple channels
    chan_rep = repeat(1:nchan, 1, ncols, neff)
    
    # for mass univariate there is no basisname
    if isnothing(basisname)
        basisname = fill(nothing,ncols)
        basisname_rep = permutedims(repeat(basisname,1,nchan,neff),[2,1,3])
    else
        basisname_rep = permutedims(repeat(basisname,1,nchan,1),[2,1,3])
    end
    

    result = Dict(   :yhat => linearize(eff'),
            :time => linearize(colnames_basis_rep),
            :channel => linearize(chan_rep),
            :basisname =>linearize(basisname_rep))

    for k = 1:neffCol
            push!(result,Symbol(names(r)[k]) =>vcat(vcat(coefs_rep[:,k]...)...))
    end

    return result
end
Effects._trmequal(t1::Unfold.uf_bsplineTerm,t2::AbstractTerm) = _symequal(t1.term,t2)
Effects._trmequal(t1::Unfold.uf_bsplineTerm,t2::Unfold.uf_bsplineTerm) = _symequal(t1.term,t2.term)

Effects._trmequal(t1::AbstractTerm,t2::Unfold.uf_bsplineTerm) = _symequal(t1,t2.term)
Effects._symequal(t1::AbstractTerm,t2::Unfold.uf_bsplineTerm) = _symequal(t1,t2.term)
#Effects._symequal(t1::AbstractTerm,t2::Unfold.TimeExpandedTerm) = _symequal(t1,t2.term)
#function Effects._replace(matrix_term::MatrixTerm{<:Tuple{<:Unfold.TimeExpandedTerm}},typicals::Dict)
    
#    replaced_term = MatrixTerm((Effects._replace.(matrix_term.terms, Ref(typicals))...,))
#    basisfunctionTerm = getfield(matrix_term,:terms)[1]
#    return TimeExpandedTerm(replaced_term,basisfunctionTerm.basisfunction,basisfunctionTerm.eventfields)

#end
