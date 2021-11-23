import Effects: effects
import Effects:_reference_grid
import Effects:typify
import Effects.typify
import Effects:_symequal


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
Effects.typify(reference_grid,form::Matrix,X;kwargs...) = typify.(Ref(reference_grid),form,Ref(X);kwargs...)

function effects(design::AbstractDict, model::UnfoldModel;typical=mean)
    reference_grid = _reference_grid(design)
    form = Unfold.formula(model) # get formula
    form_typical = typify(reference_grid, form, modelmatrix(model); typical=typical) # replace non-specified fields with "constants"
    #@show reference_grid
    @show size(form_typical)
    
    eff = yhat(model,form_typical,reference_grid)
    #eff = yhat(model,X) # apply coefficients
    @show eff
    @show reference_grid
    # because coefficients are 2D/3D arry, we have to cast it correctly to one big dataframe
    if isa(eff,Tuple)
        result = DataFrame(cast_referenceGrid(reference_grid,eff[3],eff[2] ))
        else
         result = DataFrame(cast_referenceGrid(reference_grid,eff,model ))
    end
    
return result   
end
cast_reference_grid(r,eff::Tuple{Any},model) = cast_reference_grid(r,eff[3],eff[1])


 cast_referenceGrid(r,eff,model::UnfoldLinearModel) = cast_referenceGrid(r,eff,times(model)[1])

 function cast_referenceGrid(r,eff,times)
   
    nchan = size(eff, 2)
    neff = size(r,1)
    neffCol = size(r,2)
    ncols = size(eff,1) ÷ neff # typically ntimes
   
    # replicate
    # for each predictor in r
    coefs_rep = Array{Float64}(undef,nchan,ncols,neff,neffCol)
    for k = 1:neffCol
        # repeat it for nchan + ncols
        coefs_rep[:,:,:,k] = permutedims(repeat(r[:,k], outer = [1, nchan, ncols]), [2, 3, 1])
    end
    @show times
    @show neff
    if length(times)  == neff*ncols
        # in case we have timeexpanded, times is already in long format and doesnt need to be repeated for each coefficient
        colnames_basis_rep = permutedims(repeat(times, 1, nchan, 1), [2 1 3])
        else
    colnames_basis_rep = permutedims(repeat(times, 1, nchan, neff), [2 1 3])
        end
    chan_rep = repeat(1:nchan, 1, ncols, neff)
    result = Dict(   :yhat => linearize(eff),
            :time => linearize(colnames_basis_rep),
            :channel => linearize(chan_rep))
   
    for k = 1:neffCol
            push!(result,Symbol(names(r)[k]) =>linearize(coefs_rep[:,:,:,k]))
    end
    return result
end


_symequal(t1::AbstractTerm,t2::Unfold.TimeExpandedTerm) = _symequal(t1,t2.term)