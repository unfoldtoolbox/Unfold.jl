import Effects: effects
using Effects

# Implementation based on Effects Package; likely should repackage in UnfoldEffects; but then...
function effects(design::AbstractDict, model::UnfoldModel;
    typical=mean)
    reference_grid = Effects._reference_grid(design)
    form = formula(model) # get formula
    form_typical = Effects.typify(reference_grid, form, modelmatrix(model); typical=typical) # replace non-specified fields with "constants"
    X = modelcols(form_typical, reference_grid) # get model cols
    eff = yhat(model,X) # apply coefficients
    
    # because coefficients are 2D/3D arry, we have to cast it correctly to one big dataframe
    result = DataFrame(cast_referenceGrid(reference_grid,eff,times(model)[1]) )
    
return result   
end

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
    colnames_basis_rep = permutedims(repeat(times, 1, nchan, neff), [2 1 3])
    chan_rep = repeat(1:nchan, 1, ncols, neff)
    result = Dict(   :yhat => linearize(eff),
            :time => linearize(colnames_basis_rep),
            :channel => linearize(chan_rep))
   
    for k = 1:neffCol
            push!(result,Symbol(names(r)[k]) =>linearize(coefs_rep[:,:,:,k]))
    end
    return result
end