import Effects.effects

function Effects.effects(design::AbstractDict, model::UnfoldModel;
    eff_col=nothing, err_col=:err, typical=mean,
    lower_col=:lower, upper_col=:upper)
grid = _reference_grid(design)
dv = something(eff_col, formula(model).lhs.sym)
effects!(grid, model; eff_col=dv, err_col=err_col, typical=typical)
# XXX DataFrames dependency
grid[!, lower_col] = grid[!, dv] - grid[!, err_col]
grid[!, upper_col] = grid[!, dv] + grid[!, err_col]
return grid
# up_low = let dv = getproperty(reference_grid, dv), err = getproperty(reference_grid, err_col)
#     (; lower_col => dv .- err, upper_col => dv .+ err)
# end
# return (; reference_grid..., up_low...)
end
