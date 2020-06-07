# These functions are not included yet because I am not yet sure which plotting library to put my bets on.
# Makie&StatsMakie could be supercool, but I ran into some issues
import Plots
using SparseArrays
Plots.plot(m::unfold.UnfoldModel)  = plot_results(m.results)

function plot_results(results::DataFrame)
    Plots.plot(results.time,results.estimate,
            group=(results.term,results.group),
            legend=:outerbottom)

end


Plots.heatmap(X::SparseMatrixCSC) = heatmap(Matrix(X))
