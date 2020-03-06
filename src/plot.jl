using Plots
Plots.plot(m::unfold.UnfoldModel)  = plot_results(m.results)

function plot_results(results::DataFrame)
    plot(results.time,results.estimate,
            group=(results.term,results.group),
            layout=2,legend=:outerbottom)

end


Plots.heatmap(X::SparseMatrixCSC) = heatmap(Matrix(X))
