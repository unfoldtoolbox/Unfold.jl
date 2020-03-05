using Plots
function plot(m::UnfoldModel)
    plot(m.results.time,m.results.estimate,
            group=(m.results.term,m.results.group),
            layout=2,legend=:outerbottom)

end


Plots.heatmap(X::SparseMatrixCSC) = heatmap(Matrix(X))
