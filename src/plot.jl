# These functions are not included yet because I am not yet sure which plotting library to put my bets on.
# Makie&StatsMakie could be supercool, but I ran into some issues

using AlgebraOfGraphics, Makie

#---
#Plots.plot(m::Unfold.UnfoldModel)  = plot_results(m.results)



function plot_results(results::DataFrame;y=:estimate,color=:term,layout=:group,stderror=false,pvalue = DataFrame(:from=>[],:to=>[],:pval=>[]))
    results = deepcopy(results)
    results.group[isnothing.(results.group)] .= :fixef
    m = mapping(:colname_basis,y,color=color,layout=layout)


    basic = AlgebraOfGraphics.data(results) * visual(Lines) * m

    if stderror
        res_se = copy(results)
        res_se = res_se[.!isnothing.(res_se.stderror),:]
        res_se[!,:se_low] = res_se[:,y].-res_se.stderror
        res_se[!,:se_high] = res_se[:,y].+res_se.stderror
        basic =  AlgebraOfGraphics.data(res_se)*visual(Band,alpha=0.5)*mapping(:colname_basis,:se_low,:se_high,color=:term,layout=layout) + basic
    end
    
    d = basic |> draw

    # add the pvalues
    if !isempty(pvalue)

        yval = d.grid[1,1].scales[2].extrema|>x->-0.01*(x[2]-x[1])+x[1]
        x = [Point(x,yval)=>Point(y,yval) for (x,y) in zip(pvalue.from,pvalue.to)]
        linesegments!(d.grid[1,1].axis,x,linewidth=4) # assumes first one is where we want to plot. Limitation!

    end
    return d
    
end


#using SparseArrays
#Plots.heatmap(X::SparseMatrixCSC) = heatmap(Matrix(X))
