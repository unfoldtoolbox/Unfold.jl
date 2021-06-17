# These functions are not included yet because I am not yet sure which plotting library to put my bets on.
# Makie&StatsMakie could be supercool, but I ran into some issues

using AlgebraOfGraphics, Makie

#---
#Plots.plot(m::Unfold.UnfoldModel)  = plot_results(m.results)

function plot_results(results::DataFrame;y=:estimate,color=:term,layout_x=:group,stderror=false,pvalue = DataFrame(:from=>[],:to=>[],:pval=>[]))
    m = mapping(:colname_basis,y,color=color,layout_x=layout_x)

    basic = AlgebraOfGraphics.data(results) * visual(Lines) * m

    if stderror
        res_se = copy(results)
        res_se = res_se[.!isnothing.(res_se.stderror),:]
        res_se[!,:se_low] = res_se[:,y].-res_se.stderror
        res_se[!,:se_high] = res_se[:,y].+res_se.stderror
        basic =  AlgebraOfGraphics.data(res_se)*visual(Band,alpha=0.5)*mapping(:colname_basis,:se_low,:se_high,color=:term,layout_x=layout_x) + basic
    end
    
    d = basic |> draw

    # add the pvalues
    if !isempty(pvalue)
        x = [Point(x,0.)=>Point(y,0.) for (x,y) in zip(pvalue.from,pvalue.to)]
        linesegments!(d.children[1],x,linewidth=2) # assumes first one is where we want to plot. Limitation!

    end
    return d
    
end


#using SparseArrays
#Plots.heatmap(X::SparseMatrixCSC) = heatmap(Matrix(X))
