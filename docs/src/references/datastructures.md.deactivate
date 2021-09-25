
## What is happening under the hood?
```@example Main

Xdc = designmatrix(UnfoldModel,f,evts,basisfunction)
```


Formula-Terms are wrapped with a `TimeExpandedTerm`, which upon calling `modelcols` will timeexpand the designmatrix.
There is one TimeExpandedTerm for the FixedEffects and one for each RandomEffectsTerm.
```@example Main
typeof(Xdc.formulas.rhs)
```



Visualizing the designmatrices.
Fixed Effects:
```@example Main

Plots.heatmap(Matrix(Xdc.Xs[1][1:300,:]))
```






Random Effects
```@example Main

Plots.heatmap(Matrix(Xdc.Xs[2][1:2000,1:500]))
```







And finally, generate the linear mixed model manually & fit it.
```@example Main
mf = unfoldfit(Unfold.UnfoldLinearMixedModel,Xs,data)
results = condense_long(mf)
first(results,6)
```


