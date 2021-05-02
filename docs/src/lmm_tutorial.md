---
author: "Benedikt Ehinger with help from Dave Kleinschmidt"
title: "Overlap Correction with Linear Mixed Models"
date: 2021-04-29
---

```@example Main


using StatsModels, MixedModels, DataFrames
import Plots
using Unfold
include("../../test/test_utilities.jl"); # function to load the simulated data
```





This notebook is similar to the `lm_tutorial`, but fits mass-univariate *mixed* models and time-expanded/overlap-corrected *mixed* models.

## Reading input
The data were simulated in MatLab using the `unmixed toolbox (www.unfoldtoolbox.org)` with the function`EEG_to_csv.m`.

**Limitation**: due to current implementation in MixedModels.jl, we cannot fit overlap-corrected random effects.
That is, the `(1|item)` cannot be modelled at the moment.

```@example Main

data, evts = loadtestdata("testCase3",dataPath = "../../test/")
data = data.+ 0.1*randn(size(data)) # we have to add minimal noise, else mixed models crashes.

categorical!(evts,:subject);
```



The `events` dataFrame looks like this
```@example Main

first(evts,6)
```



With the important fields being `latency`, `condA`, `condB` and `subject`.

The data are a vector.
```@example Main

println(typeof(data))
println(size(data))
```





**Limitation** Note how small it is! Only 12k samples, that is only ~5minutes of recording in total for 25 subjects. More realistic samples quickly take hours to fit.

## Without Overlap Correction
We define the formula
```@example Main
f  = @formula 0~1+condA*condB+(1+condA*condB|subject);
```






epoch the data for the mass-univariate mixed model case
```@example Main
data_r = reshape(data,(1,:))
# cut the data into epochs
data_epochs,times = Unfold.epoch(data=data_r,tbl=evts,τ=(-0.4,0.8),sfreq=50);
# missing or partially missing epochs are currenlty _only_ supported for non-mixed models!
evts,data_epochs = Unfold.dropMissingEpochs(evts,data_epochs)
```




We can now run the LinearMixedModel on each time point
```@example Main
m,results = fit(UnfoldLinearMixedModel,f,evts,data_epochs,times) 
```


### Fixed Effects
```@example Main
res_fixef = results[results.group.=="fixed",:]
Plots.plot(res_fixef.colname_basis,res_fixef.estimate,
        group=res_fixef.term,
        layout=1,legend=:outerbottom)
```


We see the condition effects and some residual overlap activity in the fixed effects

### Random Effects
And the random effect results
```@example Main
res_ranef = results[results.group.=="subject",:]
Plots.plot(res_ranef.colname_basis,res_ranef.estimate,
        group=res_ranef.term,
        layout=1,legend=:outerbottom)
```





The random effects are very high in areas where we simulated overlap. (i.e. <-0.1 and >0.2)

## With Overlap Correction
For overlap correction, we have to use a basis function once again.
```@example Main
basisfunction = firbasis(τ=(-0.05,.4),sfreq=40)
f  = @formula 0~1+condA*condB+(1+condA*condB|subject);
```



**Limitation:** Currently we cannot model correlation between time-points or random slopes.

**Limitation:** See the low sampling frequency? This is because the modelfit increases quadratically with the number of predictors

We can now run the mixed model.

Easy syntax: Specify formula, events, EEG-data & the basis function
```@example Main
@time mm,results = fit(UnfoldLinearMixedModel,f,evts,data,basisfunction) 
```





We receive an object containing the (very large) mixed model:
```@example Main
show(coeftable(mm.modelinfo))
```



But again, we also get a *tidy*-dataframe with the results
```@example Main
first(results,6)
```



and thus we can easily plot the fixed effect results.
```@example Main

res_fixef = results[results.group.=="fixed",:]
Plots.plot(res_fixef.colname_basis,res_fixef.estimate,
        group=res_fixef.term,
        layout=1,legend=:outerbottom)
```





And the random effect results.
```@example Main

res_ranef = results[results.group.=="subject",:]
Plots.plot(res_ranef.colname_basis,res_ranef.estimate,
        group=res_ranef.term,
        layout=1,legend=:outerbottom)
```

## What is happening under the hood?
```@example Main

Xdc = designmatrix(UnfoldLinearMixedModel,f,evts,basisfunction)
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




## Summary
There are four different model types currently "fitable"

1. Timeexpansion **No**, Mixed **No**  : `fit(UnfoldLinearModel,f,evts,data_epoch,times)`
1. Timeexpansion **Yes**, Mixed **No** : `fit(UnfoldLinearModel,f,evts,data,basisfunction)`
1. Timeexpansion **No**, Mixed **Yes** : `fit(UnfoldLinearMixedModel,f,evts,data_epoch,times)`
1. Timeexpansion **Yes**, Mixed **Yes**: `fit(UnfoldLinearMixedModel,f,evts,data,basisfunction)`


