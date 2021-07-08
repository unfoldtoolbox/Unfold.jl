# Overlap Correction with Linear Mixed Models

```@example Main
using StatsModels, MixedModels, DataFrames,CategoricalArrays
import Plots
using Unfold
include("../../../test/test_utilities.jl"); # function to load the simulated data
nothing;#hide
```


This notebook is similar to the `lm_tutorial`, but fits mass-univariate *mixed* models 



## Mass Univariate **Mixed** Models
Again we have 4 steps:
1. epoch the data
2. specify a formula 
3. fit a linear model to each time point & channel
4. visualize the results.


#### 1. Epoching
The data were simulated in MatLab using the `unmixed toolbox (www.unfoldtoolbox.org)` with the function`EEG_to_csv.m`.
```@example Main

data, evts = loadtestdata("testCase3",dataPath = "../../../test/data/")
data = data.+ 0.1*randn(size(data)) # we have to add minimal noise, else mixed models crashes.

transform!(evts,:subject=>categorical=>:subject);
nothing #hide
```

The `events` dataFrame has an additional column (besides being much taller): `subject`
```@example Main
first(evts,6)
```

!!! note 
        Note how small the data is! Only 12k samples, that is only ~5minutes of recording in total for 25 subjects. More realistic samples quickly take hours to fit.
        ```@example Main
        size(data)
        ```

Now we are ready to epoch the data - same as for the mass univariate, but we have more trials (nsubject more)
```@example Main
data_r = reshape(data,(1,:))
# cut the data into epochs
data_epochs,times = Unfold.epoch(data=data_r,tbl=evts,τ=(-0.4,0.8),sfreq=50);
# missing or partially missing epochs are currenlty _only_ supported for non-mixed models!
evts,data_epochs = Unfold.dropMissingEpochs(evts,data_epochs); nothing #hide
```



#### 2. Specify the formula
We define the formula. Importantly we need to specify a random effect. We are using `zerocorr` to speed up the calculation and show off that we can use it.
```@example Main
f  = @formula 0~1+condA*condB+zerocorr(1+condA*condB|subject);
```


#### 3. Fit the model
We can now run the LinearMixedModel on each time point
```@example Main
m,results = fit(UnfoldLinearMixedModel,f,evts,data_epochs,times) 
```


#### 4. Visualize results

!!! note We are working on UnfoldMakie.jl - a library to make these plots automatic and beautiful. This tutorial will be updated

Let's start with the **fixed** Effects
```@example Main
res_fixef = results[results.group.==:fixed,:]
Plots.plot(res_fixef.colname_basis,res_fixef.estimate,
        group=res_fixef.term,
        layout=1,legend=:outerbottom)
```


We see the condition effects and some residual overlap activity in the fixed effects


And now the **random** effect results
```@example Main
res_ranef = results[results.group.==:subject,:]
Plots.plot(res_ranef.colname_basis,res_ranef.estimate,
        group=res_ranef.term,
        layout=1,legend=:outerbottom)
```





The random effects are very high in areas where we simulated overlap. (i.e. <-0.1 and >0.2)

## With Overlap Correction
For overlap correction, we have to use a basis function once again.
```@example Main
basisfunction = firbasis(τ=(-0.05,.4),sfreq=40,name="basis1")
f  = @formula 0~1+condA*condB+(1+condA*condB|subject);
```



!!! warning
See the low sampling frequency? This is because the modelfit increases quadratically with the number of predictors

We can now run the mixed model.

Easy syntax: Specify formula, events, EEG-data & the basis function
```@example Main
@time mm,results = fit(UnfoldLinearMixedModel,Dict(Any=>(f,basisfunction)),evts,data); nothing #hide 
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

res_fixef = results[results.group.==:fixed,:]
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


