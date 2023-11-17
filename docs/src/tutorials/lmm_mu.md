# [Mass Univariate Linear Mixed Models](@id lmm_massunivariate)

```@example Main

using Unfold
using UnfoldSim
using MixedModels # important to load to activate the UnfoldMixedModelsExtension
using UnfoldMakie,CairoMakie # plotting
using DataFrames
using CategoricalArrays
nothing;#hide
```
!!! important
    You have to run `using MixedModels` before or after loading Unfold to activate the MixedModels abilities!

This notebook is similar to the [Mass Univariate Linear Models (no overlap correction) tutorial](@ref lm_massunivariate) , but fits mass-univariate *mixed* models - that is, one model over all subjects, instead one model per subject. This allows to incorporate e.g. Item-effects.



## Mass Univariate **Mixed** Models
Again we have 4 steps:
1. epoch the data
2. specify a formula 
3. fit a linear model to each time point & channel
4. visualize the results.


#### 1. Epoching

```@example Main

data, evts = UnfoldSim.predef_eeg(10)
transform!(evts,:subject=>categorical=>:subject); # has to be categorical, else MixedModels.jl complains
nothing #hide
```

The `events` dataFrame has an additional column (besides being much taller): `subject`
```@example Main
first(evts,6)
```        


Now we are ready to epoch the data - same as for the mass univariate, but we have more trials (times `nsubject` more)
```@example Main
data_r = reshape(data,(1,:))
# cut the data into epochs
data_epochs,times = Unfold.epoch(data=data_r,tbl=evts,τ=(-0.4,0.8),sfreq=50);
# missing or partially missing epochs are currenlty _only_ supported for non-mixed models!
evts,data_epochs = Unfold.dropMissingEpochs(evts,data_epochs);

nothing #hide
```

#### 2. Specify the formula
We define the formula. Importantly we need to specify a random effect. We are using `zerocorr` to speed up the calculation and show off that we can use it.

```@example Main
f  = @formula 0~1+condition*continuous+zerocorr(1+condition*continuous|subject);
nothing #hide
```


#### 3. Fit the model
We can now run the LinearMixedModel on each time point
```@example Main
m = fit(UnfoldModel,f,evts,data_epochs,times)
nothing #hide
```


#### 4. Visualize results

Let's start with the **fixed** Effects
```@example Main
results = coeftable(m)

res_fixef = results[isnothing.(results.group),:]
plot_erp(res_fixef)
```


We see the condition effects and some residual overlap activity in the fixed effects


And now the **random** effect results
```@example Main
res_ranef = results[results.group.==:subject,:]
plot_erp(res_ranef)
```


### Statistics
Check out the [LMM p-value tutorial](@ref lmm_pvalues)