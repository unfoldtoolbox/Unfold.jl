# Mass Univariate Linear Mixed Models

```@example Main
using StatsModels, MixedModels, DataFrames,CategoricalArrays

using Unfold
using UnfoldMakie,CairoMakie
using DataFrames
include(joinpath(dirname(pathof(Unfold)), "../test/test_utilities.jl") ) # to load data
nothing;#hide
```


This notebook is similar to the [Mass Univariate Linear Models (no overlap correction) tutorial](@ref lm_massunivariate) , but fits mass-univariate *mixed* models 



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

transform!(evts,:subject=>categorical=>:subject); # has to be categorical, else MixedModels.jl complains
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
evts,data_epochs = Unfold.dropMissingEpochs(evts,data_epochs);
nothing #hide
```

#### 2. Specify the formula
We define the formula. Importantly we need to specify a random effect. We are using `zerocorr` to speed up the calculation and show off that we can use it.
```@example Main
f  = @formula 0~1+condA*condB+zerocorr(1+condA*condB|subject);
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


The random effects are very high in areas where we simulated overlap. (i.e. <-0.1 and >0.2)

### Statistics
Check out the [LMM p-value tutorial](@ref lmm_pvalues)