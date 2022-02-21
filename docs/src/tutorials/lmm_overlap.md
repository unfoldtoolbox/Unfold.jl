# Overlap Correction with Linear Mixed Models

```@example Main
using StatsModels, MixedModels, DataFrames,CategoricalArrays

using Unfold
using UnfoldMakie,CairoMakie
using DataFrames
include(joinpath(dirname(pathof(Unfold)), "../test/test_utilities.jl") ) # to load data
nothing;#hide
```


This notebook is similar to the [Linear Model with Overlap Correction](@ref) tutorial, but fits **mixed** models with overlap correction

!!! warning 
    **Limitation**: This is not production ready at all. Still lot's of things to find out and tinker with. Don't use this if you did not look under the hood of the toolbox!

## Load the data

```@example Main
dat, evts = loadtestdata("testCase3",dataPath = "../../../test/data/")
dat = dat .+ 0.1*randn(size(dat)) # we have to add minimal noise, else mixed models crashes.
dat = dat'
evts.subject  = categorical(Array(evts.subject))
nothing #hide
```


## Linear **Mixed** Model Continuous Time
Again we have 4 steps:
1. specify a temporal basisfunction
2. specify a formula
3. fit a linear model for each channel (one for all timepoints!)
4. visualize the results.

#### 1. specify a temporal basisfunction
By default, we would want to use a FIR basisfunction. See [Basis Functions](@ref) for more details.
```@example Main
basisfunction = firbasis(τ=(-0.4,.8),sfreq=50,name="stimulus")
nothing #hide
```




#### 2. Specify the formula
We define the formula. Importantly we need to specify a random effect. 

!!! note
    We are using `zerocorr` because we need it here, else the model will try to model all correlations between all timepoints and all factors!

```@example Main
f  = @formula 0~1+condA*condB+zerocorr(1+condA*condB|subject);
```


#### 3. Fit the model
```@example Main
bfDict = Dict(Any=>(f,basisfunction))
# for some reason this results in a big error. Skipping this tutorial right now
m = fit(UnfoldModel,bfDict,evts,dat) 
results = coeftable(m)
first(results,6)
```


#### 4. Visualize results

```@example Main
plot_results(results)
```
