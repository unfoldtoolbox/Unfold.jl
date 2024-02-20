# [Overlap Correction with Linear Mixed Models](@id lmm_overlap)

```@example Main

using Unfold
using UnfoldSim

using CategoricalArrays
using MixedModels
using UnfoldMakie,CairoMakie
using DataFrames

nothing;#hide
```


This notebook is similar to the Linear Model with Overlap Correction tutorial, but fits **mixed** models with overlap correction

!!! warning 
    **Limitation**: This is not production ready at all. Still lot's of things to find out and tinker with. Don't use this if you did not look under the hood of the toolbox!

!!! important
    Even worse - it is right now not working. Some new bug appeared

## Get some  data

```@example Main
dat,evts = UnfoldSim.predef_2x2(;signalsize=20,n_items=16,n_subjects=16)
dat = dat[:] # we need all data concatenated over subjects
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
basisfunction = firbasis(τ=(-0.4,.8),sfreq=20,name="stimulus")
nothing #hide
```




#### 2. Specify the formula
We define the formula. Importantly we need to specify a random effect. 

!!! note
    We are using `zerocorr` because we need it here, else the model will try to model all correlations between all timepoints and all factors!

```@example Main
f  = @formula 0~1+A*B+zerocorr(1+A*B|subject);
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
results.group = string.(results.group) # this will not be necessary in future versions
plot_erp(results;mapping=(;col=:group))
```
