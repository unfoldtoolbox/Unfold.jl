# [Overlap Correction with Linear Mixed Models](@id lmm_overlap)

```@example Main

using Unfold
using UnfoldSim

using CategoricalArrays
using MixedModels
using UnfoldMakie, CairoMakie
using DataFrames

nothing;#hide
```

This notebook is similar to the Linear Model with Overlap Correction tutorial, but fits **mixed** models with overlap correction

!!! warning
    **Limitation**: This functionality is not ready for general use. There are still a lot of things to find out and tinker with. Don't use this if you haven't looked under the hood of the toolbox! Be aware of crashes / timeouts for non-trivial problems

## Get some data

```@example Main
dat, evts = UnfoldSim.predef_2x2(; signalsize=20, n_items=16, n_subjects=16)

# We also need to fix the latencies, they are now relative to 1:size(data, 1), but we want a continuous long EEG.
subj_idx = [parse(Int, split(string(s), 'S')[2]) for s in evts.subject]
evts.latency .+= size(dat, 1) .* (subj_idx .- 1)

dat = dat[:] # we need all data concatenated over subjects
evts.subject  = categorical(Array(evts.subject))
nothing #hide
```

## Linear **Mixed** Model Continuous Time

Again we have 4 steps:

1. Specify a temporal basisfunction
2. Specify a formula
3. Fit a linear model for each channel (one model for all timepoints!)
4. Visualize the results.

### 1. Specify a temporal basisfunction

By default, we would want to use a FIR basis function. See [Basis Functions](@ref) for more details.

```@example Main
basisfunction = firbasis(τ=(-0.4, .8), sfreq=20, name="stimulus")
nothing #hide
```

### 2. Specify the formula

Define the formula and specify a random effect.

!!! note
    We use `zerocorr` to prevent the model from computing all correlations between all timepoints and factors.

```@example Main
f  = @formula 0 ~ 1 + A  *B + zerocorr(1 + A*B|subject);
```

### 3. Fit the model

```@example Main
bfDict = Dict(Any=>(f, basisfunction))
# Skipping this tutorial for now due to a significant error.
m = fit(UnfoldModel, bfDict, evts, dat)

results = coeftable(m)
first(results, 6)
```

### 4. Visualize results

```@example Main
plot_erp(results; mapping=(; col = :group))
```
