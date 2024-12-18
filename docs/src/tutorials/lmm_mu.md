# [Mass Univariate Linear Mixed Models](@id lmm_massunivariate)

```@example Main

using Unfold
using UnfoldSim
using MixedModels # important to load to activate the UnfoldMixedModelsExtension
using UnfoldMakie, CairoMakie # plotting
using DataFrames
using CategoricalArrays
nothing;#hide
```

!!! important
    You have to run `using MixedModels` before or after loading Unfold to activate the MixedModels abilities!

This notebook is similar to the [Mass Univariate Linear Models (no overlap correction) tutorial](@ref lm_massunivariate), but fits mass-univariate *mixed* models - that is, one model over all subjects, instead of one model per subject. This allows to include item effects, for example.

## Mass Univariate **Mixed** Models

Again we have 4 steps:

1. Split data into epochs
2. Specify a formula
3. Fit a linear model to each time point & channel
4. Visualize the results.

### 1. Epoching

```@example Main
data, evts = UnfoldSim.predef_eeg(10; return_epoched = true) # simulate 10 subjects
data = reshape(data, 1, size(data, 1), :) # concatenate the data into a long EEG dataset
times = range(0, length = size(data, 2), step = 1 / 100)
transform!(evts, :subject => categorical => :subject); # :subject must be categorical, otherwise MixedModels.jl complains
nothing #hide
```

The `events` dataFrame has an additional column (besides being much taller): `subject`

```@example Main
first(evts, 6)
```

### 2. Formula specification

We define the formula. Importantly, we need to specify a random effect. We use `zerocorr` to speed up the calculation.

```@example Main
f = @formula 0 ~ 1 + condition * continuous + zerocorr(1 + condition * continuous | subject);
nothing #hide
```

### 3. Model fitting

We can now run the LinearMixedModel at each time point.

```@example Main
m = fit(UnfoldModel, f, evts, data, times)
nothing #hide
```

### 4. Visualization of results

Let's start with the **fixed** effects.
We see the condition effects and some residual overlap activity in the fixed effects.

```@example Main
results = coeftable(m)

res_fixef = results[isnothing.(results.group), :]
plot_erp(res_fixef)
```

And now comes the **random** effect:

```@example Main
res_ranef = results[results.group .== :subject, :]
plot_erp(res_ranef)
```

### Statistics

Check out the [LMM p-value tutorial](@ref lmm_pvalues)
