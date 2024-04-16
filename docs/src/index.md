# Unfold Documentation
Documentation is currently written.

If you want to follow the **tutorials**, best to start with the [mass-univariate approach](@ref lm_massunivariate), which should be familiar to you if you did ERPs before. Then the [overlap-correction tutorial](@ref lm_overlap), [mixed mass univariate](@ref lmm_massunivariate), [mixed overlap (tricky!)](@ref lmm_overlap). If you are then not satisfied, check out more advanced topics: [effects-interface (aka what to do after fitting)](@ref effects), or non-linear effects.

In case you want to understand the tools better, check out our **explanations**.

Once you are familiar with the tools, check out further **how-to guides** for specific applications.

In case you want to understand the toolbox better, we plan to offer **technical references**. This includes Benchmarks & Explorations.


## Quick start
There are four different model types currently "fitable"

1. Timeexpansion **No**, Mixed **No**  : `fit(UnfoldModel, [Any=>(f, -0.1:0.01:0.5)], evts, data_epoch)`
1. Timeexpansion **Yes**, Mixed **No** : `fit(UnfoldModel, [Any=>(f, basisfunction)], evts, data)`
1. Timeexpansion **No**, Mixed **Yes** : `fit(UnfoldModel, [Any=>(fLMM, -0.1:0.01:0.5)], evts, data_epoch)`
1. Timeexpansion **Yes**, Mixed **Yes**: `fit(UnfoldModel, [Any=>(fLMM, basisfunction)], evts, data)`

With

```julia
f = @formula 0 ~ 1 + condition
fLMM = @formula 0 ~ 1 + condition + (1|subject) + (1|item)
basisfunction = firbasis(τ = (-0.1,0.5), sfreq = 100))
```

