# Unfold Documentation
Documentation is currently written.

Best to start with [installation](@Ref), then go for the mass-univariate approach [lm_mu.md](@Ref), which should be familiar to you if you did ERPs before. Then the overlap-correction [lm_overlap.md](@Ref), LMM (to-be-done) and non-linear (to-be-done) would be good stations.

In case you want to understand the tools better, check out our explanations.

Once you are familiar with the tools, check out further how-to guides for specific applications.

In case you want to understand the toolbox better, we plan to offer technical references. This includes Benchmarks & Explorations.


## Summary
There are four different model types currently "fitable"

1. Timeexpansion **No**, Mixed **No**  : `fit(UnfoldModel,Dict(Any=>(f,-0.1:0.01:0.5)),evts,data_epoch)`
1. Timeexpansion **Yes**, Mixed **No** : `fit(UnfoldModel,Dict(Any=>(f,basisfunction)),evts,data)`
1. Timeexpansion **No**, Mixed **Yes** : `fit(UnfoldModel,Dict(Any=>(fLMM,-0.1:0.01:0.5)),evts,data_epoch)`
1. Timeexpansion **Yes**, Mixed **Yes**: `fit(UnfoldModel,Dict(Any=>(fLMM,basisfunction)),evts,data)`

With
`f = @formula 0~1+condition`
`fLMM = @formula 0~1+condition+(1|subject) + (1|item)`
`basisfunction = firbasis(τ=(-0.1,0.5),sfreq=100"))`




## TOC
```@contents
```

