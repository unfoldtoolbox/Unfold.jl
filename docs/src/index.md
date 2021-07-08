# Unfold Documentation
Documentation is currently written.

Best to start with [installation](@Ref), then go for the mass-univariate approach [lm_mu.md](@Ref), which should be familiar to you if you did ERPs before. Then the overlap-correction [lm_overlap.md](@Ref), LMM (to-be-done) and non-linear (to-be-done) would be good stations.

In case you want to understand the tools better, check out our explanations.

Once you are familiar with the tools, check out further how-to guides for specific applications.

In case you want to understand the toolbox better, we plan to offer technical references. This includes Benchmarks & Explorations.

## Summary


## Summary
There are four different model types currently "fitable"

1. Timeexpansion **No**, Mixed **No**  : `fit(UnfoldLinearModel,f,evts,data_epoch,times)`
1. Timeexpansion **Yes**, Mixed **No** : `fit(UnfoldLinearModel,f,evts,data,basisfunction)`
1. Timeexpansion **No**, Mixed **Yes** : `fit(UnfoldLinearMixedModel,f,evts,data_epoch,times)`
1. Timeexpansion **Yes**, Mixed **Yes**: `fit(UnfoldLinearMixedModel,f,evts,data,basisfunction)`




## TOC
```@contents
```

