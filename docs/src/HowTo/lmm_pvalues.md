## [How To get P-Values for Mass-Univariate LMM](@id lmm_pvalues)    
There are currently two ways to obtain p-values for LMMs: Wald's t-test and likelihood ratio tests (mass univariate only).

#### Setup
```@example Main
using MixedModels, Unfold # we require to load MixedModels to load the PackageExtension
using DataFrames
using UnfoldSim
using CairoMakie
data_epoch, evts =
    UnfoldSim.predef_2x2(; n_items = 52, n_subjects = 40, return_epoched = true)
data_epoch = reshape(data_epoch, size(data_epoch, 1), :) # 
times = range(0, 1, length = size(data_epoch, 1))
```

#### Define f0 & f1 and fit!
```@example Main
f0 = @formula 0 ~ 1 + A + (1 + A | subject);
f1 = @formula 0 ~ 1 + A + B + (1 + A | subject); # could also differ in random effects

m0 = fit(UnfoldModel, Dict(Any => (f0, times)), evts, data_epoch);
m1 = fit(UnfoldModel, Dict(Any => (f1, times)), evts, data_epoch);
```

## Likelihood ratio
```@example Main
uf_lrt = likelihoodratiotest(m0, m1)
uf_lrt[1]
```
As you can see, we have some likelihood ratio outcomes, exciting!

#### Extract p-values
```@example Main
pvalues(uf_lrt)
```
We have extracted the p-values and now need to make them usable.     The solution can be found in the documentation under `?pvalues`.
```@example Main
pvals_lrt = vcat(pvalues(uf_lrt)...)
nchan = 1
ntime = length(times)
reshape(pvals_lrt, ntime, nchan)' # note the last transpose via ' !
```

Perfecto, these are the LRT p-values of a model `condA` vs. `condA+condB` with same random effect structure.

## Walds T-Test
This method is easier to calculate but has limitations in accuracy and scope. It may also be less accurate due to estimation of degrees of freedom. Testing is limited in this case, as random effects cannot be tested and only single predictors can be used, which may not be appropriate for spline effects. It is important to note that this discussion is beyond the scope of this LMM package. 

```@example Main
res = coeftable(m1)
# only fixed effects: what is not in a ranef group is a fixef.
res = res[isnothing.(res.group), :] 
# calculate t-value
res[:, :tvalue] = res.estimate ./ res.stderror
``` 

We obtained Walds t, but how to translate to a p-value?

Determining the necessary degrees of freedom for the t-distribution is a complex issue with much debate surrounding it. 
One approach is to use the number of subjects as an upper bound for the p-value (your df will be between $n_{subject}$ and $\sum{n_{trials}}$).

```@example Main
df = length(unique(evts.subject))
```
Plug it into the t-distribution. 

```@example Main
using Distributions
res.pvalue = pdf.(TDist(df),res.tvalue)
```

## Comparison of methods
Cool! Let's compare both methods of p-value calculation!
```@example Main
df = DataFrame(:walds => res[res.coefname.=="B: b_tiny", :pvalue], :lrt => pvals_lrt)
f = Figure()
scatter(f[1, 1], times, res[res.coefname.=="B: b_tiny", :estimate])
scatter(f[1, 2], df.walds, df.lrt)
scatter(f[2, 1], times, df.walds)
scatter(f[2, 2], times, df.lrt)

f
``` 
The figures appears to be similar! Note that the Wald's T-test is considered more liberal than the LRT. It is recommended to use the upcoming `MixedModelsPermutations.jl` package or, alternatively, to use `KenwardRoger` in R (not yet published).
