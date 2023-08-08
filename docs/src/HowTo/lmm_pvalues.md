## [How To get P-Values for Mass-Univariate LMM](@id lmm_pvalues)
There are currently two ways to get p-values for LMMs: Walds t-test & likelihood ratio tests (mass univariate only).

#### Setup
```@example Main
using MixedModels, Unfold # we require to load MixedModels to load the PackageExtension
using DataFrames
using UnfoldSim
data_epoch,evts = UnfoldSim.predef_2x2(;n_items=52,n_subjects=40,return_epoched=true)
data_epoch = reshape(data_epoch,1,size(data_epoch)... )
times = range(0,1,length=size(data_epoch,1))
```

#### Define f0 & f1 and fit!
```@example Main
        f0  = @formula 0~1+A + (1+A|subject);
        f1  = @formula 0~1+A+B + (1+A|subject); # could also differ in random effects
            
        m0 = fit(UnfoldModel,Dict(Any=>(f0,times)),evts,data_epoch);
        m1 = fit(UnfoldModel,Dict(Any=>(f1,times)),evts,data_epoch);
```

## Likelihoodratio
```@example Main
        uf_lrt = likelihoodratiotest(m0,m1)
        uf_lrt[1]
```
As you can see, we have some lrts, exciting!

#### extract pvalues
```@example Main
    pvalues(uf_lrt)
```
Now we extracted the p-values. 

Let's look at the structure. How do we get them to something usable? The solution is in the docs `?pvalues`

```@example Main
    pvals_lrt = vcat(pvalues(uf_lrt)...)
    nchan = 1
    ntime = length(times)
    reshape(pvals_lrt,ntime,nchan)' # note the last transpose via ' !
```

Perfecto, these are the LRT pvalues of a model `condA` vs. `condA+condB` with same random effect structure.

## Walds T-Test
This method is easier to calculate, but more limited and likely less accurate (due to degrees of freedom estimation, but that is a discussion for a proper LMM package). It is further limited, as e.g. random effects cannot be tested, and you can only test single predictors, which might make no sense for e.g. spline-effects

```@example Main
res = coeftable(m1)
# only fixef (what is not in a ranef group is a fixef)
res = res[isnothing.(res.group),:] 
#calculate t-value
res[:,:tvalue] = res.estimate ./ res.stderror
``` 

Okay cool, we have walds-t, but how to translate to a p-value?

This is by no means trivial, as the necessary degrees of freedom of the t-distribution are not well specified. There is a huge debate how to do it.
One simple way is to simply assume the number of subjects as an upper bound to your p-value (your df will be between $n_{subject}$ and $\sum{n_{trials}}$)

```@example Main
    df = length(unique(evts.subject))
```
Plug it into the t-distribution
```@example Main
    using Distributions
    res.pvalue = pdf.(TDist(df),res.tvalue)

```

Cool! Let's compare!
```@example Main
df = DataFrame(:walds=>res[res.coefname .== "B: b_tiny",:pvalue],:lrt=>pvals_lrt)
f = Figure()
scatter(f[1,1],times,res[res.coefname .== "B: b_tiny",:estimate])
scatter(f[1,2],df.walds,df.lrt)
scatter(f[2,1],times,df.walds)
scatter(f[2,2],times,df.lrt)

f
``` 
Look pretty similar! note that the Walds-T is typically too liberal (LRT also, but to a lesser exted). Best is to use the forthcoming MixedModelsPermutations.jl or go the route via R and use KenwardRoger (data not yet published)
