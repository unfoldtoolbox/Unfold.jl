## How To get P-Values for Mass-Univariate LMM
There are currently two ways to get p-values for LMMs: Walds t-test & likelihood ratio tests (mass univariate only).

#### Setup
```@example Main
    using Unfold, DataFrames
   include(joinpath(dirname(pathof(Unfold)), "../test/test_utilities.jl") ) # to load data
     data, evts = loadtestdata("testCase3",dataPath="../../../test/data/");
```

#### Get some 3D data / epoched
```@example Main
        # convert subject to categorical
        evts[!,:subject] .= string.(evts.subject); 
   
        # we need to add noise, else LMM crashes
        data=data.+2 .*randn(MersenneTwister(1),size(data))
        
        # epoch
        data_epoch,times = Unfold.epoch(data = data,tbl=evts,τ=(-0.1,0.3),sfreq=10);
        
        # LMMs don't like missings
        evts_epoch,data_epoch = Unfold.dropMissingEpochs(evts,data_epoch) 
```

#### Define f0 & f1 and fit!
```@example Main
        f0  = @formula 0~1+condA + (1|subject);
        f1  = @formula 0~1+condA+condB + (1+condB|subject); # could also differ in random effects
            
        m0 = fit(UnfoldModel,Dict(Any=>(f0,times)),evts_epoch,data_epoch);
        m1 = fit(UnfoldModel,Dict(Any=>(f1,times)),evts_epoch,data_epoch);
```

## Likelihoodratio
```@example Main
        uf_lrt = likelihoodratiotest(m0,m1)
        uf_lrt
```
As you can see, we have some lrts, exciting!

#### extract pvalues
```@example Main
    pvalues(uf_lrt)
```
Now we extracted the p-values. Oops, lot's of NaN's, this is the case because we do not have realistic variability in our simulation. You can interpret them as "1's" or "undefined".

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
DataFrame(:walds=>res[res.coefname .== "condB",:pvalue],:lrt=>pvals_lrt)
``` 
At least qualitatively similar - we should look for a better example here; working on UnfoldSim :-)
