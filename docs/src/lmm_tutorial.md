---
author: "Benedikt Ehinger with help from Dave Kleinschmidt"
title: "Overlap Correction with Linear Mixed Models (aka unmixed.jl)"
date: 2020-02-17
---
~~~~{.julia}

using StatsModels, MixedModels, DataFrames
import Plots
using unfold
include("../test/test_utilities.jl") # to load the simulated data
~~~~~~~~~~~~~


~~~~
loadtestdata (generic function with 2 methods)
~~~~





This notebook is similar to the `lm_tutorial`, but fits mass-univariate *mixed* models and time-expanded/overlap-corrected *mixed* models.

## Reading input
The data were simulated in MatLab using the `unmixed toolbox (www.unfoldtoolbox.org)` with the function`EEG_to_csv.m`.

**Limitation**: due to current implementation in MixedModels.jl, we cannot fit overlap-corrected random effects.
That is, the `(1|item)` cannot be modelled at the moment.

~~~~{.julia}

data, evts = loadtestdata("testcase3","../test/")
data = data.+ 0.1*randn(size(data)) # we have to add minimal noise, else mixed models crashes.

categorical!(evts,:subject);
~~~~~~~~~~~~~


~~~~
7050×8 DataFrame. Omitted printing of 1 columns
│ Row  │ type   │ latency │ trialnum │ condA │ condB │ stimulus │ subject │
│      │ String │ Int64   │ Int64    │ Int64 │ Int64 │ Int64    │ Cat…    │
├──────┼────────┼─────────┼──────────┼───────┼───────┼──────────┼─────────┤
│ 1    │ sim    │ 16      │ 1        │ 0     │ 0     │ 4        │ 1       │
│ 2    │ sim    │ 36      │ 2        │ 1     │ 1     │ 10       │ 1       │
│ 3    │ sim    │ 51      │ 3        │ 0     │ 1     │ 6        │ 1       │
│ 4    │ sim    │ 66      │ 4        │ 0     │ 1     │ 8        │ 1       │
│ 5    │ sim    │ 77      │ 5        │ 1     │ 1     │ 9        │ 1       │
│ 6    │ sim    │ 91      │ 6        │ 1     │ 0     │ 13       │ 1       │
│ 7    │ sim    │ 101     │ 7        │ 1     │ 1     │ 13       │ 1       │
⋮
│ 7043 │ sim    │ 106425  │ 113      │ 0     │ 0     │ 6        │ 100     │
│ 7044 │ sim    │ 106439  │ 114      │ 0     │ 1     │ 5        │ 100     │
│ 7045 │ sim    │ 106455  │ 115      │ 0     │ 0     │ 8        │ 100     │
│ 7046 │ sim    │ 106471  │ 116      │ 1     │ 1     │ 20       │ 100     │
│ 7047 │ sim    │ 106491  │ 117      │ 0     │ 0     │ 3        │ 100     │
│ 7048 │ sim    │ 106507  │ 118      │ 1     │ 0     │ 20       │ 100     │
│ 7049 │ sim    │ 106518  │ 119      │ 0     │ 1     │ 8        │ 100     │
│ 7050 │ sim    │ 106534  │ 120      │ 0     │ 1     │ 10       │ 100     │
~~~~




The `events` dataFrame looks like this
~~~~{.julia}

first(evts,6)
~~~~~~~~~~~~~


~~~~
6×8 DataFrame. Omitted printing of 1 columns
│ Row │ type   │ latency │ trialnum │ condA │ condB │ stimulus │ subject │
│     │ String │ Int64   │ Int64    │ Int64 │ Int64 │ Int64    │ Cat…    │
├─────┼────────┼─────────┼──────────┼───────┼───────┼──────────┼─────────┤
│ 1   │ sim    │ 16      │ 1        │ 0     │ 0     │ 4        │ 1       │
│ 2   │ sim    │ 36      │ 2        │ 1     │ 1     │ 10       │ 1       │
│ 3   │ sim    │ 51      │ 3        │ 0     │ 1     │ 6        │ 1       │
│ 4   │ sim    │ 66      │ 4        │ 0     │ 1     │ 8        │ 1       │
│ 5   │ sim    │ 77      │ 5        │ 1     │ 1     │ 9        │ 1       │
│ 6   │ sim    │ 91      │ 6        │ 1     │ 0     │ 13       │ 1       │
~~~~




With the important fields being `latency`, `condA`, `condB` and `subject`.

The data are a vector.
~~~~{.julia}

println(typeof(data))
println(size(data))
~~~~~~~~~~~~~


~~~~
Array{Float64,1}
(106540,)
~~~~




**Limitation** Note how small it is! Only 12k samples, that is only ~5minutes of recording in total for 25 subjects. More realistic samples quickly take hours to fit.

## Without Overlap Correction
We define the formula
~~~~{.julia}

f  = @formula 0~1+condA*condB+(1+condA*condB|subject);
~~~~~~~~~~~~~


~~~~
FormulaTerm
Response:
  0
Predictors:
  1
  condA(unknown)
  condB(unknown)
  (condA,condB,subject)->(1 + condA * condB) | subject
  condA(unknown) & condB(unknown)
~~~~





epoch the data for the mass-univariate mixed model case
~~~~{.julia}

data_r = reshape(data,(1,:))
# cut the data into epochs
data_epochs,times = unfold.epoch(data=data_r,tbl=evts,τ=(-0.4,0.8),sfreq=50);
# missing or partially missing epochs are currenlty _only_ supported for non-mixed models!
evts,data_epochs = unfold.dropMissingEpochs(evts,data_epochs)
~~~~~~~~~~~~~


~~~~
(1, 1, 7050)(7046×8 DataFrame. Omitted printing of 1 columns
│ Row  │ type   │ latency │ trialnum │ condA │ condB │ stimulus │ subject │
│      │ String │ Int64   │ Int64    │ Int64 │ Int64 │ Int64    │ Cat…    │
├──────┼────────┼─────────┼──────────┼───────┼───────┼──────────┼─────────┤
│ 1    │ sim    │ 36      │ 2        │ 1     │ 1     │ 10       │ 1       │
│ 2    │ sim    │ 51      │ 3        │ 0     │ 1     │ 6        │ 1       │
│ 3    │ sim    │ 66      │ 4        │ 0     │ 1     │ 8        │ 1       │
│ 4    │ sim    │ 77      │ 5        │ 1     │ 1     │ 9        │ 1       │
│ 5    │ sim    │ 91      │ 6        │ 1     │ 0     │ 13       │ 1       │
│ 6    │ sim    │ 101     │ 7        │ 1     │ 1     │ 13       │ 1       │
│ 7    │ sim    │ 113     │ 8        │ 0     │ 1     │ 3        │ 1       │
⋮
│ 7039 │ sim    │ 106377  │ 110      │ 1     │ 0     │ 13       │ 100     │
│ 7040 │ sim    │ 106393  │ 111      │ 0     │ 1     │ 6        │ 100     │
│ 7041 │ sim    │ 106413  │ 112      │ 0     │ 1     │ 8        │ 100     │
│ 7042 │ sim    │ 106425  │ 113      │ 0     │ 0     │ 6        │ 100     │
│ 7043 │ sim    │ 106439  │ 114      │ 0     │ 1     │ 5        │ 100     │
│ 7044 │ sim    │ 106455  │ 115      │ 0     │ 0     │ 8        │ 100     │
│ 7045 │ sim    │ 106471  │ 116      │ 1     │ 1     │ 20       │ 100     │
│ 7046 │ sim    │ 106491  │ 117      │ 0     │ 0     │ 3        │ 100     │
, [-2.1036185545900064 -0.09620616101249467 … 0.03457492319936805 0.1594141
7693406884]

[-0.010428599490538766 0.11276513561412636 … -0.10155139247255301 13.138691
65144944]

[0.003689672020482346 -0.11497699754834818 … -0.07875501263475712 -0.102698
97300242908]

...

[0.17073360929659176 -0.032826413925233666 … -0.11495917097108686 0.0399040
39006940256]

[-0.02614791264332993 -0.008244069659384892 … 0.11096995533891257 -0.139675
377005149]

[4.810295203468854 -0.12468795475277622 … -0.16685554419753457 -0.019240882
427578737])
~~~~





We can now run the LinearMixedModel on each time point
~~~~{.julia}

m,results = unfold.fit(UnfoldLinearMixedModel,f,evts,data_epochs,times) # just "fit" without unfold should also work, but didnt in the Notebook
~~~~~~~~~~~~~


~~~~
(1, 61, 7046)(7046, 4)(1, 61, 4)(Unfold object
formula: 0 ~ 1 + condA + condB + condA & condB + (1 + condA + condB + condA
 & condB | subject)
Fields:	.modelinfo (potentially contains info on fitting procedure) 
	.beta extracted parameters 
	.X designmatrix (with fields .X.Xs, .X.events .X.formulas
, 488×7 DataFrame. Omitted printing of 2 columns
│ Row │ term          │ channel │ basisname       │ colname_basis │ estimat
e  │
│     │ Any           │ Int64   │ String          │ Float64       │ Float64
   │
├─────┼───────────────┼─────────┼─────────────────┼───────────────┼────────
───┤
│ 1   │ (Intercept)   │ 1       │ mass-univariate │ -0.4          │ 0.55582
   │
│ 2   │ (Intercept)   │ 1       │ mass-univariate │ -0.38         │ 0.90467
1  │
│ 3   │ (Intercept)   │ 1       │ mass-univariate │ -0.36         │ 1.31158
   │
│ 4   │ (Intercept)   │ 1       │ mass-univariate │ -0.34         │ 0.94732
   │
│ 5   │ (Intercept)   │ 1       │ mass-univariate │ -0.32         │ 1.07633
   │
│ 6   │ (Intercept)   │ 1       │ mass-univariate │ -0.3          │ 1.1854 
   │
│ 7   │ (Intercept)   │ 1       │ mass-univariate │ -0.28         │ 1.28006
   │
⋮
│ 481 │ condA & condB │ 1       │ mass-univariate │ 0.66          │ 1.36589
   │
│ 482 │ condA & condB │ 1       │ mass-univariate │ 0.68          │ 1.74169
   │
│ 483 │ condA & condB │ 1       │ mass-univariate │ 0.7           │ 0.68496
9  │
│ 484 │ condA & condB │ 1       │ mass-univariate │ 0.72          │ 0.48884
3  │
│ 485 │ condA & condB │ 1       │ mass-univariate │ 0.74          │ 0.01856
26 │
│ 486 │ condA & condB │ 1       │ mass-univariate │ 0.76          │ 0.60035
6  │
│ 487 │ condA & condB │ 1       │ mass-univariate │ 0.78          │ 0.19045
5  │
│ 488 │ condA & condB │ 1       │ mass-univariate │ 0.8           │ 0.85210
1  │)
~~~~





### Fixed Effects
~~~~{.julia}

res_fixef = results[results.group.=="fixed",:]
Plots.plot(res_fixef.colname_basis,res_fixef.estimate,
        group=res_fixef.term,
        layout=1,legend=:outerbottom)
~~~~~~~~~~~~~


~~~~
Error: MethodError: no method matching union()
Closest candidates are:
  union(!Matched::DataStructures.IntSet, !Matched::Any) at C:\Users\behinge
r\.julia\packages\DataStructures\DLSxi\src\int_set.jl:96
  union(!Matched::DataStructures.SparseIntSet, !Matched::Any) at C:\Users\b
ehinger\.julia\packages\DataStructures\DLSxi\src\sparse_int_set.jl:151
  union(!Matched::BitSet, !Matched::Any...) at bitset.jl:312
  ...
~~~~




We see the condition effects and some residual overlap activity in the fixed effects

### Random Effects
And the random effect results
~~~~{.julia}

res_ranef = results[results.group.=="subject",:]
Plots.plot(res_ranef.colname_basis,res_ranef.estimate,
        group=res_ranef.term,
        layout=1,legend=:outerbottom)
~~~~~~~~~~~~~


~~~~
Error: MethodError: no method matching union()
Closest candidates are:
  union(!Matched::DataStructures.IntSet, !Matched::Any) at C:\Users\behinge
r\.julia\packages\DataStructures\DLSxi\src\int_set.jl:96
  union(!Matched::DataStructures.SparseIntSet, !Matched::Any) at C:\Users\b
ehinger\.julia\packages\DataStructures\DLSxi\src\sparse_int_set.jl:151
  union(!Matched::BitSet, !Matched::Any...) at bitset.jl:312
  ...
~~~~




The random effects are very high in areas where we simulated overlap. (i.e. <-0.1 and >0.2)

## With Overlap Correction
For overlap correction, we have to use a basis function once again.

~~~~{.julia}

basisfunction = firbasis(τ=(-0.05,.4),sfreq=40)
f  = @formula 0~1+condA*condB+(1+condA*condB|subject);
~~~~~~~~~~~~~


~~~~
Error: UndefKeywordError: keyword argument name not assigned
~~~~





**Limitation:** Currently we cannot model correlation between time-points or random slopes.

**Limitation:** See the low sampling frequency? This is because the modelfit increases quadratically with the number of predictors

We can now run the mixed model.

Easy syntax: Specify formula, events, EEG-data & the basis function
~~~~{.julia}

@time mm,results = unfold.fit(UnfoldLinearMixedModel,f,evts,data,basisfunction) # just "fit" without unfold should also work, but didnt in the Notebook
~~~~~~~~~~~~~


~~~~
Error: UndefVarError: basisfunction not defined
~~~~




We receive an object containing the (very large) mixed model:
~~~~{.julia}

show(coeftable(mm.modelinfo))
~~~~~~~~~~~~~


~~~~
Error: UndefVarError: mm not defined
~~~~





But again, we also get a *tidy*-dataframe with the results
~~~~{.julia}

first(results,6)
~~~~~~~~~~~~~


~~~~
6×7 DataFrame. Omitted printing of 2 columns
│ Row │ term        │ channel │ basisname       │ colname_basis │ estimate 
│
│     │ Any         │ Int64   │ String          │ Float64       │ Float64  
│
├─────┼─────────────┼─────────┼─────────────────┼───────────────┼──────────
┤
│ 1   │ (Intercept) │ 1       │ mass-univariate │ -0.4          │ 0.55582  
│
│ 2   │ (Intercept) │ 1       │ mass-univariate │ -0.38         │ 0.904671 
│
│ 3   │ (Intercept) │ 1       │ mass-univariate │ -0.36         │ 1.31158  
│
│ 4   │ (Intercept) │ 1       │ mass-univariate │ -0.34         │ 0.94732  
│
│ 5   │ (Intercept) │ 1       │ mass-univariate │ -0.32         │ 1.07633  
│
│ 6   │ (Intercept) │ 1       │ mass-univariate │ -0.3          │ 1.1854   
│
~~~~





and thus we can easily plot the fixed effect results.
~~~~{.julia}

res_fixef = results[results.group.=="fixed",:]
Plots.plot(res_fixef.colname_basis,res_fixef.estimate,
        group=res_fixef.term,
        layout=1,legend=:outerbottom)
~~~~~~~~~~~~~


~~~~
Error: MethodError: no method matching union()
Closest candidates are:
  union(!Matched::DataStructures.IntSet, !Matched::Any) at C:\Users\behinge
r\.julia\packages\DataStructures\DLSxi\src\int_set.jl:96
  union(!Matched::DataStructures.SparseIntSet, !Matched::Any) at C:\Users\b
ehinger\.julia\packages\DataStructures\DLSxi\src\sparse_int_set.jl:151
  union(!Matched::BitSet, !Matched::Any...) at bitset.jl:312
  ...
~~~~





And the random effect results.
~~~~{.julia}

res_ranef = results[results.group.=="subject",:]
Plots.plot(res_ranef.colname_basis,res_ranef.estimate,
        group=res_ranef.term,
        layout=1,legend=:outerbottom)
~~~~~~~~~~~~~


~~~~
Error: MethodError: no method matching union()
Closest candidates are:
  union(!Matched::DataStructures.IntSet, !Matched::Any) at C:\Users\behinge
r\.julia\packages\DataStructures\DLSxi\src\int_set.jl:96
  union(!Matched::DataStructures.SparseIntSet, !Matched::Any) at C:\Users\b
ehinger\.julia\packages\DataStructures\DLSxi\src\sparse_int_set.jl:151
  union(!Matched::BitSet, !Matched::Any...) at bitset.jl:312
  ...
~~~~






## What is happening under the hood?
~~~~{.julia}

Xdc = designmatrix(UnfoldLinearMixedModel,f,evts,basisfunction)
~~~~~~~~~~~~~


~~~~
Error: UndefVarError: basisfunction not defined
~~~~





Formula-Terms are wrapped with a `TimeExpandedTerm`, which upon calling `modelcols` will timeexpand the designmatrix.
There is one TimeExpandedTerm for the FixedEffects and one for each RandomEffectsTerm.

~~~~{.julia}

typeof(Xdc.formulas.rhs)
~~~~~~~~~~~~~


~~~~
Error: UndefVarError: Xdc not defined
~~~~





Visualizing the designmatrices.
Fixed Effects:
~~~~{.julia}

Plots.heatmap(Matrix(Xdc.Xs[1][1:300,:]))
~~~~~~~~~~~~~


~~~~
Error: UndefVarError: Xdc not defined
~~~~





Random Effects
~~~~{.julia}

Plots.heatmap(Matrix(Xdc.Xs[2][1:2000,1:500]))
~~~~~~~~~~~~~


~~~~
Error: UndefVarError: Xdc not defined
~~~~






And finally, generate the linear mixed model manually & fit it.
~~~~{.julia}

mf = unfoldfit(unfold.UnfoldLinearMixedModel,Xs,data)
results = condense_long(mf)
first(results,6)
~~~~~~~~~~~~~




## Summary
There are four different model types currently "fitable"

1. Timeexpansion **No**, Mixed **No**  : `fit(UnfoldLinearModel,f,evts,data_epoch,times)`
1. Timeexpansion **Yes**, Mixed **No** : `fit(UnfoldLinearModel,f,evts,data,basisfunction)`
1. Timeexpansion **No**, Mixed **Yes** : `fit(UnfoldLinearMixedModel,f,evts,data_epoch,times)`
1. Timeexpansion **Yes**, Mixed **Yes**: `fit(UnfoldLinearMixedModel,f,evts,data,basisfunction)`


