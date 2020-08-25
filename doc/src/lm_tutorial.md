---
author: "Benedikt Ehinger, with help Dave Kleinschmidt"
title: "Overlap Correction with Linear Models (aka unfold.jl)"
date: 2020-06-07
---
~~~~{.julia}

using StatsModels, MixedModels, DataFrames
import DSP.conv
import Plots
using unfold
include("../test/test_utilities.jl") # to load the simulated data
~~~~~~~~~~~~~


~~~~
loadtestdata (generic function with 2 methods)
~~~~





In this notebook we will fit regression models to (simulated) EEG data. We will see that we need some type of overlap correction, as the events are close in time to each other, so that the respective brain responses overlap.
If you want more detailed introduction to this topic check out my paper: https://peerj.com/articles/7838/
~~~~{.julia}

println(pwd())
data, evts = loadtestdata("testcase2","../../test/");
~~~~~~~~~~~~~


~~~~
C:\Users\behinger\.julia\dev\unfold\doc\src
([0.1631473347, 0.3657296693, 0.6774393828, -0.1076610889, 0.2682471728, 0.2871913762, -0.00
77087311, 0.0299771656, -0.4588290227, -0.2621458047  …  -0.2578017734, 0.0061671647, -0.338
6812029, 0.3918692339, 0.0189132161, -0.1214769364, -0.4293630494, 0.2533207626, 0.264371720
9, 0.4062427968], 299×5 DataFrame
│ Row │ latency │ type      │ intercept │ conditionA │ conditionB │
│     │ Int64   │ String    │ Int64     │ Int64      │ Int64      │
├─────┼─────────┼───────────┼───────────┼────────────┼────────────┤
│ 1   │ 20      │ stimulus2 │ 1         │ 1          │ 0          │
│ 2   │ 40      │ stimulus2 │ 1         │ 1          │ 1          │
│ 3   │ 69      │ stimulus2 │ 1         │ 1          │ 0          │
│ 4   │ 90      │ stimulus2 │ 1         │ 0          │ 1          │
│ 5   │ 109     │ stimulus2 │ 1         │ 1          │ 0          │
│ 6   │ 128     │ stimulus2 │ 1         │ 0          │ 1          │
│ 7   │ 147     │ stimulus2 │ 1         │ 0          │ 0          │
⋮
│ 292 │ 5839    │ stimulus2 │ 1         │ 1          │ 0          │
│ 293 │ 5860    │ stimulus2 │ 1         │ 0          │ 1          │
│ 294 │ 5879    │ stimulus2 │ 1         │ 1          │ 1          │
│ 295 │ 5910    │ stimulus2 │ 1         │ 1          │ 1          │
│ 296 │ 5940    │ stimulus2 │ 1         │ 1          │ 0          │
│ 297 │ 5960    │ stimulus2 │ 1         │ 0          │ 1          │
│ 298 │ 5979    │ stimulus2 │ 1         │ 0          │ 1          │
│ 299 │ 5999    │ stimulus2 │ 1         │ 0          │ 1          │)
~~~~



~~~~{.julia}

show(first(evts,6,),allcols=true)
~~~~~~~~~~~~~


~~~~
6×5 DataFrame
│ Row │ latency │ type      │ intercept │ conditionA │ conditionB │
│     │ Int64   │ String    │ Int64     │ Int64      │ Int64      │
├─────┼─────────┼───────────┼───────────┼────────────┼────────────┤
│ 1   │ 20      │ stimulus2 │ 1         │ 1          │ 0          │
│ 2   │ 40      │ stimulus2 │ 1         │ 1          │ 1          │
│ 3   │ 69      │ stimulus2 │ 1         │ 1          │ 0          │
│ 4   │ 90      │ stimulus2 │ 1         │ 0          │ 1          │
│ 5   │ 109     │ stimulus2 │ 1         │ 1          │ 0          │
│ 6   │ 128     │ stimulus2 │ 1         │ 0          │ 1          │
~~~~





The data has little noise and the underlying signal is a pos-neg spike pattern
~~~~{.julia}

Plots.plot(range(1/50,length=300,step=1/50),data[1:300])
Plots.vline!(evts[evts.latency.<=300,:latency]./50) # show events
~~~~~~~~~~~~~


![](figures/lm_tutorial_4_1.png)\ 





## Traditional Mass Univariate Analysis
In order to demonstrate why overlap correction is important, we will first epoch the data and fit a linear model to each time point.
This is a "traditional mass-univariate analysis".
~~~~{.julia}

# we have multi channel support
data_r = reshape(data,(1,:))
# cut the data into epochs
data_epochs,times = unfold.epoch(data=data_r,tbl=evts,τ=(-0.4,0.8),sfreq=50)
~~~~~~~~~~~~~


~~~~
(Union{Missing, Float64}[missing 0.1631473347 … 0.4391097799 -0.0399733701]

Union{Missing, Float64}[1.1812906796 2.9424269262 … -0.329011023 0.0539261965]

Union{Missing, Float64}[-0.0319555568 0.7806779454 … 0.0234600983 1.4829347135]

...

Union{Missing, Float64}[1.1693038647 4.3205835742 … 0.2643717209 0.4062427968]

Union{Missing, Float64}[-0.0618761139 0.3113625732 … missing missing]

Union{Missing, Float64}[-0.6249721672 1.0476676356 … missing missing], -0.4:0.02:0.8)
~~~~





We define a formula that we want to apply to each point in time
~~~~{.julia}

f  = @formula 0~1+conditionA+conditionB # 0 as a dummy, we will combine wit data later
~~~~~~~~~~~~~


~~~~
FormulaTerm
Response:
  0
Predictors:
  1
  conditionA(unknown)
  conditionB(unknown)
~~~~





We fit the `UnfoldLinearModel` to the data
~~~~{.julia}

m,results = unfold.fit(UnfoldLinearModel,f,evts,data_epochs,times)
~~~~~~~~~~~~~


~~~~
(Unfold object
formula: 0 ~ 1 + conditionA + conditionB
Fields:	.modelinfo (potentially contains info on fitting procedure) 
	.beta extracted parameters 
	.X designmatrix (with fields .X.Xs, .X.events .X.formulas
, 183×7 DataFrame. Omitted printing of 2 columns
│ Row │ term        │ channel │ basisname       │ colname_basis │ estimate  │
│     │ Any         │ Int64   │ String          │ Float64       │ Float64   │
├─────┼─────────────┼─────────┼─────────────────┼───────────────┼───────────┤
│ 1   │ (Intercept) │ 1       │ mass-univariate │ -0.4          │ 1.24531   │
│ 2   │ (Intercept) │ 1       │ mass-univariate │ -0.38         │ 2.08041   │
│ 3   │ (Intercept) │ 1       │ mass-univariate │ -0.36         │ 1.90501   │
│ 4   │ (Intercept) │ 1       │ mass-univariate │ -0.34         │ 1.23681   │
│ 5   │ (Intercept) │ 1       │ mass-univariate │ -0.32         │ -0.104443 │
│ 6   │ (Intercept) │ 1       │ mass-univariate │ -0.3          │ -1.24538  │
│ 7   │ (Intercept) │ 1       │ mass-univariate │ -0.28         │ -2.01211  │
⋮
│ 176 │ conditionB  │ 1       │ mass-univariate │ 0.66          │ 0.157426  │
│ 177 │ conditionB  │ 1       │ mass-univariate │ 0.68          │ -0.574824 │
│ 178 │ conditionB  │ 1       │ mass-univariate │ 0.7           │ -1.01565  │
│ 179 │ conditionB  │ 1       │ mass-univariate │ 0.72          │ -1.21997  │
│ 180 │ conditionB  │ 1       │ mass-univariate │ 0.74          │ -0.889418 │
│ 181 │ conditionB  │ 1       │ mass-univariate │ 0.76          │ -0.331796 │
│ 182 │ conditionB  │ 1       │ mass-univariate │ 0.78          │ 0.118075  │
│ 183 │ conditionB  │ 1       │ mass-univariate │ 0.8           │ 0.15781   │)
~~~~




The object has the following fields
~~~~{.julia}

println(typeof(m))
m
~~~~~~~~~~~~~


~~~~
UnfoldLinearModel
Unfold object
formula: 0 ~ 1 + conditionA + conditionB
Fields:	.modelinfo (potentially contains info on fitting procedure) 
	.beta extracted parameters 
	.X designmatrix (with fields .X.Xs, .X.events .X.formulas
~~~~




Which contain the model, the original formula, the original events and returns extra a *tidy*-dataframe with the results
~~~~{.julia}

first(results,6)
~~~~~~~~~~~~~


~~~~
6×7 DataFrame. Omitted printing of 2 columns
│ Row │ term        │ channel │ basisname       │ colname_basis │ estimate  │
│     │ Any         │ Int64   │ String          │ Float64       │ Float64   │
├─────┼─────────────┼─────────┼─────────────────┼───────────────┼───────────┤
│ 1   │ (Intercept) │ 1       │ mass-univariate │ -0.4          │ 1.24531   │
│ 2   │ (Intercept) │ 1       │ mass-univariate │ -0.38         │ 2.08041   │
│ 3   │ (Intercept) │ 1       │ mass-univariate │ -0.36         │ 1.90501   │
│ 4   │ (Intercept) │ 1       │ mass-univariate │ -0.34         │ 1.23681   │
│ 5   │ (Intercept) │ 1       │ mass-univariate │ -0.32         │ -0.104443 │
│ 6   │ (Intercept) │ 1       │ mass-univariate │ -0.3          │ -1.24538  │
~~~~





We can also plot it:
~~~~{.julia}

Plots.plot(results.colname_basis,results.estimate,
        group=results.term,
        layout=1,legend=:outerbottom)
# equivalent: plot(m)
~~~~~~~~~~~~~


![](figures/lm_tutorial_10_1.png)\ 



(:colname_basis is used instead of :time, this might change. But the reason is that not all basisfunctions have a time dimension)
As can be seen a lot is going on here. As we will see later, most of the activity is due to overlap with the next event


## Basis Functions
#### HRF / BOLD
We are now ready to define a basisfunction. There are currently only two basisfunction implemented, so not much choice.
We first have a look at the BOLD-HRF basisfunction:

~~~~{.julia}

TR = 1.5
bold = hrfbasis(TR) # using default SPM parameters
eventonset = 1.3
Plots.plot(bold.kernel(eventonset))
~~~~~~~~~~~~~


![](figures/lm_tutorial_11_1.png)\ 



Classically, we would convolve this HRF function with a impulse-vector, with impulse at the event onsets
~~~~{.julia}

y = zeros(100)
y[[10,30,37,45]] .=1
y_conv = conv(y,bold.kernel(0))
Plots.plot(y_conv)
~~~~~~~~~~~~~


![](figures/lm_tutorial_12_1.png)\ 



Which one would use as a regressor against the recorded BOLD timecourse.

Note that events could fall inbetween TR (the sampling rate). Some packages subsample the time signal, but in `unfold` we can directly call the `bold.kernel` function at a given event-time, which allows for non-TR-multiples to be used.

### FIR Basis Function

Okay, let's have a look at a different basis function: The FIR basisfunction.

~~~~{.julia}

basisfunction = firbasis(τ=(-0.4,.8),sfreq=50,name="myFIRbasis")
Plots.plot(basisfunction.kernel(0))
~~~~~~~~~~~~~


![](figures/lm_tutorial_13_1.png)\ 




Not very clear, better show it in 2D
~~~~{.julia}

basisfunction.kernel(0)[1:10,1:10]
~~~~~~~~~~~~~


~~~~
10×10 SparseArrays.SparseMatrixCSC{Int64,Int64} with 19 stored entries:
  [1 ,  1]  =  1
  [2 ,  1]  =  0
  [2 ,  2]  =  1
  [3 ,  2]  =  0
  [3 ,  3]  =  1
  [4 ,  3]  =  0
  [4 ,  4]  =  1
  [5 ,  4]  =  0
  [5 ,  5]  =  1
  [6 ,  5]  =  0
  [6 ,  6]  =  1
  [7 ,  6]  =  0
  [7 ,  7]  =  1
  [8 ,  7]  =  0
  [8 ,  8]  =  1
  [9 ,  8]  =  0
  [9 ,  9]  =  1
  [10,  9]  =  0
  [10, 10]  =  1
~~~~




The FIR basisset consists of multiple basisfunctions. That is, each event will now be *timeexpanded* to multiple predictors, each with a certain time-delay to the event onset.
This allows to model any arbitrary linear overlap shape, and doesn't force us to make assumptions on the convolution kernel (like we had to do in the BOLD case)


## Timeexpanded / Deconvolved ModelFit
Remember our formula from above
~~~~{.julia}

f
~~~~~~~~~~~~~


~~~~
FormulaTerm
Response:
  0
Predictors:
  1
  conditionA(unknown)
  conditionB(unknown)
~~~~





For the left-handside we use "0" as the data is separated from the events. This will in the future allow us to fit multiple channels easily.

And fit a `UnfoldLinearModel`. Not that instead of `times` as in the mass-univariate case, we have a `BasisFunction` object now.
~~~~{.julia}

m,results = unfold.fit(UnfoldLinearModel,f,evts,data,basisfunction)
~~~~~~~~~~~~~


~~~~
(Unfold object
formula: 0 ~ ["myFIRbasis : (Intercept) : -0.4", "myFIRbasis : (Intercept) : -0.38", "myFIRb
asis : (Intercept) : -0.36000000000000004", "myFIRbasis : (Intercept) : -0.34", "myFIRbasis 
: (Intercept) : -0.32", "myFIRbasis : (Intercept) : -0.30000000000000004", "myFIRbasis : (In
tercept) : -0.28", "myFIRbasis : (Intercept) : -0.26", "myFIRbasis : (Intercept) : -0.240000
00000000002", "myFIRbasis : (Intercept) : -0.22000000000000003", "myFIRbasis : (Intercept) :
 -0.2", "myFIRbasis : (Intercept) : -0.18000000000000002", "myFIRbasis : (Intercept) : -0.16
000000000000003", "myFIRbasis : (Intercept) : -0.14", "myFIRbasis : (Intercept) : -0.12", "m
yFIRbasis : (Intercept) : -0.10000000000000003", "myFIRbasis : (Intercept) : -0.080000000000
00002", "myFIRbasis : (Intercept) : -0.06", "myFIRbasis : (Intercept) : -0.04000000000000003
6", "myFIRbasis : (Intercept) : -0.020000000000000018", "myFIRbasis : (Intercept) : 0.0", "m
yFIRbasis : (Intercept) : 0.019999999999999962", "myFIRbasis : (Intercept) : 0.0399999999999
9998", "myFIRbasis : (Intercept) : 0.06", "myFIRbasis : (Intercept) : 0.07999999999999996", 
"myFIRbasis : (Intercept) : 0.09999999999999998", "myFIRbasis : (Intercept) : 0.12", "myFIRb
asis : (Intercept) : 0.14", "myFIRbasis : (Intercept) : 0.16000000000000003", "myFIRbasis : 
(Intercept) : 0.17999999999999994", "myFIRbasis : (Intercept) : 0.19999999999999996", "myFIR
basis : (Intercept) : 0.21999999999999997", "myFIRbasis : (Intercept) : 0.24", "myFIRbasis :
 (Intercept) : 0.26", "myFIRbasis : (Intercept) : 0.28", "myFIRbasis : (Intercept) : 0.30000
000000000004", "myFIRbasis : (Intercept) : 0.31999999999999995", "myFIRbasis : (Intercept) :
 0.33999999999999997", "myFIRbasis : (Intercept) : 0.36", "myFIRbasis : (Intercept) : 0.38",
 "myFIRbasis : (Intercept) : 0.4", "myFIRbasis : (Intercept) : 0.42000000000000004", "myFIRb
asis : (Intercept) : 0.43999999999999995", "myFIRbasis : (Intercept) : 0.45999999999999996",
 "myFIRbasis : (Intercept) : 0.48", "myFIRbasis : (Intercept) : 0.5", "myFIRbasis : (Interce
pt) : 0.52", "myFIRbasis : (Intercept) : 0.54", "myFIRbasis : (Intercept) : 0.55999999999999
99", "myFIRbasis : (Intercept) : 0.58", "myFIRbasis : (Intercept) : 0.6", "myFIRbasis : (Int
ercept) : 0.62", "myFIRbasis : (Intercept) : 0.64", "myFIRbasis : (Intercept) : 0.66", "myFI
Rbasis : (Intercept) : 0.68", "myFIRbasis : (Intercept) : 0.7000000000000001", "myFIRbasis :
 (Intercept) : 0.7200000000000001", "myFIRbasis : (Intercept) : 0.7400000000000001", "myFIRb
asis : (Intercept) : 0.7599999999999999", "myFIRbasis : (Intercept) : 0.7799999999999999", "
myFIRbasis : (Intercept) : 0.7999999999999999", "myFIRbasis : conditionA : -0.4", "myFIRbasi
s : conditionA : -0.38", "myFIRbasis : conditionA : -0.36000000000000004", "myFIRbasis : con
ditionA : -0.34", "myFIRbasis : conditionA : -0.32", "myFIRbasis : conditionA : -0.300000000
00000004", "myFIRbasis : conditionA : -0.28", "myFIRbasis : conditionA : -0.26", "myFIRbasis
 : conditionA : -0.24000000000000002", "myFIRbasis : conditionA : -0.22000000000000003", "my
FIRbasis : conditionA : -0.2", "myFIRbasis : conditionA : -0.18000000000000002", "myFIRbasis
 : conditionA : -0.16000000000000003", "myFIRbasis : conditionA : -0.14", "myFIRbasis : cond
itionA : -0.12", "myFIRbasis : conditionA : -0.10000000000000003", "myFIRbasis : conditionA 
: -0.08000000000000002", "myFIRbasis : conditionA : -0.06", "myFIRbasis : conditionA : -0.04
0000000000000036", "myFIRbasis : conditionA : -0.020000000000000018", "myFIRbasis : conditio
nA : 0.0", "myFIRbasis : conditionA : 0.019999999999999962", "myFIRbasis : conditionA : 0.03
999999999999998", "myFIRbasis : conditionA : 0.06", "myFIRbasis : conditionA : 0.07999999999
999996", "myFIRbasis : conditionA : 0.09999999999999998", "myFIRbasis : conditionA : 0.12", 
"myFIRbasis : conditionA : 0.14", "myFIRbasis : conditionA : 0.16000000000000003", "myFIRbas
is : conditionA : 0.17999999999999994", "myFIRbasis : conditionA : 0.19999999999999996", "my
FIRbasis : conditionA : 0.21999999999999997", "myFIRbasis : conditionA : 0.24", "myFIRbasis 
: conditionA : 0.26", "myFIRbasis : conditionA : 0.28", "myFIRbasis : conditionA : 0.3000000
0000000004", "myFIRbasis : conditionA : 0.31999999999999995", "myFIRbasis : conditionA : 0.3
3999999999999997", "myFIRbasis : conditionA : 0.36", "myFIRbasis : conditionA : 0.38", "myFI
Rbasis : conditionA : 0.4", "myFIRbasis : conditionA : 0.42000000000000004", "myFIRbasis : c
onditionA : 0.43999999999999995", "myFIRbasis : conditionA : 0.45999999999999996", "myFIRbas
is : conditionA : 0.48", "myFIRbasis : conditionA : 0.5", "myFIRbasis : conditionA : 0.52", 
"myFIRbasis : conditionA : 0.54", "myFIRbasis : conditionA : 0.5599999999999999", "myFIRbasi
s : conditionA : 0.58", "myFIRbasis : conditionA : 0.6", "myFIRbasis : conditionA : 0.62", "
myFIRbasis : conditionA : 0.64", "myFIRbasis : conditionA : 0.66", "myFIRbasis : conditionA 
: 0.68", "myFIRbasis : conditionA : 0.7000000000000001", "myFIRbasis : conditionA : 0.720000
0000000001", "myFIRbasis : conditionA : 0.7400000000000001", "myFIRbasis : conditionA : 0.75
99999999999999", "myFIRbasis : conditionA : 0.7799999999999999", "myFIRbasis : conditionA : 
0.7999999999999999", "myFIRbasis : conditionB : -0.4", "myFIRbasis : conditionB : -0.38", "m
yFIRbasis : conditionB : -0.36000000000000004", "myFIRbasis : conditionB : -0.34", "myFIRbas
is : conditionB : -0.32", "myFIRbasis : conditionB : -0.30000000000000004", "myFIRbasis : co
nditionB : -0.28", "myFIRbasis : conditionB : -0.26", "myFIRbasis : conditionB : -0.24000000
000000002", "myFIRbasis : conditionB : -0.22000000000000003", "myFIRbasis : conditionB : -0.
2", "myFIRbasis : conditionB : -0.18000000000000002", "myFIRbasis : conditionB : -0.16000000
000000003", "myFIRbasis : conditionB : -0.14", "myFIRbasis : conditionB : -0.12", "myFIRbasi
s : conditionB : -0.10000000000000003", "myFIRbasis : conditionB : -0.08000000000000002", "m
yFIRbasis : conditionB : -0.06", "myFIRbasis : conditionB : -0.040000000000000036", "myFIRba
sis : conditionB : -0.020000000000000018", "myFIRbasis : conditionB : 0.0", "myFIRbasis : co
nditionB : 0.019999999999999962", "myFIRbasis : conditionB : 0.03999999999999998", "myFIRbas
is : conditionB : 0.06", "myFIRbasis : conditionB : 0.07999999999999996", "myFIRbasis : cond
itionB : 0.09999999999999998", "myFIRbasis : conditionB : 0.12", "myFIRbasis : conditionB : 
0.14", "myFIRbasis : conditionB : 0.16000000000000003", "myFIRbasis : conditionB : 0.1799999
9999999994", "myFIRbasis : conditionB : 0.19999999999999996", "myFIRbasis : conditionB : 0.2
1999999999999997", "myFIRbasis : conditionB : 0.24", "myFIRbasis : conditionB : 0.26", "myFI
Rbasis : conditionB : 0.28", "myFIRbasis : conditionB : 0.30000000000000004", "myFIRbasis : 
conditionB : 0.31999999999999995", "myFIRbasis : conditionB : 0.33999999999999997", "myFIRba
sis : conditionB : 0.36", "myFIRbasis : conditionB : 0.38", "myFIRbasis : conditionB : 0.4",
 "myFIRbasis : conditionB : 0.42000000000000004", "myFIRbasis : conditionB : 0.4399999999999
9995", "myFIRbasis : conditionB : 0.45999999999999996", "myFIRbasis : conditionB : 0.48", "m
yFIRbasis : conditionB : 0.5", "myFIRbasis : conditionB : 0.52", "myFIRbasis : conditionB : 
0.54", "myFIRbasis : conditionB : 0.5599999999999999", "myFIRbasis : conditionB : 0.58", "my
FIRbasis : conditionB : 0.6", "myFIRbasis : conditionB : 0.62", "myFIRbasis : conditionB : 0
.64", "myFIRbasis : conditionB : 0.66", "myFIRbasis : conditionB : 0.68", "myFIRbasis : cond
itionB : 0.7000000000000001", "myFIRbasis : conditionB : 0.7200000000000001", "myFIRbasis : 
conditionB : 0.7400000000000001", "myFIRbasis : conditionB : 0.7599999999999999", "myFIRbasi
s : conditionB : 0.7799999999999999", "myFIRbasis : conditionB : 0.7999999999999999"]

Fields:	.modelinfo (potentially contains info on fitting procedure) 
	.beta extracted parameters 
	.X designmatrix (with fields .X.Xs, .X.events .X.formulas
, 183×7 DataFrame. Omitted printing of 2 columns
│ Row │ term        │ channel │ basisname  │ colname_basis │ estimate    │
│     │ SubString…  │ Int64   │ SubString… │ Float64       │ Float64     │
├─────┼─────────────┼─────────┼────────────┼───────────────┼─────────────┤
│ 1   │ (Intercept) │ 1       │ myFIRbasis │ -0.4          │ 0.0780395   │
│ 2   │ (Intercept) │ 1       │ myFIRbasis │ -0.38         │ 0.0659579   │
│ 3   │ (Intercept) │ 1       │ myFIRbasis │ -0.36         │ -0.00965051 │
│ 4   │ (Intercept) │ 1       │ myFIRbasis │ -0.34         │ 0.193663    │
│ 5   │ (Intercept) │ 1       │ myFIRbasis │ -0.32         │ 0.0136937   │
│ 6   │ (Intercept) │ 1       │ myFIRbasis │ -0.3          │ 0.0527402   │
│ 7   │ (Intercept) │ 1       │ myFIRbasis │ -0.28         │ 0.0543905   │
⋮
│ 176 │ conditionB  │ 1       │ myFIRbasis │ 0.66          │ -0.0430643  │
│ 177 │ conditionB  │ 1       │ myFIRbasis │ 0.68          │ -0.0623931  │
│ 178 │ conditionB  │ 1       │ myFIRbasis │ 0.7           │ 0.0784401   │
│ 179 │ conditionB  │ 1       │ myFIRbasis │ 0.72          │ 0.04224     │
│ 180 │ conditionB  │ 1       │ myFIRbasis │ 0.74          │ 0.0124627   │
│ 181 │ conditionB  │ 1       │ myFIRbasis │ 0.76          │ -0.0437433  │
│ 182 │ conditionB  │ 1       │ myFIRbasis │ 0.78          │ 0.000654015 │
│ 183 │ conditionB  │ 1       │ myFIRbasis │ 0.8           │ 0.0470511   │)
~~~~



~~~~{.julia}

Plots.plot(results.colname_basis,results.estimate,
        group=results.term,
        layout=1,legend=:outerbottom)
~~~~~~~~~~~~~


![](figures/lm_tutorial_17_1.png)\ 



Cool! All overlapping activity has been removed and we recovered the simulated underlying signal.



