# [Mass Univariate Linear Models (no overlap correction)](@id lm_massunivariate)

In this notebook we will fit regression models to simulated EEG data. We will see that we need some type of overlap correction, as the events are close in time to each other, so that the respective brain responses overlap.
If you want more detailed introduction to this topic check out [our paper](https://peerj.com/articles/7838/).



## Setting up & loading the data
```@example Main
using DataFrames
using Unfold
using UnfoldMakie, CairoMakie # for plotting
using UnfoldSim


nothing # hide
```


## Load Data
```@example Main
data, evts = UnfoldSim.predef_eeg()
nothing # hide
```
## Inspection
The data has only little noise. The underlying signal pattern is a positive-negative-positive spike.
```@example Main
times = range(1/50, length=200, step=1/50)
f = Figure()
plot(f[1, 1], data[1:200], times)
vlines!(evts[evts.latency .<= 200, :latency] ./ 50) # show events, latency in samples!
f
```

To inspect the event dataframe we use
```@example Main
show(first(evts, 6), allcols = true)
```
Every row is an experimental event. Note that `:latency` refers to time in samples, whereas `:onset` would typically refer to seconds.


## Traditional Mass Univariate Analysis
To perform a mass univariate analysis, you must complete the following steps:

1. Split data into epochs 
2. Specify a formula 
3. Fit a linear model to each time point & channel
4. Visualize the results.


#### 1. Split data into epochs 
Initially, you have data with a duration that represents the whole experimental trial. You need to cut the data into small regular epochs related to the some event, e.g. start of fixation.

```@example Main
# we have multi channel support
data_r = reshape(data, (1,:))
# cut the data into epochs
data_epochs, times = Unfold.epoch(data = data_r, tbl = evts, τ = (-0.4, 0.8), sfreq = 50);
size(data_epochs)
```
- `τ` specifies the epoch size.
- `sfreq` - sampling rate, converts `τ` to samples.


```@example Main
typeof(data_epochs)
```
!!! note
   Partial trials could be included or ignored by the corresponding functions due to `missing` data type. Check functions such as the Julia-based `disallowmissing` and the internal `Unfold.drop_missing_epochs`.

#### 2. Specify a formula
Define a formula to be applied to each time point.

```@example Main
f = @formula 0 ~ 1 + condition + continuous # 0 as a dummy, we will combine with data later
nothing # hide
```

#### 3. Fit a linear model to each time point & channel

Fit the `UnfoldLinearModel`.
```@example Main
m = fit(UnfoldModel, f, evts, data_epochs, times); 
nothing #hide
```

Alternative way to call this model is below. This syntax allows you to fit multiple events at once. For example, replacing `Any` with `:fixation =>...` will fit this model specifically to the fixation event type.
```@example Main
m = fit(UnfoldModel, Dict(Any=>(f, times)), evts, data_epochs); 
nothing #hide
```

Inspect the fitted model:
```@example Main
m
```
Note these functions to discover the model: `design`, `designmatrix`, `modelfit` and most importantly, `coeftable`. 

!!! info
        There are of course further methods, e.g. `coef`, `ranef`, `Unfold.formula`, `modelmatrix` which might be helpful at some point, but not important now.

Using `coeftable`, we can get a *tidy* DataFrames, very useful for your further analysis.

```@example Main
first(coeftable(m), 6)
```

#### 4. Visualize the results
Tidy DataFrames are easy to visualize using e.g. AlgebraOfGraphics.jl. Function `plot_erp` from `UnfoldMakie`makes it even easier.  

```@example Main
results = coeftable(m)
plot_erp(results)
```
As you can see, there is a lot going on, even in the baseline period! This is because the signal was simulated with overlapping events. In the next tutorial you will learn how to fix this.
