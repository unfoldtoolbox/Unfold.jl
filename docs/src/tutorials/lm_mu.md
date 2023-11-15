# [Mass Univariate Linear Models (no overlap correction)](@id lm_massunivariate)

In this notebook we will fit regression models to (simulated) EEG data. We will see that we need some type of overlap correction, as the events are close in time to each other, so that the respective brain responses overlap.
If you want more detailed introduction to this topic check out [our paper](https://peerj.com/articles/7838/).



## Setting up & loading the data
```@example Main
using DataFrames
using Unfold
using UnfoldMakie,CairoMakie # for plotting
using UnfoldSim


nothing # hide
```


## Load Data
```@example Main
data, evts = UnfoldSim.predef_eeg()
nothing # hide
```
## Inspection
The data has little noise and the underlying signal is a pos-neg-pos spike pattern
```@example Main
times = range(1/50,length=200,step=1/50)
plot(times,data[1:200])
vlines!(current_axis(),evts[evts.latency.<=200,:latency]./50) # show events, latency in samples!
```

To inspect the event dataframe we use
```@example Main
show(first(evts,6,),allcols=true)
```
Every row is an event. Note that `:latency` is commonly the timestamp in samples, whereas `:onset` would typically refer to seconds.


## Traditional Mass Univariate Analysis
For a mass univariate analysis, will 
1. epoch the data
2. specify a formula 
3. fit a linear model to each time point & channel
4. visualize the results.


#### 1. Epoch the data
We have to cut the data into small regular trials.

```@example Main
# we have multi channel support
data_r = reshape(data,(1,:))
# cut the data into epochs
data_epochs,times = Unfold.epoch(data=data_r,tbl=evts,τ=(-0.4,0.8),sfreq=50);
size(data_epochs)
```

`τ` specifies the epoch size. To convert `τ` to samples, we need to specify the sampling rate.


```@example Main
typeof(data_epochs)
```
!!! note
    In julia, `missing` is supported throughout the ecosystem. Thus, we can have partial trials and they will be incorporated / ignored at the respective functions. Helpful functions are the julia-base `disallowmissing` and the internal `Unfold.dropMissingEpochs` functions



#### 2. Specify a formula
We define a formula that we want to apply to each point in time
```@example Main
f  = @formula 0~1+condition+continuous # 0 as a dummy, we will combine wit data later
nothing # hide
```

#### 3. fit a linear model to each time-point & channel

We fit the `UnfoldLinearModel`.
```@example Main
m = fit(UnfoldModel,f,evts,data_epochs,times); 
nothing #hide
```

Alternatively, we could also call it using:
```@example Main
m = fit(UnfoldModel,Dict(Any=>(f,times)),evts,data_epochs); 
nothing #hide
```
A syntax which allows for fitting multiple events at the same time (by replacing Any with e.g. `:fixation =>...` it will fit that model specifically to the fixation event type)

We can inspect the object
```@example Main
m
```
And see that there are several helpful functions to recover `design`, `designmatrix`, raw-`modelfit` and most importantly, `coeftable`. 

!!! info
        There are of course further methods, e.g. `coef`, `ranef`, `Unfold.formula`, `modelmatrix` which might be helpful at some point, but not important now

Using `coeftable`, we can get a *tidy*-dataframe, very useful for your further analysis.

```@example Main
first(coeftable(m),6)
```

#### 4. Visualize the results
Tidy-Dataframes make them easy to visualize using e.g. AlgebraOfGraphics.jl. We simpliefied this even more, by implementing a wrapper in `UnfoldMakie`
```@example Main
results = coeftable(m)
plot_erp(results)
```
As you can see here, a lot is going on, even in the baseline-period! This is because the signal was simulated with overlapping events. Head over to the next tutorial to find out how to remedy this.
