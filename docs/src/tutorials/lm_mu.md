# [Mass Univariate Linear Models (no overlap correction)](@id lm_massunivariate)

In this notebook we will fit regression models to simulated EEG data. We will see that we need some type of overlap correction, as the events are close in time to each other, so that the respective brain responses overlap.
If you want more detailed introduction to this topic check out [our paper](https://peerj.com/articles/7838/).

## Setting up & loading the data

```@example Main
using DataFrames
using Unfold
using UnfoldMakie, CairoMakie # for plotting
using UnfoldSim
using DisplayAs # hide

nothing # hide
```

## Load Data

We'll start with some predefined simulated continuos EEG data. We have 2000 events, 1 channel and one condition with two levels

```@example Main
data, evts = UnfoldSim.predef_eeg()
nothing # hide
```

## Inspection

The data has only little noise. The underlying signal pattern is a positive-negative-positive spike.

```@example Main
times_cont = range(0,length=200,step=1/100) # we simulated with 100hz for 0.5 seconds

f,ax,h = plot(times_cont,data[1:200])
vlines!(evts[evts.latency .<= 200, :latency] ./ 100;color=:black) # show events, latency in samples!
ax.xlabel = "time [s]"
ax.ylabel = "voltage [µV]"
f
```

To inspect the event dataframe we use

```@example Main
show(first(evts, 6), allcols = true)
```

Every row is an experimental event. Note that `:latency` refers to time in samples, (in BIDS-specification,  `:onset` would typically refer to seconds).

## Traditional Mass Univariate Analysis

To perform a mass univariate analysis, you must complete the following steps:

1. Split data into epochs
2. Specify a formula
3. Fit a linear model to each time point & channel
4. Visualize the results.

### 1. Split data into epochs

Initially, you have data with a duration that represents the whole experimental trial. You need to cut the data into small regular epochs related to the some event, e.g. start of fixation.

```@example Main
# Unfold supports multi-channel, so we could provide matrix ch x time, which we can create like this from a vector:
data_r = reshape(data, (1,:))
# cut the data into epochs
data_epochs, times = Unfold.epoch(data = data, tbl = evts, τ = (-0.4, 0.8), sfreq = 100); # channel x timesteps x trials
size(data_epochs)
```

- `τ` specifies the epoch size.
- `sfreq` - sampling rate, converts `τ` to samples.

```@example Main
typeof(data_epochs)
```

!!! note
    In julia, `missing` is supported throughout the ecosystem. Thus, we can have partial trials and they will be incorporated / ignored at the respective functions. Helpful functions are the julia-base `disallowmissing` and the internal `Unfold.drop_missing_epochs` functions

### 2. Specify a formula

Define a formula to be applied to each time point (and each channel) relative to the event. `condition` and `continuous` are the names of the event-describing columns in `evts` that we want to use for modelling.

```@example Main
f = @formula 0 ~ 1 + condition + continuous # note the formulas left side is `0 ~ ` for technical reasons`
nothing # hide
```

### 3. Fit a linear model to each time point & channel

Fit the "`UnfoldModel`" (the `fit` syntax is used throughout the Julia ecosystem, with the first element indicating what kind of model to fit)

```@example Main
m = fit(UnfoldModel, f, evts, data_epochs, times);
nothing #hide
```

Alternative way to call this model is below. This syntax allows you to fit multiple events at once. For example, replacing `Any` with `:fixation =>...` will fit this model specifically to the fixation event type.

```@example Main
m = fit(UnfoldModel, [Any=>(f, times)], evts, data_epochs);
nothing #hide
```

Inspect the fitted model:

```@example Main
m
m|> DisplayAs.withcontext(:is_pluto=>true) # hide
```

Note these functions to discover the model: `design`, `designmatrix`, `modelfit` and most importantly, `coeftable`.

!!! info
        There are of course further methods, e.g. `coef`, `ranef`, `Unfold.formula`, `modelmatrix` which might be helpful at some point, but not important now.

Using `coeftable`, we can get a *tidy* DataFrames, very useful for your further analysis.

```@example Main
first(coeftable(m), 6)
```

### 4. Visualize the results

Tidy DataFrames are easy to visualize using e.g. AlgebraOfGraphics.jl. Function `plot_erp` from `UnfoldMakie`makes it even easier.

```@example Main
results = coeftable(m)
plot_erp(results)
```

As you can see, there is a lot going on, even in the baseline period! This is because the signal was simulated with overlapping events. In the next tutorial you will learn how to fix this.
