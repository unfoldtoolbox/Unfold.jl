# Basis Functions
```@setup main
using CairoMakie
```

This document will give you an explanation of basis functions. We start with basis functions for fMRI because they are very popular.

#### HRF / BOLD
We want to define a basis function. There are currently only few basisfunctions implemented in Unfold.jl, but your imagination knows no borders!

We first have a look at the BOLD-HRF basisfunction aka [Blood Oxygenation Level Dependent Hemodynamic Response Function](https://en.wikipedia.org/wiki/Blood-oxygen-level-dependent_imaging):

```@example main
using Unfold, DSP

TR = 1.5 # the sampling rate
bold = hrfbasis(TR) # using default SPM parameters
eventonset = 1.3
bold_kernel = e -> Unfold.kernel(bold, e)
lines(bold_kernel(eventonset)[:,1]) # returns a matrix, thus [:, 1]
```
This is the shape that is assumed to reflect the activity for an event. Generally, we would like to know how much to scale this response shape per condition, e.g. in `condA` we might scale it by 0.7, in `condB` by 1.2.

But let's start at the beginning and first simulate an fMRI signal. Then you will also appreciate why we need to deconvolve it later.

### Convolving a response shape to get a "recorded" fMRI signal
We start by convolving this HRF function with an impulse vector at event onsets.

```@example main
y = zeros(100) # signal length = 100
y[[10, 30, 45]] .= 0.7 # 3 events at given for condition A
y[[37]] .= 1.2 # 1 events at given for condition B

y_conv = conv(y, bold_kernel(0)) # convolve!
lines(y_conv[:,1])
```
Next, we would add some noise:

```@example main
using Random
y_conv += randn(size(y_conv))
lines(y_conv[:,1])
```
🎉 - we did it, we simulated fMRI data.

Now you can see that the conditions overlap in time. To get back to the original amplitude values, we need to specify a basis function and use Unfold to deconvolve the signals.

!!! note 
    Events can fall between TR (the sampling rate). Some packages subsample the time signal, but in `Unfold` we can call the `bold.kernel` function directly at a given event time, which allows us to use non-TR multiples.


### FIR Basis Function

Okay, let's have a look at a different basis function: The FIR basisfunction. FIR stands for [Finite-Impulse-Response](https://en.wikipedia.org/wiki/Finite_impulse_response) and is a term taken from the filtering literature.

```@example main
using Unfold #hide

basisfunction = firbasis(τ=(-0.4,.8), sfreq=50, name="myFIRbasis")
fir_kernel = e -> Unfold.kernel(basisfunction, e)
m = fir_kernel(0)
f = Figure()
f[1,1] = Axis(f)
for col = 1:size(m, 2)
    lines!(m[:,col])
end
current_figure()
```

The first thing to notice is that it is not a single basisfunction, but a set of basisfunctions. So every condition is explained by several basis functions!

To make it clear better show it in 2D:

```@example main
fir_kernel(0)[1:10,1:10]
```
(all `.` are `0`'s)

The FIR basis set consists of multiple basis functions. That is, each event is now *time-expanded* to multiple predictors, each with a certain time delay to the event onset.
This allows us to model any linear overlap shape, and doesn't force us to make assumptions about the convolution kernel, as we had to do in the BOLD case.
