
# Basis Functions

This document provides you an explanation of basisfunctions. We start with fMRI because they are very popular.

#### HRF / BOLD
We are now ready to define a basisfunction. There are currently only few basisfunction implemented.
We first have a look at the BOLD-HRF basisfunction:

```@example Main


TR = 1.5
bold = hrfbasis(TR) # using default SPM parameters
eventonset = 1.3
Plots.plot(bold.kernel(eventonset))
```



Classically, we would convolve this HRF function with a impulse-vector, with impulse at the event onsets
```@example Main


y = zeros(100)
y[[10,30,37,45]] .=1
y_conv = conv(y,bold.kernel(0))
Plots.plot(y_conv)
```

Which one would use as a regressor against the recorded BOLD timecourse.

note that events could fall inbetween TR (the sampling rate). Some packages subsample the time signal, but in `Unfold` we can directly call the `bold.kernel` function at a given event-time, which allows for non-TR-multiples to be used.

### FIR Basis Function

Okay, let's have a look at a different basis function: The FIR basisfunction.

```@example Main


basisfunction = firbasis(τ=(-0.4,.8),sfreq=50,name="myFIRbasis")
Plots.plot(basisfunction.kernel(0))
```



Not very clear, better show it in 2D
```@example Main
basisfunction.kernel(0)[1:10,1:10]
```
(all `.` are `0`'s)



The FIR basisset consists of multiple basisfunctions. That is, each event will now be *timeexpanded* to multiple predictors, each with a certain time-delay to the event onset.
This allows to model any arbitrary linear overlap shape, and doesn't force us to make assumptions on the convolution kernel (like we had to do in the BOLD case)
