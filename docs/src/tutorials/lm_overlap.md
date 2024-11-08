# [Linear Model with Overlap Correction](@id lm_overlap)

!!! note 
      We recommend you briefly go over the mass-univariate linear modelling tutorial


In this notebook we will fit regression models to (simulated) EEG data. We will see that we need some type of overlap correction, as the events are close in time to each other, so that the respective brain responses overlap.
If you want more detailed introduction to this topic check out [our paper](https://peerj.com/articles/7838/).

## Setting up & loading the data
```@example Main
using Unfold
using UnfoldSim
using UnfoldMakie,CairoMakie
using DataFrames
using DisplayAs # hide

data, evts = UnfoldSim.predef_eeg()
nothing # hide
```

## Overlap Correction
For an overlap correction analysis we will do one additional step: define a temporal basisfunction. The steps are as following:

1. specify a temporal basisfunction
2. specify a formula
3. fit a linear model for each channel (one for all timepoints!)
4. visualize the results.



## Timeexpanded / Deconvolved ModelFit
#### 1. specify a temporal basisfunction
By default, we would want to use a FIR basisfunction. See [Basis Functions](@ref) for more details.
```@example Main
basisfunction = firbasis(τ=(-0.4,.8),sfreq=100)
nothing #hide
```


#### 2. specify a formula
We specify the same formula as before
```@example Main
f  = @formula 0~1+condition+continuous
nothing #hide
```



#### 3. fit the linear model

The formula and basisfunction is not enough on their own. We also need to specify which event and which formula matches - this is important in cases where there are multiple events with different formulas
```@example Main
bf_vec = [Any=>(f,basisfunction)]
bf_vec|> DisplayAs.withcontext(:is_pluto=>true) # hide
```
!!! note
      The `Any` means to use all rows in `evts`. In case you have multiple events, you'd want to specify multiple basisfunctions e.g. 
      ```
      bfDict = ["stimulus"=>(f1,basisfunction1),
                    "response"=>(f2,basisfunction2)]
      ```

      You likely have to specify a further argument to `fit`: `eventcolumn="type"` with `type` being the column in `evts` that codes for the event (stimulus / response in this case)


Now we are ready to fit a `UnfoldLinearModel`. Not that instead of `times` as in the mass-univariate case, we have to provide the `BasisFunction` type now.

```@example Main
m = fit(UnfoldModel,bf_vec,evts,data); 
nothing #hide
```

#### 4. Visualize the model
Similarly to the previous tutorial, we can visualize the model
```@example Main
results = coeftable(m)
plot_erp(results)
```
Cool! All overlapping activity has been removed and we recovered the simulated underlying signal.



