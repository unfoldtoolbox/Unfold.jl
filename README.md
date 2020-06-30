# unfold.jl

**Beta** Toolbox to perform linear regression on biological signals.

This tool combines mass-univariate linear (mixed) models with overlap correction.

This kind of overlap correction is also known as encoding modeling, linear deconvolution, Temporal Response Functions (TRFs) and probably under other names. Typical fMRI models with HRF-basis functions are also supported.

![](https://github.com/unfoldtoolbox/unfold.jl/workflows/CI/badge.svg)

## Install
```
Pkg.add("https://github.com/unfoldtoolbox/unfold.jl")
Pkg.add("https://github.com/JuliaStats/MixedModels.jl") #to install latest mixed model 3.0.0 - will not be necessary once released
using unfold

```

## Usage
For a quickstart:

1. Timeexpansion **No**, Mixed **No**  : `fit(UnfoldLinearModel,f,evts,data_epoch,times)`
1. Timeexpansion **Yes**, Mixed **No** : `fit(UnfoldLinearModel,f,evts,data,basisfunction)`
1. Timeexpansion **No**, Mixed **Yes** : `fit(UnfoldLinearMixedModel,f,evts,data_epoch,times)`
1. Timeexpansion **Yes**, Mixed **Yes**: `fit(UnfoldLinearMixedModel,f,evts,data,basisfunction)`

With **formula** e.g. `@formula 0~1+condA`, **evts** a `DataFrame` with events, **data** an `Array{Number,}` and  **basisfunction**: `firbasis(τ=[-0.1,0.5],srate=50)` or **times** `range(-0.1,0.5,step=1/50)`


## Documentation
Still being written. Tutorials see `doc/lmm_tutorial.html` & `doc/lm_tutorial.html`


## Acknowledgements
Thanks to **Dave Kleinschmidt** for discussing the formula interface and thanks to him & **Phillip Alday** who answered all my (potentially naive) Julia questions.
This work was supported by the Center for Interdisciplinary Research, Bielefeld (ZiF) Cooperation Group "Statistical models for psychological and linguistic data".
