# unfold.jl
**Alpha** Toolbox to perform linear regression on biological signals.

This tool combines mass-univariate linear (mixed) models with overlap correction.

This kind of overlap correction is also known as encoding modeling, linear deconvolution, Temporal Response Functions (TRFs) and probably under other names.


## Install
```
Pkd.add("https://github.com/unfoldtoolbox/unfold.jl")
using unfold
```
Currently nothing is exported, thus all function calls need to be: `ùnfold.functioncall`

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
Thanks to **Dave Kleinschmidt** for discussing the formula interface and thanks to him & **Phillip Alday** who answered all my naive Julia questions.
This work was supported by the Center for Interdisciplinary Research, Bielefeld (ZiF) Cooperation Group "Statistical models for psychological and linguistic data".
