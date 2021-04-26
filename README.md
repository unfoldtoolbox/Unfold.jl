# unfold.jl

**Beta** Toolbox to perform linear regression on biological signals. ![](https://github.com/unfoldtoolbox/unfold.jl/workflows/CI/badge.svg)

This tool combines mass-univariate linear (mixed) models with overlap correction.

This kind of overlap correction is also known as encoding modeling, linear deconvolution, Temporal Response Functions (TRFs) and probably under other names. Typical fMRI models with HRF-basis functions are also supported.

## Relation to unfold (matlab)
The matlab toolbox is recommended for research work. It is richer in features, better documented and tested.

The julia toolbox is a type of playground and aspires to combine unfold & unmixed in a single toolbox.


| Feature                 | unfold | unmixed | unfold.jl |
|-------------------------|--------|---------|-----------|
| overlap correction      | x      | x       | x         |
| non-linear splines      | x      | x       | x         |
| plotting tools          | x      |         |           |
| sanity checks           | x      |         |           |
| tutorials               | x      |         |           |
| speed                   | x      |         | x         |
| unittests               | x      |         | x         |
| HRF (fMRI) basis        |        |         | x         |
| mix different basisfunctions      |        |         | x         |
| different timewindows per event   |        |         | x         |
| mixed models            |        | x       | x         |
| item & subject effects  |        | x       | x         |

## Install
```julia
using Pkg;
Pkg.add(url = "https://github.com/mclements/Splines2.jl")
Pkg.add(url = "https://github.com/unfoldtoolbox/unfold.jl")
```

For some of the testing functionality in the `test/` path, you will also need

```julia
 Pkg.add("Makie") # use CairoMakie if in a headless environment
 Pkg.add("StatsMakie")
 Pkg.add("MAT")
 Pkg.add("HDF5")
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


## Contributors (alphabetically)
**Phillip Alday**
**Benedikt Ehinger**
**Dave Kleinschmidt**

## Acknowledgements
This work was supported by the Center for Interdisciplinary Research, Bielefeld (ZiF) Cooperation Group "Statistical models for psychological and linguistic data".
