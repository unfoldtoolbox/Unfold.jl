# Unfold.jl

**Beta** Toolbox to perform linear regression on biological signals. ![](https://github.com/unfoldtoolbox/Unfold.jl/workflows/CI/badge.svg)

This tool combines mass-univariate linear (mixed) models with overlap correction.

This kind of overlap correction is also known as encoding modeling, linear deconvolution, Temporal Response Functions (TRFs) and probably under other names. fMRI models with HRF-basis functions are also supported.

## Relation to Unfold (matlab)
The matlab toolbox is recommended for research work. It is richer in features, better documented and tested.

The julia toolbox is a type of research-playground, but offers LinearMixedModel support.


| Feature                 | Unfold | unmixed | Unfold.jl |
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
Pkg.add(url = "https://github.com/unfoldtoolbox/Unfold.jl")
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

```julia
f = @formula 0~1+condA
events::DataFrame
data::Array{Float64,2}
epochs::Array{Float64,3} # channel x time x epochs (n-epochs == nrows(events))
times = range(0,length=size(epochs,3),step=1/sampling_rate)

basisfunction::Unfold.BasisFunction
basis = firbasis(τ=(-0.3,0.5),srate=250)
```

1. Timeexpansion **No**, Mixed **No**  : `fit(UnfoldLinearModel,formula,events,epochs,times)`
1. Timeexpansion **No**, Mixed **Yes** : `fit(UnfoldLinearMixedModel,formula,events,epochs,times)`
1. Timeexpansion **Yes**, Mixed **No** : `fit(UnfoldLinearModel,Dict("eventname"=>(formula,basisfunction)),events,data)`
1. Timeexpansion **Yes**, Mixed **Yes**: `fit(UnfoldLinearMixedModel,Dict("eventname"=>(formula,basisfunction),"event2"=>(formula2,basis2)),events,data)`


## Documentation
Most functions have documentation, e.g. `?Unfold.fit`

Tutorials see `doc/lmm_tutorial.html` & `doc/lm_tutorial.html` - more to come. Contributions very welcome!


## Contributors (alphabetically)
- **Phillip Alday**
- **Benedikt Ehinger**
- **Dave Kleinschmidt**

## Acknowledgements
This work was supported by the Center for Interdisciplinary Research, Bielefeld (ZiF) Cooperation Group "Statistical models for psychological and linguistic data".
