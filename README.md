# Unfold.jl

**Beta** Toolbox to perform linear regression on biological signals. 

![](https://github.com/unfoldtoolbox/Unfold.jl/workflows/CI/badge.svg)

[link to tutorials / documentation](https://unfoldtoolbox.github.io/Unfold.jl/dev/)

This tool combines mass-univariate linear (mixed) models with overlap correction.

This kind of overlap correction is also known as encoding modeling, linear deconvolution, Temporal Response Functions (TRFs) and probably under other names. fMRI models with HRF-basis functions are also supported.

## Relation to Unfold (matlab)
The matlab toolbox is recommended for research work. It is richer in features, better documented and tested.

The julia toolbox is a type of research-playground, but offers LinearMixedModel support.


| Feature                 | Unfold | unmixed | Unfold.jl |
|-------------------------|--------|---------|-----------|
| overlap correction      | x      | x       | x         |
| non-linear splines      | x      | x       | x         |
| plotting tools          | x      |         | UnfoldMakie.jl - beta        |
| sanity checks           | x      |         |           |
| tutorials               | x      |         | x       |
| speed                   | x      |         | x         |
| unittests               | x      |         | x         |
| HRF (fMRI) basis        |        |         | x         |
| mix different basisfunctions      |        |         | x         |
| different timewindows per event   |        |         | x         |
| mixed models            |        | x       | x         |
| item & subject effects  |        | x       | x         |
| decoding  |        |        | back2back regression         |

## Install
```julia
]add Unfold
```

## Usage
For a quickstart:

```julia
f = @formula 0~1+condA
fLMM = @formula 0~1+condA+(1|subject) + (1|item)
events::DataFrame
data::Array{Float64,2}
epochs::Array{Float64,3} # channel x time x epochs (n-epochs == nrows(events))
times = range(0,length=size(epochs,3),step=1/sampling_rate)

basisfunction::Unfold.BasisFunction
basis = firbasis(τ=(-0.3,0.5),srate=250)
```


1. Timeexpansion **No**, Mixed **No**  : `fit(UnfoldModel,Dict(Any=>(f,times)),evts,data_epoch)`
1. Timeexpansion **Yes**, Mixed **No** : `fit(UnfoldModel,Dict(Any=>(f,basis)),evts,data)`
1. Timeexpansion **No**, Mixed **Yes** : `fit(UnfoldModel,Dict(Any=>(fLMM,times)),evts,data_epoch)`
1. Timeexpansion **Yes**, Mixed **Yes**: `fit(UnfoldModel,Dict(Any=>(fLMM,basis)),evts,data)`


## Documentation
Most functions have documentation, e.g. `?Unfold.fit`

Tutorials see [the documentation](https://unfoldtoolbox.github.io/Unfold.jl/dev/)



## Contributors (alphabetically)
- **Phillip Alday**
- **Benedikt Ehinger**
- **Dave Kleinschmidt**
- **Judith Schepers**
- **Felix Schröder**
- **René Skukies**


## Acknowledgements
This work was supported by the Center for Interdisciplinary Research, Bielefeld (ZiF) Cooperation Group "Statistical models for psychological and linguistic data".

## Research notice
Please note that this repository is participating in a study into sustainability
 of open source projects. Data will be gathered about this repository for
 approximately the next 12 months, starting from June 2021.

Data collected will include number of contributors, number of PRs, time taken to
 close/merge these PRs, and issues closed.

For more information, please visit
[our informational page](https://sustainable-open-science-and-software.github.io/) or download our [participant information sheet](https://sustainable-open-science-and-software.github.io/assets/PIS_sustainable_software.pdf).
