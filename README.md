# Unfold.jl

<img align="right" width="33%" src="https://www.unfoldtoolbox.org/_images/unfold_800x377.png">



Toolbox to perform linear regression on biological signals. 


[![Docs](https://img.shields.io/badge/docs-main-blue.svg)](https://unfoldtoolbox.github.io/Unfold.jl/dev)
![semver](https://img.shields.io/badge/semantic-versioning-green)
![](https://github.com/unfoldtoolbox/Unfold.jl/workflows/CI/badge.svg)

This tool can model event related time series with mass-univariate linear (mixed) models, with optional non-linear effects and overlap correction.

This kind of modelling is also known as encoding modeling, linear deconvolution, Temporal Response Functions (TRFs), linear system identification, and probably under other names. fMRI models with HRF-basis functions are also supported.

## Citation
For now, please cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6423476.svg)](https://doi.org/10.5281/zenodo.6423476) or [Ehinger & Dimigen](https://peerj.com/articles/7838/)

## Relation to Unfold (matlab)
I regularly use Unfold.jl in my research, but there are likely more hard corners there compared to the matlab version (but julia is more flexible and has some unique features)

| Feature                 | Unfold | unmixed | Unfold.jl |
|-------------------------|--------|---------|-----------|
| overlap correction      | x      | x       | x         |
| non-linear splines      | x      | x       | x         |
| plotting tools          | x      |         | UnfoldMakie.jl  |
| sanity checks           | x      |         | x         |
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
Please check out [the documentation](https://unfoldtoolbox.github.io/Unfold.jl/dev) for extensive tutorials, explanations...

Here a quick overview what to expect.

#### What you need
```julia
events::DataFrame

# formula with or without random effects
f = @formula 0~1+condA
fLMM = @formula 0~1+condA+(1|subject) + (1|item)

# in case of [overlap-correction] we need continuous data plus per-eventtype one basisfunction (typically firbasis)
data::Array{Float64,2}
basis = firbasis(τ=(-0.3,0.5),srate=250)

# in case of [mass univariate] we need to epoch the data into trials, and a accompanying time vector
epochs::Array{Float64,3} # channel x time x epochs (n-epochs == nrows(events))
times = range(0,length=size(epochs,3),step=1/sampling_rate)
```

To fit any of the models, Unfold.jl offers a unified syntax:

| Overlap-Correction | Mixed Modelling | julia syntax |
|:---:|:---:|---|
|  |  | `fit(UnfoldModel,Dict(Any=>(f,times)),evts,data_epoch)` |
| x |  | `fit(UnfoldModel,Dict(Any=>(f,basis)),evts,data)` |
|  | x | `fit(UnfoldModel,Dict(Any=>(fLMM,times)),evts,data_epoch)` |
| x | x | `fit(UnfoldModel,Dict(Any=>(fLMM,basis)),evts,data)` | 


## Documentation
Many functions have documentation from the Julia REPL by typing e.g. `julia>?Unfold.fit`

For tutorials see [the documentation](https://unfoldtoolbox.github.io/Unfold.jl/dev/)

## Contributions
Contributions are very welcome. These could be typos, bugreports, feature-requests, speed-optimization, new solvers, better code, better documentation.

#### How-to contribute
You are very welcome to raise issues and start pull requests!

#### Adding Documentation
1. We recommend to write a Literate.jl document and place it in `docs/_literate/FOLDER/FILENAME.jl` with `FOLDER` being `HowTo`, `Explanation`, `Tutorial` or `Reference` ([recommended reading on the 4 categories](https://documentation.divio.com/)).
2. Literate.jl converts the `.jl` file to a `.md` automatically and places it in `doc/src/_literate/FILENAME.jl`.
3. Edit [make.jl](https://github.com/unfoldtoolbox/Unfold.jl/blob/main/docs/make.jl) with a reference to `doc/src/_literate/FILENAME.jl`



## Contributors (alphabetically)
- **Phillip Alday**
- **Benedikt Ehinger**
- **Dave Kleinschmidt**
- **Judith Schepers**
- **Felix Schröder**
- **René Skukies**


## Acknowledgements
This work was supported by the Center for Interdisciplinary Research, Bielefeld (ZiF) Cooperation Group "Statistical models for psychological and linguistic data".

Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany´s Excellence Strategy – EXC 2075 – 390740016
