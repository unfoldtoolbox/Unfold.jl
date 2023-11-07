# [![Unfold.jl EEG toolbox](https://github.com/unfoldtoolbox/Unfold.jl/assets/10183650/3cbe57c1-e1a7-4150-817a-ce3dcc844485)](https://github.com/unfoldtoolbox/Unfold.jl)

[![Docs][Doc-img]][Doc-url] ![semver][semver-img] [![Build Status][build-img]][build-url]

[Doc-img]: https://img.shields.io/badge/docs-main-blue.svg
[Doc-url]: https://unfoldtoolbox.github.io/Unfold.jl/dev
[semver-img]: https://img.shields.io/badge/semantic-versioning-green
[build-img]: https://github.com/unfoldtoolbox/UnfoldSim.jl/workflows/CI/badge.svg
[build-url]: https://github.com/unfoldtoolbox/UnfoldSim.jl/workflows/CI.yml

|Estimation|Visualisation|Simulation|BIDS pipeline|Decoding|Statistics|
|---|---|---|---|---|---|
| <a href="https://github.com/unfoldtoolbox/Unfold.jl/tree/main"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623787-757575d0-aeb9-4d94-a5f8-832f13dcd2dd.png"></a> | <a href="https://github.com/unfoldtoolbox/UnfoldMakie.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623793-37af35a0-c99c-4374-827b-40fc37de7c2b.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldSim.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623795-328a4ccd-8860-4b13-9fb6-64d3df9e2091.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldBIDS.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622460-2956ca20-9c48-4066-9e50-c5d25c50f0d1.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldDecode.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622487-802002c0-a1f2-4236-9123-562684d39dcf.png"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldStats.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623799-4c8f2b5a-ea84-4ee3-82f9-01ef05b4f4c6.png"></a>|

Toolbox to perform linear / GAM / hierarchical / deconvolution regression on biological signals.

This kind of modelling is also known as encoding modeling, linear deconvolution, Temporal Response Functions (TRFs), linear system identification, and probably under other names. fMRI models with HRF-basis functions and pupil-dilation bases are also supported.

## Getting started

### 🐍Python User?
We clearly recommend Julia 😉 - but [Python users can use juliacall/Unfold directly from python!](https://unfoldtoolbox.github.io/Unfold.jl/dev/generated/HowTo/juliacall_unfold/)

### Julia installation
<details>
<summary>Click to expand</summary>

The recommended way to install julia is [juliaup](https://github.com/JuliaLang/juliaup).
It allows you to, e.g., easily update Julia at a later point, but also test out alpha/beta versions etc.

TL:DR; If you dont want to read the explicit instructions, just copy the following command

#### Windows

AppStore -> JuliaUp,  or `winget install julia -s msstore` in CMD

#### Mac & Linux

`curl -fsSL https://install.julialang.org | sh` in any shell
</details>

### Unfold.jl installation

```julia
using Pkg
Pkg.add("Unfold")
```

## Usage

Please check out [the documentation](https://unfoldtoolbox.github.io/Unfold.jl/dev) for extensive tutorials, explanations and more!

Here is a quick overview on what to expect.

### What you need

```julia
using Unfold

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

## Comparison to Unfold (matlab)
<details>
<summary>Click to expand</summary>

The matlab version is still maintained, but active development happens in Julia.

| Feature                 | Unfold | unmixed (defunct) | Unfold.jl |
|-------------------------|--------|---------|-----------|
| overlap correction      | x      | x       | x         |
| non-linear splines      | x      | x       | x         |
| speed |       |  🐌      | ⚡ 2-100x        |
| GPU support | | | 🚀|
| plotting tools          | x      |         | [UnfoldMakie.jl](https://unfoldtoolbox.github.io/UnfoldMakie.jl/dev/)  |
| Interactive plotting  |       |         | stay tuned - coming soon! |
| simulation tools          | x      |         | [UnfoldSim.jl](https://unfoldtoolbox.github.io/UnfoldSim.jl)  |
| BIDS support          | x      |         | alpha: [UnfoldBIDS.jl](https://github.com/ReneSkukies/UnfoldBIDS.jl/))  |
| sanity checks           | x      |         | x         |
| tutorials               | x      |         | x       |
| unittests               | x      |         | x         |
| Alternative bases e.g. HRF (fMRI)        |        |         | x         |
| mix different basisfunctions      |        |         | x         |
| different timewindows per event   |        |         | x         |
| mixed models            |        | x       | x         |
| item & subject effects  |        | (x)       | x         |
| decoding  |        |        | back2back regression         |
| outlier-robust fits  |        |        |  [many options (but slower)](https://unfoldtoolbox.github.io/Unfold.jl/dev/HowTo/custom_solvers/#Robust-Solvers)   |
| 🐍Python support | | | [via juliacall](https://unfoldtoolbox.github.io/Unfold.jl/dev/generated/HowTo/pyjulia_unfold/)|
</details>

## Contributions

Contributions are very welcome. These could be typos, bugreports, feature-requests, speed-optimization, new solvers, better code, better documentation.

### How-to Contribute

You are very welcome to raise issues and start pull requests!

### Adding Documentation

1. We recommend to write a Literate.jl document and place it in `docs/literate/FOLDER/FILENAME.jl` with `FOLDER` being `HowTo`, `Explanation`, `Tutorial` or `Reference` ([recommended reading on the 4 categories](https://documentation.divio.com/)).
2. Literate.jl converts the `.jl` file to a `.md` automatically and places it in `docs/src/generated/FOLDER/FILENAME.md`.
3. Edit [make.jl](https://github.com/unfoldtoolbox/Unfold.jl/blob/main/docs/make.jl) with a reference to `docs/src/generated/FOLDER/FILENAME.md`.

## Contributors (alphabetically)

- **Phillip Alday**
- **Benedikt Ehinger**
- **Dave Kleinschmidt**
- **Judith Schepers**
- **Felix Schröder**
- **René Skukies**

## Citation

For now, please cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6423476.svg)](https://doi.org/10.5281/zenodo.6423476) or [Ehinger & Dimigen](https://peerj.com/articles/7838/)

## Acknowledgements

This work was initially supported by the Center for Interdisciplinary Research, Bielefeld (ZiF) Cooperation Group "Statistical models for psychological and linguistic data".

Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany´s Excellence Strategy – EXC 2075 – 390740016
