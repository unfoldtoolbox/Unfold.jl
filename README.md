# ![Unfold.jl EEG toolbox](https://github.com/unfoldtoolbox/Unfold.jl/assets/10183650/3cbe57c1-e1a7-4150-817a-ce3dcc844485)

Toolbox to perform linear regression on biological signals. 

[![Docs](https://img.shields.io/badge/docs-main-blue.svg)](https://unfoldtoolbox.github.io/Unfold.jl/dev)
![semver](https://img.shields.io/badge/semantic-versioning-green)
![](https://github.com/unfoldtoolbox/Unfold.jl/workflows/CI/badge.svg)

This tool can model event related time series with mass-univariate linear (mixed) models, with optional non-linear effects and overlap correction.

This kind of modelling is also known as encoding modeling, linear deconvolution, Temporal Response Functions (TRFs), linear system identification, and probably under other names. fMRI models with HRF-basis functions are also supported.

## Discover the Unfold-Family & Friends
<details>
<summary>Click to expand</summary>
  
| | Link to Toolbox|
|---|---|
You are here | ![logo_Unfold jl_30px](https://github.com/unfoldtoolbox/Unfold.jl/assets/10183650/757575d0-aeb9-4d94-a5f8-832f13dcd2dd)
Pipeline |![[](unfoldbids.com)](https://github.com/unfoldtoolbox/Unfold.jl/assets/10183650/2956ca20-9c48-4066-9e50-c5d25c50f0d1)
Decode | ![logo_UnfoldDecode jl_30px](https://github.com/unfoldtoolbox/Unfold.jl/assets/10183650/802002c0-a1f2-4236-9123-562684d39dcf)
![logo_UnfoldStats jl_30px](https://github.com/unfoldtoolbox/Unfold.jl/assets/10183650/4c8f2b5a-ea84-4ee3-82f9-01ef05b4f4c6)
![logo_UnfoldSim jl_30px](https://github.com/unfoldtoolbox/Unfold.jl/assets/10183650/328a4ccd-8860-4b13-9fb6-64d3df9e2091)
![logo_UnfoldMakie jl_30px](https://github.com/unfoldtoolbox/Unfold.jl/assets/10183650/37af35a0-c99c-4374-827b-40fc37de7c2b)
PyMNE.jl
</details>


## Citation
For now, please cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6423476.svg)](https://doi.org/10.5281/zenodo.6423476) or [Ehinger & Dimigen](https://peerj.com/articles/7838/)

## Relation to Unfold (matlab)
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
| 🐍Python support | | | [via Pycall, link to notebook](https://github.com/unfoldtoolbox/Unfold.jl/blob/main/docs/src/HowTo/pyjulia_unfold.ipynb)|

## Install
```julia
]add Unfold
```

## Usage
Please check out [the documentatio)n](https://unfoldtoolbox.github.io/Unfold.jl/dev) for extensive tutorials, explanations...

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
