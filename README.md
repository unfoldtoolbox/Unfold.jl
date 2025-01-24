# [![Unfold.jl EEG toolbox](https://github.com/unfoldtoolbox/Unfold.jl/assets/10183650/3cbe57c1-e1a7-4150-817a-ce3dcc844485)](https://github.com/unfoldtoolbox/Unfold.jl)

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://unfoldtoolbox.github.io/Unfold.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://unfoldtoolbox.github.io/Unfold.jl/dev)
[![Test workflow status](https://github.com/unfoldtoolbox/Unfold.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/Unfold.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/unfoldtoolbox/Unfold.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/Unfold.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/unfoldtoolbox/Unfold.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/unfoldtoolbox/Unfold.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/unfoldtoolbox/Unfold.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/unfoldtoolbox/Unfold.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5759066.svg)](https://doi.org/10.5281/zenodo.5759066)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![All Contributors](https://img.shields.io/github/all-contributors/unfoldtoolbox/Unfold.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)

|Estimation|Visualisation|Simulation|BIDS pipeline|Decoding|Statistics|MixedModelling|
|---|---|---|---|---|---|---|
| <a href="https://github.com/unfoldtoolbox/Unfold.jl/tree/main"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623787-757575d0-aeb9-4d94-a5f8-832f13dcd2dd.png" alt="Unfold.jl Logo"></a> | <a href="https://github.com/unfoldtoolbox/UnfoldMakie.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623793-37af35a0-c99c-4374-827b-40fc37de7c2b.png" alt="UnfoldMakie.jl Logo"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldSim.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623795-328a4ccd-8860-4b13-9fb6-64d3df9e2091.png" alt="UnfoldSim.jl Logo"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldBIDS.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622460-2956ca20-9c48-4066-9e50-c5d25c50f0d1.png" alt="UnfoldBIDS.jl Logo"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldDecode.jl"><img src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277622487-802002c0-a1f2-4236-9123-562684d39dcf.png" alt="UnfoldDecode.jl Logo"></a>|<a href="https://github.com/unfoldtoolbox/UnfoldStats.jl"><img  src="https://github-production-user-asset-6210df.s3.amazonaws.com/10183650/277623799-4c8f2b5a-ea84-4ee3-82f9-01ef05b4f4c6.png" alt="UnfoldStats.jl Logo"></a>|<a href=""><img src="https://github.com/user-attachments/assets/ffb2bba6-3a30-48b7-9849-7d4e7195b297" alt="UnfoldMixedModels.jl logo"></a>|

Package (-family) to perform linear / GAM / hierarchical / deconvolution regression on biological signals.

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

### Tipp on Docs

You can read the docs online: [![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://unfoldtoolbox.github.io/Unfold.jl/stable)  - or use the `?fit`, `?effects` julia-REPL feature. To filter docs, use e.g. `?fit(::UnfoldModel)`

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
basis = firbasis(τ=(-0.3,0.5),srate=250) # for "timeexpansion" / deconvolution

# in case of [mass univariate] we need to epoch the data into trials, and a accompanying time vector
epochs::Array{Float64,3} # channel x time x epochs (n-epochs == nrows(events))
times = range(0,length=size(epochs,3),step=1/sampling_rate)
```

To fit any of the models, Unfold.jl offers a unified syntax:

| Overlap-Correction | Mixed Modelling | julia syntax |
|:---:|:---:|---|
|  |  | `fit(UnfoldModel,[Any=>(f,times)),evts,data_epoch]` |
| x |  | `fit(UnfoldModel,[Any=>(f,basis)),evts,data]` |
|  | x | `fit(UnfoldModel,[Any=>(fLMM,times)),evts,data_epoch]` |
| x | x | `fit(UnfoldModel,[Any=>(fLMM,basis)),evts,data]` |

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
| decoding  |        |        | UnfoldDecode.jl         |
| outlier-robust fits  |        |        |  [many options (but slower)](https://unfoldtoolbox.github.io/Unfold.jl/dev/HowTo/custom_solvers/#Robust-Solvers)   |
| 🐍Python support | | | [via juliacall](https://unfoldtoolbox.github.io/Unfold.jl/dev/generated/HowTo/juliacall_unfold/)|

</details>

## Contributions

Contributions are very welcome. These could be typos, bugreports, feature-requests, speed-optimization, new solvers, better code, better documentation.

### How-to Contribute

You are very welcome to raise issues and start pull requests!

### Adding Documentation

1. We recommend to write a Literate.jl document and place it in `docs/literate/FOLDER/FILENAME.jl` with `FOLDER` being `HowTo`, `Explanation`, `Tutorial` or `Reference` ([recommended reading on the 4 categories](https://documentation.divio.com/)).
2. Literate.jl converts the `.jl` file to a `.md` automatically and places it in `docs/src/generated/FOLDER/FILENAME.md`.
3. Edit [make.jl](https://github.com/unfoldtoolbox/Unfold.jl/blob/main/docs/make.jl) with a reference to `docs/src/generated/FOLDER/FILENAME.md`.

## Contributors
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/jschepers"><img src="https://avatars.githubusercontent.com/u/22366977?v=4?s=100" width="100px;" alt="Judith Schepers"/><br /><sub><b>Judith Schepers</b></sub></a><br /><a href="#bug-jschepers" title="Bug reports">🐛</a> <a href="#code-jschepers" title="Code">💻</a> <a href="#doc-jschepers" title="Documentation">📖</a> <a href="#tutorial-jschepers" title="Tutorials">✅</a> <a href="#ideas-jschepers" title="Ideas, Planning, & Feedback">🤔</a> <a href="#test-jschepers" title="Tests">⚠️</a></td>
      <td align="center" valign="top" width="14.28%"><a href="http://www.benediktehinger.de"><img src="https://avatars.githubusercontent.com/u/10183650?v=4?s=100" width="100px;" alt="Benedikt Ehinger"/><br /><sub><b>Benedikt Ehinger</b></sub></a><br /><a href="#bug-behinger" title="Bug reports">🐛</a> <a href="#code-behinger" title="Code">💻</a> <a href="#doc-behinger" title="Documentation">📖</a> <a href="#tutorial-behinger" title="Tutorials">✅</a> <a href="#ideas-behinger" title="Ideas, Planning, & Feedback">🤔</a> <a href="#test-behinger" title="Tests">⚠️</a> <a href="#infra-behinger" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#test-behinger" title="Tests">⚠️</a> <a href="#maintenance-behinger" title="Maintenance">🚧</a> <a href="#review-behinger" title="Reviewed Pull Requests">👀</a> <a href="#question-behinger" title="Answering Questions">💬</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://reneskukies.de/"><img src="https://avatars.githubusercontent.com/u/57703446?v=4?s=100" width="100px;" alt="René Skukies"/><br /><sub><b>René Skukies</b></sub></a><br /><a href="#bug-ReneSkukies" title="Bug reports">🐛</a> <a href="#doc-ReneSkukies" title="Documentation">📖</a> <a href="#tutorial-ReneSkukies" title="Tutorials">✅</a> <a href="#code-ReneSkukies" title="Code">💻</a> <a href="#ideas-ReneSkukies" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://reboreexplore.github.io/"><img src="https://avatars.githubusercontent.com/u/43548330?v=4?s=100" width="100px;" alt="Manpa Barman"/><br /><sub><b>Manpa Barman</b></sub></a><br /><a href="#infra-ReboreExplore" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://www.phillipalday.com"><img src="https://avatars.githubusercontent.com/u/1677783?v=4?s=100" width="100px;" alt="Phillip Alday"/><br /><sub><b>Phillip Alday</b></sub></a><br /><a href="#code-palday" title="Code">💻</a> <a href="#infra-palday" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a></td>
      <td align="center" valign="top" width="14.28%"><a href="http://davekleinschmidt.com"><img src="https://avatars.githubusercontent.com/u/135920?v=4?s=100" width="100px;" alt="Dave Kleinschmidt"/><br /><sub><b>Dave Kleinschmidt</b></sub></a><br /><a href="#doc-kleinschmidt" title="Documentation">📖</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/ssaket"><img src="https://avatars.githubusercontent.com/u/27828189?v=4?s=100" width="100px;" alt="Saket Saurabh"/><br /><sub><b>Saket Saurabh</b></sub></a><br /><a href="#bug-ssaket" title="Bug reports">🐛</a></td>
    </tr>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/suddha-bpn"><img src="https://avatars.githubusercontent.com/u/7974144?v=4?s=100" width="100px;" alt="suddha-bpn"/><br /><sub><b>suddha-bpn</b></sub></a><br /><a href="#bug-suddha-bpn" title="Bug reports">🐛</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/vladdez"><img src="https://avatars.githubusercontent.com/u/33777074?v=4?s=100" width="100px;" alt="Vladimir Mikheev"/><br /><sub><b>Vladimir Mikheev</b></sub></a><br /><a href="#bug-vladdez" title="Bug reports">🐛</a> <a href="#doc-vladdez" title="Documentation">📖</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/carmenamme"><img src="https://avatars.githubusercontent.com/u/100191854?v=4?s=100" width="100px;" alt="carmenamme"/><br /><sub><b>carmenamme</b></sub></a><br /><a href="#doc-carmenamme" title="Documentation">📖</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/maxvanmigem"><img src="https://avatars.githubusercontent.com/u/115144441?v=4?s=100" width="100px;" alt="Maximilien Van Migem"/><br /><sub><b>Maximilien Van Migem</b></sub></a><br /><a href="#bug-maxvanmigem" title="Bug reports">🐛</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/Till223"><img src="https://avatars.githubusercontent.com/u/29772145?v=4?s=100" width="100px;" alt="Till Prölß"/><br /><sub><b>Till Prölß</b></sub></a><br /><a href="#doc-Till223" title="Documentation">📖</a> <a href="#bug-Till223" title="Bug reports">🐛</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://allcontributors.org/docs/en/specification) specification.

Contributions of any kind welcome!

## Citation

For now, please cite

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5759066.svg)](https://doi.org/10.5281/zenodo.5759066) and/or [Ehinger & Dimigen](https://peerj.com/articles/7838/)

## Acknowledgements

This work was initially supported by the Center for Interdisciplinary Research, Bielefeld (ZiF) Cooperation Group "Statistical models for psychological and linguistic data".

Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany´s Excellence Strategy – EXC 2075 – 390740016
