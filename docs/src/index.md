```@meta
CurrentModule = Unfold
```

# Unfold.jl Documentation

Welcome to [Unfold.jl](https://unfoldtoolbox.github.io/UnfoldDocs/Unfold.jl/stable/): a Julia Package for regression and event-based time series analysis, with a focus on regression ERPs for EEG analysis. The modular approach allows for easy modification to other context, like iEEG, pupil dilation, fMRI etc. - while maintaining the speed of Julia!

```@raw html
<div style="width:100%; margin: auto;">

<img src="https://cloud.s-ccs.de/public.php/dav/files/nDQXYteFgFXrAjj/" style="width:45%;">
<img src="https://cloud.s-ccs.de/public.php/dav/files/gAAaaRdSebCY4fd"  style="width:45%;">

</div>
```

## Key features

- **Overlap correction:** Multiple ways to model overlap between temporally close events
- **📈 Regression ERPs:** Fit linear and non-linear predictors, mass univariate models, define contrasts, calculate marginal effects
- **🧠 Intuitive:** Easy to specify models (`w~i+lcox` formulas), easy to get 🧹 tidy results
- **⚡ Fast & modular:** Many solvers, GPU support, easily extensible
- **🌍 Ecosystem:** A diverse ecosystem allows for mixed-models, decoding, statistics, plotting and simulation

## Installation

```julia-repl
julia> using Pkg; Pkg.add("Unfold")
```

For more detailed instructions please refer to [Installing Julia & Unfold.jl](@ref install_instruct).

## Quick start

There are four main model types

1. Timeexpansion **No**, Mixed **No**  : `fit(UnfoldModel, [Any=>(f, -0.1:0.01:0.5)], evts, data_epoch)`
1. Timeexpansion **Yes**, Mixed **No** : `fit(UnfoldModel, [Any=>(f, basisfunction)], evts, data)`
1. Timeexpansion **No**, Mixed **Yes** : `fit(UnfoldModel, [Any=>(fLMM, -0.1:0.01:0.5)], evts, data_epoch)`
1. Timeexpansion **Yes**, Mixed **Yes**: `fit(UnfoldModel, [Any=>(fLMM, basisfunction)], evts, data)`

## Usage example

### rERP model

```julia
using Unfold
using UnfoldSim
data, evts = UnfoldSim.predef_eeg()

f = @formula 0 ~ 1 + condition
basisfunction = firbasis(τ = (-0.1,0.5), sfreq = 100)
fit(UnfoldModel, [Any=>(f, basisfunction)], evts, data)
nothing #hide
```

### MixedModels

It is also possible to fit Linear Mixed Models using the sister-package [UnfoldMixedModels.jl](https://unfoldtoolbox.github.io/UnfoldDocs/UnfoldMixedModels.jl/stable/)

```julia
using UnfoldMixedModels
using UnfoldSim
data, evts = UnfoldSim.predef_eeg(10;return_epoched=true) # 10 subjects
data = reshape(data,size(data,1),:) # concatenate subjects

times = range(-0.1,0.5,size(data,1)) # arbitrary time-vector

fLMM = @formula 0 ~ 1 + condition + (1|subject) + (1|item)
fit(UnfoldModel, [Any=>(f, times)], evts, data)
nothing #hide
```

## Where to start: Learning roadmap

### 1. First step

📌 Goal: First run a mass-univariate analysis, similar to ERPs. Then add the overlap correction. Also very common is to have multiple-events (Stimulus, Response, Fixations etc.)\
🔗 [first the mass-univariate approach](@ref lm_massunivariate) - then [theoverlap-correction tutorial](@ref lm_overlap) -  finally [multiple events](@ref multievent)

### 2. Intermediate topics

📌 Goal: Next familiarize yoursel with marginal effects, and potentially non-linear spline modelling. Defining contrasts can also be helpful \
🔗  [marginal effects](@ref effects) - [non linear effects](@ref nonlinear) - [contrast coding](@ref contrasts)

### 3. Advanced topics

📌 Goal: There are a lot of advanced topics in Unfold.jl, learn how to use the GPU or outlier-robust solvers, or define your own solver \
🔗 [GPU and robust models](@ref custom_solvers) - [solver definition](@ref solver_implementation)

## Statement of need

```@raw html
<!---
Note: The statement of need is also used in the `README.md`. Make sure that they are synchronized.
-->
```
