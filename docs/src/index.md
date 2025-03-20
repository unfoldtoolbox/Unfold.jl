```@meta
CurrentModule = Unfold
```

# Unfold.jl Documentation

Welcome to [Unfold.jl](https://github.com/unfoldtoolbox/Unfold.jl): 

If you want to follow the **tutorials**, best to start with the [mass-univariate approach](@ref lm_massunivariate), which should be familiar to you if you did ERPs before. Then the [overlap-correction tutorial](@ref lm_overlap). If you are then not satisfied, check out more advanced topics: [effects-interface (aka what to do after fitting)](@ref effects), or non-linear effects.

In case you want to understand the tools better, check out our **explanations**.

Once you are familiar with the tools, check out further **how-to guides** for specific applications.

In case you want to understand the toolbox better, we plan to offer **technical references**. This includes Benchmarks & Explorations.


```@raw html
<div style="width:60%; margin: auto;">

</div>
```

## Key features
- **Modularity:**
- **:**


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
📌 Goal: 
🔗 

### 2. Intermediate topics
📌 Goal: 
🔗

### 3. Advanced topics
📌 Goal: 
🔗


## Statement of need
```@raw html
<!---
Note: The statement of need is also used in the `README.md`. Make sure that they are synchronized.
-->
```