# How to model multiple events

When dealing with overlapping data, it is often necessary to plot multiple events.

### Load Example Data
```@example main
using Unfold
using UnfoldMakie, CairoMakie
using DataFrames
using StatsModels
using MixedModels

include(joinpath(dirname(pathof(Unfold)), "../test/test_utilities.jl")) # to load data
dat, evts = loadtestdata("test_case_4b");

evts[1:5,:]
```
The `type` column of table `evts` contains two conditions: EventA and EventB. By default, the toolbox will search for these conditions.

### Specify formulas and basisfunctions

```@example main
bf1 = firbasis(τ = (-0.4, 0.8), sfreq = 50, name = "stimulusA")
bf2 = firbasis(τ = (-0.2, 1.2), sfreq = 50, name = "stimulusB")
```
For each event, a basis function and formula must be specified. The same basis and formulas may be used.
```@example main
f  = @formula 0 ~ 1
```

For each event, we must specify the formula and basis function to be used. 
```@example main
bfDict = Dict("eventA" => (f, bf1), "eventB" => (f, bf2))
```

Finally, fitting & plotting works the same way as always

```@example main
m = Unfold.fit(
    UnfoldModel,
    bfDict,
    evts,
    dat,
    solver = (x, y) -> Unfold.solver_default(x, y; stderror = true),
    eventcolumn = "type",
)
results = coeftable(m)
plot_erp(results; stderror = true, mapping = (; col = :group))
``` 
