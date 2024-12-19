# How to model multiple events

When dealing with overlapping data, it is often necessary to model multiple eventtypes (e.g. fixations, stimuli, responses).

## Load Example Data

```@example main
using Unfold
using UnfoldMakie, CairoMakie
using DataFrames
using StatsModels
using MixedModels
using DisplayAs # hide

include(joinpath(dirname(pathof(Unfold)), "../test/test_utilities.jl")) # to load data
dat, evts = loadtestdata("test_case_4b");

evts[1:5,:]
```

The `type` column of table `evts` contains two conditions: `eventA`` and`eventB` (if your eventstypes are specified in a different column, you need to define the keywordargument `eventcolumn` in the `fit` command below)

## Specify formulas and basisfunctions

```@example main

bf1 = firbasis(τ = (-0.4, 0.8), sfreq = 50)
bf2 = firbasis(τ = (-0.2, 1.2), sfreq = 50)
bf2|> DisplayAs.withcontext(:is_pluto=>true) # hide
```

For each event, a basis function and formula must be specified. The same basis and formulas may be used.

```@example main
f  = @formula 0 ~ 1
```

For each event, we must specify the formula and basis function to be used.

```@example main

bfDict = [ "eventA" => (f, bf1),
           "eventB" => (f, bf2) ]

bfDict |> DisplayAs.withcontext(:is_pluto=>true) # hide
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
plot_erp(results; stderror = true, mapping = (; col = :eventname))
```
