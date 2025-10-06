# # Circular Splines

# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages

using Unfold
using UnfoldSim
using BSplineKit
using StableRNGs
using CairoMakie, UnfoldMakie
# </details>
# ```

data, evts = UnfoldSim.predef_eeg(; return_epoched = true)
evts.cycle = rand(StableRNG(1), size(evts, 1))
times = range(0, 1, length = size(data, 1))
m = fit(UnfoldModel, @formula(0 ~ 1 + circspl(cycle, 5, 0, 1)), evts, data, times)

eff = effects(Dict(:cycle => 0:0.2:1), m)
plot_erp(eff; mapping = (; color = :cycle, group = :cycle))


data, evts = UnfoldSim.predef_eeg(; return_epoched = false)
evts.cycle = rand(StableRNG(1), size(evts, 1))
#times = range(0,1,length=size(data,1))
m = fit(
    UnfoldModel,
    [Any => (@formula(0 ~ 1 + circspl(cycle, 5, 0, 1)), firbasis((0, 1), 10))],
    evts,
    data,
)

eff = effects(Dict(:cycle => 0:0.2:1), m)
plot_erp(eff; mapping = (; color = :cycle, group = :cycle))
