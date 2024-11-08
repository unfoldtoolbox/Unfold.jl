# [Standard Errors](@id standard_errors)

### Setup some data

```@Example main
using Unfold
using UnfoldMakie, CairoMakie
using UnfoldSim
dat, evts = UnfoldSim.predef_eeg(; noiselevel = 10, return_epoched = true)

f = @formula 0 ~ 1 + condition + continuous
designDict = Dict(Any => (f, range(0, 1, length = size(dat, 1))))
```

It is possible to specify a solver that calculates the standard errors of the estimates for a single subject as it possible for [custom solvers](@ref custom_solvers).

```@Example main
se_solver = (x, y) -> Unfold.solver_default(x, y, stderror = true)
m = Unfold.fit(UnfoldModel, designDict, evts, dat, solver = se_solver)
results = coeftable(m)
plot_erp(results; stderror = true)
```
!!! warning
    **In case of overlap-correction:** Use single-subject standard errors on your own risk. EEG data is autocorrelated, which means that standard errors are typically too small.

