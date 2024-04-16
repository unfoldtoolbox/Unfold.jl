# [Alternative Solvers](@id custom_solvers)

### Setup some data

```@Example main
using Unfold
using UnfoldMakie, CairoMakie
using UnfoldSim
dat, evts = UnfoldSim.predef_eeg(; noiselevel = 10, return_epoched = true)

f = @formula 0 ~ 1 + condition + continuous
designDict = Dict(Any => (f, range(0, 1, length = size(dat, 1))))
```

### GPU Solvers
GPU solvers can significantly speed up your model fitting, with observed improvements of up to a factor of 30!

```julia
using Krylov, CUDA # necessary to load the right package extension
gpu_solver =(x, y) -> Unfold.solver_krylov(x, y; GPU = true)
m = Unfold.fit(UnfoldModel, designDict, evts, dat, solver = gpu_solver)
```
To test it, you will need to run it yourself as we cannot run it on the docs. If you require a different graphicscard vendor than NVIDA/CUDA, please create an issue. Currently, we are unable to test it due to lack of hardware.

### Robust Solvers
Robust solvers automatically account for outlier trials, but they come at a significant computational cost.
```@Example main
using RobustModels # necessary to load the Unfold package extension
se_solver = (x, y) -> Unfold.solver_robust(x, y)
m = Unfold.fit(UnfoldModel, designDict, evts, dat, solver = se_solver)
results = coeftable(m)
plot_erp(results; stderror = true)
```

### Back2Back regression
```@Example main
b2b_solver = (x, y) -> Unfold.solver_b2b(x, y; ross_val_reps = 5)
dat_3d = permutedims(repeat(dat, 1, 1, 20), [3 1 2])
m = Unfold.fit(UnfoldModel, designDict, evts, dat_3d; solver = b2b_solver)
results = coeftable(m)

plot_erp(results)
```
These are the decoding results for `conditionA` while considering `conditionB`, and vice versa. 


