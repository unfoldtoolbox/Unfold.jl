# [Alternative Solvers](@id custom_solvers)

A solver takes an Unfold-specified DesignMatrix and the data, and typically solves the equation system `y = Xb` (in the case of Linear Models). There are many different ways how one can approach this problem, depending if the matrix is sparse, if it is 2D or 3D, if one wants to use GPU etc.

Most implemented solvers ultimately make use of `solver_main` for their main loop. See the `reference` tutorial for more information if that is interesting to you.

## Setup some data

```@Example main
using Unfold
using UnfoldMakie, CairoMakie
using UnfoldSim
dat, evts = UnfoldSim.predef_eeg(; noiselevel = 10, return_epoched = true)

f = @formula 0 ~ 1 + condition + continuous
designDict = Dict(Any => (f, range(0, 1, length = size(dat, 1))))
```

## GPU Solvers

GPU solvers can significantly speed up your model fitting, with observed improvements of up to a factor of 30-100!

## fastest GPU solver

Empirically we found that solving `X'Xb = X'y` is the fastest way to solve for `b`. To achieve this, you can run:

```julia
using CUDA
gpu_solver =(x, y) -> Unfold.solver_predefined(x, y; solver=:qr)
m = Unfold.fit(UnfoldModel, designDict, evts, cu(dat), solver = gpu_solver)
```

Where the `cu` is the magic that moves the data to the GPU. Internatlly, the solver function will move the matrix as well and pre-calculate some matrices (especially `X'X`, `X'` and allocate `X'y`).

## lsmr GPU solver

the `Krylov.lsmr` implementation directly solves `y = Xb`, but allows for running on the GPU.

```julia
using Krylov, CUDA # necessary to load the right package extension
gpu_solver =(x, y) -> Unfold.solver_krylov(x, y; GPU = true)
m = Unfold.fit(UnfoldModel, designDict, evts, dat, solver = gpu_solver)
```

To test it, you will need to run it yourself as we cannot run it on the docs. If you require a different graphicscard vendor than NVIDA/CUDA, please create an issue. Currently, we are unable to test it due to lack of hardware.

## Robust Solvers

Robust solvers automatically adjust for outlier trials, but they come at a significant computational cost.

```@Example main
using RobustModels # necessary to load the Unfold package extension
se_solver = (x, y) -> Unfold.solver_robust(x, y)
m = Unfold.fit(UnfoldModel, designDict, evts, dat, solver = se_solver)
results = coeftable(m)
plot_erp(results; stderror = true)
```

## Back2Back regression

Since 2025, this solver requires `UnfoldDecode` - please find the tutorial and explanation there, the example here is for historic reasons and will be removed at a later point.

```@Example
using UnfoldDecode
b2b_solver = (x, y) -> UnfoldDecode.solver_b2b(x, y; ross_val_reps = 5)
dat_3d = permutedims(repeat(dat, 1, 1, 20), [3 1 2])
m = Unfold.fit(UnfoldModel, designDict, evts, dat_3d; solver = b2b_solver)
results = coeftable(m)

plot_erp(results)
```

These are the decoding results for `conditionA` while considering `conditionB`, and vice versa.
