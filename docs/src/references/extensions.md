# Package-extensions

In  Julia 1.9 Package Extensions were introduced. Unfold.jl is making use of them in four ways.
Prior to using some functionality, you have to add + load specific package(s) for the functionality to be available. The reason for this is, that if you don't need e.g. GPU-support, you also will not need to install it.

## GPU: Krylov,CUDA

To use gpu support as described in @Ref(custom_solvers) you have to:

```julia
using Krylov,CUDA
using Unfold
```

## RobustSolvers.jl

To use robust (outlier-"safe") solvers support as described in @Ref(custom_solvers) you have to:

```julia
using RobustSolvers
using Unfold
```

## Non-linear effects: BSplineKit.jl

Finally to use non-linear effects/splines like in `@formula 0~1+spl(continuous,5)` you have to use:

```julia
using BSplineKit
using Unfold
```

!!! note
    In principle you should be able to load the package after loading Unfold. But sometimes this doesnt work, a `Base.retry_load_extensions()` call might help in these situations.
