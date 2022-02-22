# Custom Solvers

### Setup some data

```@Example main
using Unfold
using UnfoldMakie, CairoMakie
include(joinpath(dirname(pathof(Unfold)), "../test/test_utilities.jl") ) # to load data
dat, evts = loadtestdata("test_case_3b");

basisfunction = firbasis(τ=(-0.4,.8),sfreq=50,name="stimulus")
f  = @formula 0~1+conditionA+conditionB
bfDict = Dict(Any=>(f,basisfunction))

```

### Custom Solver with standard error
```@Example main
se_solver =(x,y)->Unfold.solver_default(x,y,stderror=true)
m = Unfold.fit(UnfoldLinearModel,bfDict,evts,dat,solver=se_solver)
results =coeftable(m)
plot_results(results)
```
!!! warning
    Use single-subject SE on your own risk. Because EEG data are autocrrelated YOUR SE WILL BE TOO SMALL!

### Back2Back regression
```@Example main
b2b_solver = (x, y) -> Unfold.solver_b2b(x, y,cross_val_reps = 5)
m = Unfold.fit(UnfoldLinearModel, f, events, dat, times, solver=b2b_solver)
results = coeftable(m)

plot_results(results)
```
These are the decoding-results for `conditionA` while taking into account `conditionB` - and vice versa. Not very exciting right now because in this simple example we only have one channel ;-)


