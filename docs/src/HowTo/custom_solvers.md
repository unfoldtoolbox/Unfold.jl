# [Custom Solvers](@id custom_solvers)

### Setup some data

```@Example main
using Unfold
using UnfoldMakie, CairoMakie
using UnfoldSim,
dat,evts = UnfoldSim.predef_eeg(;noiselevel=10,return_epoched=true)

f  = @formula 0~1+condition+continuous
designDict = Dict(Any=>(f,range(0,1,length=size(dat,1))))

```

### Custom Solver with standard error
```@Example main
se_solver =(x,y)->Unfold.solver_default(x,y,stderror=true)
m = Unfold.fit(UnfoldModel,designDict,evts,dat,solver=se_solver)
results =coeftable(m)
plot_erp(results;extra=(;stderror=true))
```
!!! warning
    Use single-subject SE on your own risk. Because EEG data are autocrrelated YOUR SE WILL BE TOO SMALL!





### GPU Solvers
GPU solvers can speed up your modelfit drastically! up to factor of 30 has been observed already
```julia
using Krylov,CUDA # necessary to load the right package extension
gpu_solver =(x,y)->Unfold.solver_krylov(x,y;GPU=true)
m = Unfold.fit(UnfoldModel,designDict,evts,dat,solver=gpu_solver)
```
We can't run it on the docs though, so try it yourself! If you need something else than CUDA, write an issue, we cant test it with something else right now...


### Robust Solvers
Robust solvers automatically account for outlying trials. They come at a severe computational cost though!
```@Example main
using RobustModels # necessary to load the Unfold st
package extension
se_solver =(x,y)->Unfold.solver_robust(x,y)
m = Unfold.fit(UnfoldModel,designDict,evts,dat,solver=se_solver)
results =coeftable(m)
plot_erp(results;extra=(;stderror=true))
```

### Back2Back regression
```@Example main
b2b_solver = (x, y) -> Unfold.solver_b2b(x, y;ross_val_reps = 5)
dat_3d = permutedims(repeat(dat,1,1,20),[3 1 2])
m = Unfold.fit(UnfoldModel, designDict, evts,dat_3d; solver=b2b_solver)
results = coeftable(m)

plot_erp(results)
```
These are the decoding-results for `conditionA` while taking into account `conditionB` - and vice versa. 


