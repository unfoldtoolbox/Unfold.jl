# Custom Solvers

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

### Back2Back regression
```@Example main
b2b_solver = (x, y) -> Unfold.solver_b2b(x, y;ross_val_reps = 5)
dat_3d = permutedims(repeat(dat,1,1,20),[3 1 2])
m = Unfold.fit(UnfoldModel, designDict, evts,dat_3d; solver=b2b_solver)
results = coeftable(m)

plot_erp(results)
```
These are the decoding-results for `conditionA` while taking into account `conditionB` - and vice versa. 


