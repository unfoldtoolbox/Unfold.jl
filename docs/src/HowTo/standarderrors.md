# [Standard Errors](@id standard_errors)

### Setup some data

```@Example main
using Unfold
using UnfoldMakie, CairoMakie
using UnfoldSim,
dat,evts = UnfoldSim.predef_eeg(;noiselevel=10,return_epoched=true)

f  = @formula 0~1+condition+continuous
designDict = Dict(Any=>(f,range(0,1,length=size(dat,1))))

```

Similar to @Ref(custom_solver), it is possible to specify a solver than also calculates the standard errors of the (single-subject) estimates.

```@Example main
se_solver =(x,y)->Unfold.solver_default(x,y,stderror=true)
m = Unfold.fit(UnfoldModel,designDict,evts,dat,solver=se_solver)
results =coeftable(m)
plot_erp(results;stderror=true)
```
!!! warning
    **In case of overlap-correction:** Use single-subject SE on your own risk. Because EEG data are autocorrelated your standarderrors will be typically too small.

