# How to model multiple events

In case of overlapping data, we often need to plot multiple events.


### Load Example Data
```@example
using Unfold
using UnfoldMakie
using StatsModels
using MixedModels

include(joinpath(dirname(pathof(MyModule)), "test/test_utilities") ) # to load data
dat, evts = loadtestdata("test_case_4b");

head(evts)
```
As you can see, there are two events here. EventA and EventB. Both are in the column `type` in which the toolbox looks for different events by default.

### Specify formulas and basisfunctions

```@example main
bf1 = firbasis(τ=(-0.4,.8),sfreq=50,name="stimulusA")
bf2 = firbasis(τ=(-0.2,1.2),sfreq=50,name="stimulusB")
```
For each event, we have to specify a basisfunction and a formula. We could use the same basis and the same formulas though
```@example main
f  = @formula 0~1
```

Next, we have to specify for each event, what is the formula and what is the basisfunction we want to use
```@example main
bfDict = Dict(:EventA=>(f,bf1),
              :EventB=>(f,bf2))
```

Finally, fitting & plotting works the same way as always
```@example main

m,results = Unfold.fit(UnfoldLinearModel,bfDict,evts,data,solver=se_solver)
plot_results(results,stderror=true)
``` 
