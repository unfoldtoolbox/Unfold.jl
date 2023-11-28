# Calling Matlab from Julia and comparing to Unfold.jl
using Pkg
Pkg.activate("./benchmark")

ENV["MATLAB_ROOT"] = "/opt/common/apps/matlab/r2021a/"
ENV["JULIA_NUM_THREADS"] = 20

includet("./comparison_helper.jl")
using MATLAB
using BSplineKit, Unfold, UnfoldSim
using Random

mat"addpath('/store/users/skukies/unfold_duration/lib/eeglab')" # Add EEGLAB?
mat"addpath('/store/users/skukies/unfold_duration/lib/unfold')" #" Add Unfold?
mat"init_unfold"

# Parameters
rep_design = 50;
sfreq = 100;
ovlap = (50,20);
jl_formula = @formula(0~1+condition+spl(continuous,5));
#mat_formula = mxarray("y~1+cat(condition+spl(continuous,5)"); # It seems the formula has to be hardcoded in the function :/

#" Setting up data
# data, events = UnfoldSim.predef_eeg()
design = SingleSubjectDesign(;
        conditions=Dict(:condition=>["car","face"],:continuous=>range(-5,5,length=10))
        ) |> x->RepeatDesign(x,rep_design);

data, events = runsim(design, ovlap, false);

# Run fit two times; once for compilation and once for actual timing
m = fit(UnfoldModel,Dict(Any=>(jl_formula,firbasis(τ=[-1,1],sfreq=100,name="basis"))),events,data);
@time m = fit(UnfoldModel,Dict(Any=>(jl_formula,firbasis(τ=[-1,1],sfreq=100,name="basis"))),events,data);

## Matlab part

calc_matlab(data, events)

# Above tested XXXX/XX/XX; Julia: XX sec ; Matlab: XX sec

## Multichannel w/ headmodel
data_mc, events_mc = runsim(design, ovlap, true)

# Fit Unfold.jl
@time m = fit(UnfoldModel,Dict(Any=>(jl_formula,firbasis(τ=[-1,1],sfreq=100,name="basis"))),events_mc,data_mc);

# Fit Matlab
calc_matlab(data_mc, events_mc)

# Above tested XXXX/XX/XX; Julia: xx sec ; Matlab: ca. xx sec
