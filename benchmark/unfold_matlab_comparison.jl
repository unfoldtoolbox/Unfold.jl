# Calling Matlab from Julia and comparing to Unfold.jl
# To use this, make sure to start Julia with multiple threads e.g. 'julia --threads 20'
#using Pkg
#Pkg.activate(".") # Make sure to use the right project (benchmark) and instantiate

ENV["MATLAB_ROOT"] = "/opt/common/apps/matlab/r2021a/" # Points to Matlab folder containing executable
#ENV["JULIA_NUM_THREADS"] = 20

include("./comparison_helper.jl")
using MATLAB
using BSplineKit, Unfold, UnfoldSim
using Random
using Krylov,CUDA

mat"addpath('/store/users/skukies/unfold_duration/lib/eeglab')" # Add path to EEGLAB here
mat"addpath('/store/users/skukies/unfold_duration/lib/unfold')" #" Add path to Unfold Matlab
mat"init_unfold"

# Parameters
rep_design = 50;
sfreq = 100;
ovlap = (50,20);
jl_formula = @formula(0~1+condition+spl(continuous,5));
design_dict = Dict(Any=>(jl_formula,
                firbasis(τ=[-1,1],sfreq=100,name="basis")));

#mat_formula = mxarray("y~1+cat(condition+spl(continuous,5)"); # It seems the formula has to be hardcoded in the function :/

#" Setting up data
# data, events = UnfoldSim.predef_eeg()
design = SingleSubjectDesign(;
        conditions=Dict(:condition=>["car","face"],:continuous=>range(-5,5,length=10))
        ) |> x->RepeatDesign(x,rep_design);

data, events = runsim(design, ovlap, false);

# Run fit two times; once for compilation and once for actual timing
m = fit(UnfoldModel,design_dict,events,data);
@time m = fit(UnfoldModel,design_dict,events,data);

## Matlab part

calc_matlab(data, events)

# Above tested XXXX/XX/XX; Julia: XX sec ; Matlab: XX sec

## Multichannel w/ headmodel
data_mc, events_mc = runsim(design, ovlap, true)

# Fit Unfold using GPU (again, run twice to deal with compilation time)
gpu_solver =(x,y)->Unfold.solver_krylov(x,y;GPU=true)
m_gpu = Unfold.fit(UnfoldModel,design_dict,events_mc,data_mc,solver=gpu_solver);
@time m_gpu = Unfold.fit(UnfoldModel,design_dict,events_mc,data_mc,solver=gpu_solver);

# Fit Unfold.jl
@time m = fit(UnfoldModel,design_dict,events_mc,data_mc);

# Fit Matlab
calc_matlab(data_mc, events_mc)


# Above tested 2023/11/29; Julia: 15 sec ; Julia GPU: 11 sec ; Matlab: ca. 68 sec
