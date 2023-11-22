# Calling Matlab from Julia and comparing to Unfold.jl

ENV["MATLAB_ROOT"] = "/opt/common/apps/matlab/r2021a/"

using MATLAB
using BSplineKit, Unfold, UnfoldSim
using Random

mat"addpath('/store/users/skukies/unfold_duration/lib/eeglab')" # Add EEGLAB?
mat"addpath('/store/users/skukies/unfold_duration/lib/unfold')" # Add Unfold?
mat"init_unfold"

## Julia
# simulate data
# data, events = UnfoldSim.predef_eeg()
design = SingleSubjectDesign(;
        conditions=Dict(:condition=>["car","face"],:continuous=>range(-5,5,length=10))
        ) |> x->RepeatDesign(x,100);

p1 =  LinearModelComponent(;
        basis = p100(),
        formula = @formula(0~1),
        β = [5]
);

n1 =  LinearModelComponent(;
        basis = n170(),
        formula = @formula(0~1+condition),
        β = [5,-3]
);

p3 =  LinearModelComponent(;
        basis = p300(),
        formula = @formula(0~1+continuous),
        β = [5,1]
);

components = [p1,n1,p3]
data,events = simulate(MersenneTwister(1),design,[p1,n1,p3],UniformOnset(;width=500,offset=200),PinkNoise());

# Run fit two times; once for compilation and once for actual timing
m = fit(UnfoldModel,Dict(Any=>(@formula(0~1+condition+continuous),firbasis(τ=[-0.1,1],sfreq=100,name="basis"))),events,data);
@time m = fit(UnfoldModel,Dict(Any=>(@formula(0~1+condition+spl(continuous,5)),firbasis(τ=[-0.1,1],sfreq=100,name="basis"))),events,data);

## Matlab part

mat"
    EEG = eeg_emptyset();
    EEG.data = repmat($data',30,1);
    EEG.data = EEG.data + randn(size(EEG.data));
    EEG.srate = 100;
"

mat_events = mxarray(Dict("continuous"=>events.continuous, "condition"=>events.condition, "latency"=>events.latency))
is_struct(mat_events)
mat"class($mat_events)"
mat"disp($mat_events)"

mat"
    for e = 1:length($mat_events.latency)
        EEG.event(end+1).latency =  $mat_events.latency(e);
        EEG.event(end).condition =  $mat_events.condition{e};
        EEG.event(end).continuous=  $mat_events.continuous(e);
        EEG.event(end).continuousB=  rand(1);
        EEG.event(end).type = 'fixation';
    end
"

mat"EEG = eeg_checkset(EEG)"

mat"    
    cfgDesign = [];
    cfgDesign.eventtypes = {'fixation'};
    cfgDesign.formula = 'y ~ 1+ cat(condition)+spl(continuous,5)';
    tic
    EEG = uf_designmat(EEG,cfgDesign);

    EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-1,1]);
    toc
    EEG = uf_glmfit(EEG);
    toc
"