# Calling Matlab from Julia and comparing to Unfold.jl

ENV["MATLAB_ROOT"] = "/opt/common/apps/matlab/r2021a/"

using MATLAB
using Unfold, UnfoldSim
using Random

mat"addpath('/store/users/skukies/unfold_duration/lib/eeglab')" # Add EEGLAB?
mat"addpath('/store/users/skukies/unfold_duration/lib/unfold')" # Add Unfold?
mat"init_unfold"

## Julia
data, events = UnfoldSim.predef_eeg()

# Run fit two times; once for compilation and once for actual timing
m = fit(UnfoldModel,Dict(Any=>(@formula(0~1+condition+continuous),firbasis(τ=[-0.1,1],sfreq=100,name="basis"))),events,data);
@time m = fit(UnfoldModel,Dict(Any=>(@formula(0~1+condition+continuous),firbasis(τ=[-0.1,1],sfreq=100,name="basis"))),events,data);

## Matlab part

mat"
    EEG = eeg_emptyset();
    EEG.data = repmat($data',30,1);
    EEG.data = EEG.data + randn(size(EEG.data));
    EEG.srate = 100;
"

mat"
    for e = 1:length($events.latency)
        EEG.event(end+1).latency =  $events.latency(e);
        EEG.event(end).condition =  $events.condition{e};
        EEG.event(end).continuous=  $events.continuous(e);
        EEG.event(end).continuousB=  rand(1);
        EEG.event(end).type = 'fixation';
    end
"

mat"EEG = eeg_checkset(EEG)"

mat"    
    cfgDesign = [];
    cfgDesign.eventtypes = {'fixation'};
    cfgDesign.formula = 'y ~ 1+ cat(condition)+spl(continuous,5)+spl(continuousB,20)';
    tic
    EEG = uf_designmat(EEG,cfgDesign);

    EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-1,1]);
    toc
    EEG = uf_glmfit(EEG);
    toc
"