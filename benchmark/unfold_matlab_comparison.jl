# Calling Matlab from Julia and comparing to Unfold.jl

ENV["MATLAB_ROOT"] = "/opt/common/apps/matlab/r2021a/"

using MATLAB
using BSplineKit, Unfold, UnfoldSim
using Random

mat"addpath('/store/users/skukies/unfold_duration/lib/eeglab')" # Add EEGLAB?
mat"addpath('/store/users/skukies/unfold_duration/lib/unfold')" #" Add Unfold?
mat"init_unfold"

# Matlab function for comparison
function calc_matlab(datajl, eventsjl)
    mat"
        EEG = eeg_emptyset();
        
        EEG.data = $datajl;
        EEG.srate = 100;
    "

    mat_events = mxarray(Dict("continuous"=>eventsjl.continuous, "condition"=>eventsjl.condition, "latency"=>eventsjl.latency))

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
end


#" Setting up data
# data, events = UnfoldSim.predef_eeg()
design = SingleSubjectDesign(;
        conditions=Dict(:condition=>["car","face"],:continuous=>range(-5,5,length=10))
        ) |> x->RepeatDesign(x,50);

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
data,events = simulate(MersenneTwister(1),design,components,UniformOnset(;width=50,offset=20),PinkNoise());

# Run fit two times; once for compilation and once for actual timing
m = fit(UnfoldModel,Dict(Any=>(@formula(0~1+condition+continuous),firbasis(τ=[-1,1],sfreq=100,name="basis"))),events,data);
@time m = fit(UnfoldModel,Dict(Any=>(@formula(0~1+condition+spl(continuous,5)),firbasis(τ=[-1,1],sfreq=100,name="basis"))),events,data);

## Matlab part

calc_matlab(data, events)

# Above tested 2023/11/22; Julia: 0.12 sec ; Matlab: 0.15sec

## Multichannel w/ headmodel

c = LinearModelComponent(;basis=p100(),formula = @formula(0~1+condition),β = [5,1]);
c2 = LinearModelComponent(;basis=p300(),formula = @formula(0~1+continuous),β = [5,-3]);

hart = headmodel(type="hartmut")
mc = UnfoldSim.MultichannelComponent(c, hart=>"Left Postcentral Gyrus")
mc2 = UnfoldSim.MultichannelComponent(c2, hart=>"Right Occipital Pole")

onset = UniformOnset(;width=50,offset=20);

data_mc,events_mc = simulate(MersenneTwister(1),design, [mc,mc2],  onset, PinkNoise())

# Fit Unfold.jl
@time m = fit(UnfoldModel,Dict(Any=>(@formula(0~1+condition+spl(continuous,5)),firbasis(τ=[-1,1],sfreq=100,name="basis"))),events_mc,data_mc);

# Fit Matlab
calc_matlab(data_mc, events_mc)

# Above tested 2023/11/22; Julia: ca.8.9 sec ; Matlab: ca. 10.5sec
