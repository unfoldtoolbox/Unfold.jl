# Calling Matlab from Julia and comparing to Unfold.jl

ENV["MATLAB_ROOT"] = "/opt/common/apps/matlab/r2024a/"

using MATLAB
using BSplineKit, Unfold, UnfoldSim
using Random
using Chairmarks
using DataFrames

mat"addpath('/store/users/skukies/unfold_duration_backup/lib/eeglab')" # Add EEGLAB?
mat"addpath('/store/users/skukies/unfold_duration_backup/lib/unfold')" #" Add Unfold?
mat"init_unfold"

# Matlab function for comparison

#---


function calc_matlab(datajl, eventsjl)
    mat"
        
        EEG = eeg_emptyset();
        
        EEG.data = $datajl;
        EEG.srate = double($srate);
    "

    mat_events = mxarray(
        Dict(
            "continuous" => eventsjl.continuous,
            "condition" => eventsjl.condition,
            "latency" => eventsjl.latency,
        ),
    )

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
        t_design = toc
        EEG = uf_glmfit(EEG);
        t_fit = toc
    "
end

#---
# New benchmark comparison to cuda/solver_comparison.jl

include("../generate_data.jl")
srate = 100
X, y, events = benchmark_data(
    n_channels = 128,
    sfreq = srate,
    n_splines = (4, 4),
    n_repeats = 100,
    return_events_too = true,
)


calc_matlab(y, events)
@mget(t_fit)

@mget t_design
#----

# julia
@time(
    m = fit(
        UnfoldModel,
        [
            Any => (
                @formula(0 ~ 1 + condition + spl(continuous, 5)),
                firbasis(τ = [-1, 1], sfreq = srate),
            ),
        ],
        events,
        y,
    )
)

#2024-11-15: 124s in matlab, 120s in julia; single threaded
#---
using CUDA
@time(
    m = fit(
        UnfoldModel,
        [
            Any => (
                @formula(0 ~ 1 + condition + spl(continuous, 5)),
                firbasis(τ = [-1, 1], sfreq = srate),
            ),
        ],
        events,
        CuArray(y);
        solver = (x, y) -> Unfold.solver_predefined(x, y; solver = :qr),
    )
)

#2024-11-15: 1.2-1.4s