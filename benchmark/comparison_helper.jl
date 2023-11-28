# comparison_helper functions

using MATLAB
using UnfoldSim
using Random


# calculate simulations

function runsim(design, ovlap::Tuple, multi_channel)
    if ~multi_channel
        p1 = LinearModelComponent(;
            basis=p100(),
            formula=@formula(0 ~ 1),
            β=[5]
        )

        n1 = LinearModelComponent(;
            basis=n170(),
            formula=@formula(0 ~ 1 + condition),
            β=[5, -3]
        )

        p3 = LinearModelComponent(;
            basis=p300(),
            formula=@formula(0 ~ 1 + continuous),
            β=[5, 1]
        )

        components = [p1, n1, p3]
        data, events = simulate(MersenneTwister(1), design, components, UniformOnset(; width=ovlap[1], offset=ovlap[2]), PinkNoise())

    else
        c = LinearModelComponent(; basis=p100(), formula=@formula(0 ~ 1 + condition), β=[5, 1])
        c2 = LinearModelComponent(; basis=p300(), formula=@formula(0 ~ 1 + continuous), β=[5, -3])

        hart = headmodel(type="hartmut")
        mc = UnfoldSim.MultichannelComponent(c, hart => "Left Postcentral Gyrus")
        mc2 = UnfoldSim.MultichannelComponent(c2, hart => "Right Occipital Pole")

        data, events = simulate(MersenneTwister(1), design, [mc, mc2], UniformOnset(; width=ovlap[1], offset=ovlap[2]), PinkNoise())
    end
    return data, events
end

# run Matlab code
function calc_matlab(datajl, eventsjl)
    mat"
        EEG = eeg_emptyset();
        
        EEG.data = $datajl;
        EEG.srate = 100;
    "

    mat_events = mxarray(Dict("continuous" => eventsjl.continuous, "condition" => eventsjl.condition, "latency" => eventsjl.latency))

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