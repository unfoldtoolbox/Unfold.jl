# # Window length effects

using Unfold, UnfoldSim
using CairoMakie, AlgebraOfGraphics, MakieThemes
using Random
using DataFrames, DataFramesMeta
using ColorSchemes, Colors

# !!! important
#     For analyzing real-world EEG data we recommend that researchers should — a priori — make an educated guess about the length of the underlying EEG activity and select this as their EW. This also suggests to use event windows with different sizes between events (as is possible with Unfold).
#     Further, as can be seen below, when choosing longer time-windows the overfit is only of moderate size, thus we additionally recommend to generally err on the longer side, to not miss any important activity. \
#     For a more in depth explanation on this, you can read our 2023 CCN paper: [Skukies & Ehinger, 2023](https://www.biorxiv.org/content/10.1101/2023.06.05.543689v1)

set_theme!(theme_ggthemr(:fresh))

# As opposed to classical averaged ERPs overlap corrected regression ERPs can be influenced by the chosen window length:
# Long estimation windows might capture all relevant event-related activity, but might introduce artifacts due to overfit, 
# short estimation windows might not overfit, but also might not capture all (overlapping) activity, and thereby introduce bias.
#
# Thus a common question we get is, how to specify the length of the estimation windows. 


# # Init functions

# First we need a function that simulates some continous data; conviently we can use UnfoldSim for this 
function gen_data(rng, noiselevel, sfreq)
    noise = PinkNoise(; noiselevel = noiselevel)

    dat, evts = UnfoldSim.predef_eeg(
        rng;
        sfreq = sfreq,
        p1 = (p100(; sfreq = sfreq), @formula(0 ~ 1 + condition), [5, 0], Dict()),
        n1 = (n170(; sfreq = sfreq), @formula(0 ~ 1 + condition), [5, 0], Dict()),
        p3 = (p300(; sfreq = sfreq), @formula(0 ~ 1 + continuous), [5, 0], Dict()),
        n_repeats = 20,
        noise = noise,
    )
    return dat, evts
end;

# Next a convience function to calculate the estimates
function calc_time_models(evts, dat, tWinList, sfreq)
    mList = []
    for twindow in tWinList
        m = fit(
            UnfoldModel,
            Dict(Any => (@formula(0 ~ 1), firbasis(twindow, sfreq))),
            evts,
            dat,
        )
        res = coeftable(m)
        res.tWin .= string.(Ref(twindow[2]))
        push!(mList, res)
    end
    return vcat(mList...)
end;


# # Init variables
tWinList = [(-0.1, x) for x in [3, 2.5, 2, 1.5, 1, 0.5]]
noiselevel = 8.5
sfreq = 250;

# # Generate data and calculate estimates

dat, evts = gen_data(MersenneTwister(2), noiselevel, sfreq);

res = calc_time_models(evts, dat, tWinList, sfreq);

# We also append some additional information to the results dataframe

# For comparison lets also generate the ground truth of our data; this is a bit cumbersome and you don't have to care (too much) about it
dat_gt, evts_gt = UnfoldSim.predef_eeg(;
    p1 = (p100(; sfreq = sfreq), @formula(0 ~ 1), [5], Dict()),
    sfreq = sfreq,
    n1 = (n170(; sfreq = sfreq), @formula(0 ~ 1), [5], Dict()),
    p3 = (p300(; sfreq = sfreq), @formula(0 ~ 1), [5], Dict()),
    n_repeats = 1,
    noiselevel = 0,
    return_epoched = true,
);
time_gt = range(0, length = length(dat_gt[:, 1]), step = 1 / sfreq)
unique_event = unique(res.tWin)
df_gt = DataFrame(
    tWin = reduce(vcat, fill.(unique_event, length(dat_gt[:, 1]))),
    eventname = Any,
    channel = repeat([1], length(dat_gt[:, 1]) * length(unique_event)),
    coefname = reduce(
        vcat,
        fill("GroundTruth", length(dat_gt[:, 1]) * length(unique_event)),
    ),
    estimate = repeat(dat_gt[:, 1], length(unique_event)),
    group = reduce(vcat, fill(nothing, length(dat_gt[:, 1]) * length(unique_event))),
    stderror = reduce(vcat, fill(nothing, length(dat_gt[:, 1]) * length(unique_event))),
    time = repeat(time_gt, length(unique_event)),
);



# And append ground truth to our results df	
res_gt = vcat(res, df_gt);


# # Plot results

# Choose which data to plot
h_t =
    AlgebraOfGraphics.data(res) * mapping(
        :time,
        :estimate,
        color = :tWin,
        group = (:tWin, :coefname) => (x, y) -> string(x[2]) * y,
    );

# We use the following to plot some length indicator lines

untWin = unique(res_gt.tWin)
segDF = DataFrame(
    :x => hcat(repeat([-0.1], length(untWin)), parse.(Float64, untWin))[:],
    :y => repeat(reverse(1:length(untWin)), outer = 2),
)
segDF.tWin .= "0.0"
segDF.tWin .= segDF.x[reverse(segDF.y .+ 6)]
segDF.y = segDF.y .* 0.2 .+ 6;


# Layer for indicator lines
h_l =
    AlgebraOfGraphics.data(@subset(segDF, :tWin .!= "3.0")) *
    mapping(:x, :y, color = :tWin, group = :tWin => x -> string.(x));

# Ground truth Layer
h_gt =
    AlgebraOfGraphics.data(df_gt) *
    mapping(:time, :estimate, group = (:tWin, :coefname) => (x, y) -> string(x) * y) *
    visual(Lines, linewidth = 5, color = Colors.Gray(0.6));

# Add all visuals together and draw
h1 =
    h_gt + visual(Lines, colormap = get(ColorSchemes.Blues, 0.3:0.01:1.2)) * (h_l + h_t) |>
    x -> draw(x, axis = (; xlabel = "time [s]", ylabel = "estimate [a.u.]"));

# Add zero grid lines
h1 = hlines!(current_axis(), [0], color = Colors.Gray(0.8));
h2 = vlines!(current_axis(), [0], color = Colors.Gray(0.8));
translate!(h1, 0, 0, -1);
translate!(h2, 0, 0, -1);

# Plot figure
current_figure()
