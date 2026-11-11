# # [Baseline regression](@id baseline_regression)
# Classical **baseline correction** subtracts the mean of a pre-stimulus window from each epoch. This assumes the baseline is pure noise - but if the baseline co-varies with the epoch (through noise correlations, slow drifts, or the overlap of the previous trial), correcting removes that shared variance from the whole epoch and can bias the estimated effects.
#
# **Baseline regression** ([Alday et al., 2021](https://doi.org/10.7717/peerj.12197)) instead treats the baseline as a *predictor*: it is entered into the model as a channel-varying covariate, because the baseline amplitude differs per channel. In Unfold.jl, channel-varying predictors are marked with [`bychannel`](@ref): the model is fit per channel, with each channel's epochs regressed on its own baseline. No data is subtracted - the `bsl` coefficient tells you how strongly the epoch co-varies with the baseline at every latency.
#
# ```@raw html
# <details>
# <summary>Set things up</summary>
# ```

using Unfold
using UnfoldSim
using UnfoldMakie
using DataFrames
using Statistics

# ```@raw html
# </details >
# ```

# # Simulate overlapping EEG data
# We simulate multichannel EEG with `UnfoldSim.jl`, with trials overlapping so that the previous trial is not yet over when the next one starts. The simulation projects the components onto all 227 channels - for this tutorial we keep three example channels:
data, evts = UnfoldSim.predef_eeg(; multichannel = true, overlap = (0.5, 0.2))
data = data[[130, 184, 50], :]

# We drop the last events, whose epochs would run past the end of the recording, and cut the continuous data into epochs:
evts = evts[evts.latency .<= size(data, 2) - 150, :]
data_e, times = Unfold.epoch(data = data, tbl = evts, τ = (-0.2, 1.5), sfreq = 100)
size(data_e) # channel × time × epoch

# # Baseline regression
# For each channel and epoch, we compute the mean baseline and enter it into the event table as one value per channel per row:
bsl_ix = findall(t -> t <= 0.0, times)
evts.bsl = [vec(mean(data_e[:, bsl_ix, i]; dims = 2)) for i in axes(data_e, 3)]

# The formula marks `bsl` as channel-varying, and the mass-univariate model is fit (currently only supported for epoched, not for time-expanded models):
f = @formula 0 ~ 1 + condition + bychannel(bsl)
m = fit(UnfoldLinearModel, [Any => (f, times)], evts, data_e)

# The `condition` effect of the simulated N1 is recovered on all channels, with channel-specific amplitude and sign - the projection of the source onto each sensor:
ct = coeftable(m)
plot_erp(ct[startswith.(ct.coefname, "condition"), :]; mapping = (; col = :channel))

# # The beta of the baseline over time
# Now the interesting part: the `bsl` coefficient over time.
plot_erp(ct[ct.coefname.=="bsl", :]; mapping = (; col = :channel))

# How to read this plot:
# - Within the baseline window, the beta averages exactly 1 for every channel - by construction, the covariate *is* the window mean. How the 1 is distributed within the window differs per channel, which is why the baseline has to enter the model per channel.
# - Beyond the window, the beta does not return to 0: it stays non-zero across the whole epoch. The baseline co-varies with the epoch - it contains signal, not just noise. Here the signal is the previous trial's activity, which leaks into the baseline window because the trials overlap.

# # Comparison: without overlap
# We repeat the analysis with trials spaced so that the previous trial is over before the next one starts (`overlap = (2, 0.2)` - the inter-trial interval exceeds the longest component). Now the baseline is pure noise, and the beta drops to ~0 right after the window:
data, evts = UnfoldSim.predef_eeg(; multichannel = true, overlap = (2, 0.2))
data = data[[130, 184, 50], :]
evts = evts[evts.latency .<= size(data, 2) - 150, :]
data_e, times = Unfold.epoch(data = data, tbl = evts, τ = (-0.2, 1.5), sfreq = 100)
evts.bsl = [vec(mean(data_e[:, bsl_ix, i]; dims = 2)) for i in axes(data_e, 3)]
m = fit(UnfoldLinearModel, [Any => (f, times)], evts, data_e)
plot_erp(
    coeftable(m)[coeftable(m).coefname.=="bsl", :];
    mapping = (; col = :channel),
)

# This difference is the baseline issue: whenever the beta is non-zero beyond the baseline window, classical baseline correction subtracts signal from the whole epoch, while baseline regression estimates (and thereby removes) exactly this co-variance per channel - keeping the effects of interest conditional on the baseline.
