# # [Baseline regression](@id baseline_regression)
# Classical **baseline correction** subtracts the mean of a pre-stimulus window from each epoch. This assumes the baseline is pure noise - but if the baseline co-varies with the epoch (through slow drifts, noise correlations, or the overlap of the previous trial), correcting removes that shared variance from the whole epoch and can bias the estimated effects.
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
using Random

# ```@raw html
# </details >
# ```

# # Simulate data with a slow drift
# We simulate multichannel EEG with `UnfoldSim.jl` - deliberately with only 200 events and little noise - and add a low-frequency oscillation to the signal, with channel-specific amplitude and phase. The simulation projects the components onto all 227 channels; for this tutorial we keep three example channels:
data, evts = UnfoldSim.predef_eeg(;
    multichannel = true, n_repeats = 10, overlap = (2, 0.2), noiselevel = 0.1)
data = data[[130, 184, 50], :]

sfreq = 100
t = (0:size(data, 2)-1) ./ sfreq
amps = [5.0, 4.0, 3.0]
phases = 2π .* rand(MersenneTwister(1), 3)
for ch in 1:3
    data[ch, :] .+= amps[ch] .* sin.(2π .* 0.5 .* t .+ phases[ch])
end

# We drop the last events, whose epochs would run past the end of the recording, and cut the continuous data into epochs:
evts = evts[evts.latency .<= size(data, 2) - 150, :]
data_e, times = Unfold.epoch(data = data, tbl = evts, τ = (-0.2, 1.5), sfreq = sfreq)

# # Baseline regression
# For each channel and epoch, we compute the mean baseline and enter it into the event table as one value per channel per row:
bsl_ix = findall(t -> t <= 0.0, times)
evts.bsl = [vec(mean(data_e[:, bsl_ix, i]; dims = 2)) for i in axes(data_e, 3)]

# The formula marks `bsl` as channel-varying, and the mass-univariate model is fit (currently only supported for epoched, not for time-expanded models):
f = @formula 0 ~ 1 + condition + bychannel(bsl)
m = fit(UnfoldLinearModel, [Any => (f, times)], evts, data_e)

# # The beta of the baseline over time
ct = coeftable(m)
plot_erp(ct[ct.coefname.=="bsl", :]; mapping = (; col = :channel))

# How to read this plot:
# - Within the baseline window, the beta averages exactly 1 for every channel - by construction, the covariate *is* the window mean.
# - Beyond the window, the beta follows the oscillation: it is large in absolute value wherever the drift is in phase with its baseline value, and flips sign where the phase is opposite. A slow oscillation therefore co-varies with the **entire** epoch - the baseline is not just noise, and this plot makes its reach visible.

# # The baseline issue
# What does this mean for the classical approach? We baseline-correct the same epochs (subtract the baseline mean) and fit a plain model:
data_corr = data_e .- mean(data_e[:, bsl_ix, :]; dims = 2)
f_plain = @formula 0 ~ 1 + condition
m_corr = fit(UnfoldLinearModel, [Any => (f_plain, times)], evts, data_corr)

# Comparing the estimated `condition` effects, baseline correction is contaminated by the drift: with only 200 events, the drift does not average out between the two conditions, so subtracting a condition-specific baseline leaves an oscillation-shaped spurious effect. Baseline regression, which estimates the drift's influence per channel instead, keeps the effect flat:
ct_corr = coeftable(m_corr)
cmp = vcat(ct_corr, ct)
cmp.baseline = vcat(fill("correction", nrow(ct_corr)), fill("regression", nrow(ct)))
plot_erp(
    cmp[startswith.(cmp.coefname, "condition"), :];
    mapping = (; col = :channel, color = :baseline),
)

# This is the baseline issue: whenever the baseline co-varies with the epoch, classical baseline correction subtracts signal along with the noise, and with few events the leftover drift enters the effects. Baseline regression instead estimates (and thereby removes) exactly this co-variance per channel - keeping the effects of interest conditional on the baseline.
