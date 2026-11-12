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
# We simulate multichannel EEG with `UnfoldSim.jl` - deliberately with only 200 events and little noise - and add a slow drift to the signal: exponentially correlated noise (`ExponentialNoise`) with channel-specific amplitude. We also present the two conditions in blocks (`event_order_function`), so each condition is recorded against a different level of the drift - like a real block experiment where the environment changes over the session. The simulation projects the components onto all 227 channels; for this tutorial we keep three example channels:
data, evts = UnfoldSim.predef_eeg(;
    multichannel = true, n_repeats = 10, overlap = (2, 0.2), noiselevel = 0.1,
    # conditions in blocks: each condition sees a different level of the drift
    event_order_function = (rng, x) -> sort(x, :condition),
)
data = data[[130, 184, 50], :]

sfreq = 100
amps = [5.0, 4.0, 3.0]
rng = MersenneTwister(1)
for ch in 1:3
    # τ = 250 samples = 2.5 s: comparable to the 1.7 s epochs, so the drift varies within the epoch
    data[ch, :] .+= simulate_noise(
        rng, ExponentialNoise(noiselevel = amps[ch], τ = 250), size(data, 2),
    )
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
# - Beyond the window, the beta decays with the drift's autocorrelation reach - to roughly 0.7 beyond 500 ms. A slow drift therefore co-varies with the **entire** epoch - the baseline is not just noise, and this plot makes its reach visible.

# # The baseline issue
# What does this mean for the classical approach? We baseline-correct the same epochs (subtract the baseline mean) and fit a plain model:
data_corr = data_e .- mean(data_e[:, bsl_ix, :]; dims = 2)
f_plain = @formula 0 ~ 1 + condition
m_corr = fit(UnfoldLinearModel, [Any => (f_plain, times)], evts, data_corr)

# Comparing the estimated `condition` effects: the conditions come in blocks, so each condition is seen against a different drift level - but with a stationary drift these differences average out over the trials of each block, and both approaches estimate similar effects. Neither removes the part of the drift the baseline cannot predict, late in the epoch - which is why baseline correction is no substitute for high-pass filtering. The difference between the methods is what they assume: correction fixes the baseline's influence at β = 1 at every latency, regression estimates it - and the coefficient plot above shows when that assumption holds.
ct_corr = coeftable(m_corr)
cmp = vcat(ct_corr, ct)
cmp.baseline = vcat(fill("correction", nrow(ct_corr)), fill("regression", nrow(ct)))
plot_erp(
    cmp[startswith.(cmp.coefname, "condition"), :];
    mapping = (; col = :channel, color = :baseline),
)

# This is the baseline issue: whenever the baseline co-varies with the epoch, classical baseline correction subtracts signal along with the noise, and assumes its influence is β = 1 at every latency. Baseline regression instead estimates (and thereby removes) exactly this influence per channel - keeping the effects of interest conditional on the baseline.
