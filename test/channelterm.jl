# Tests for channel-varying (bychannel) predictors, e.g. channel-wise baseline regression
# à la Alday et al. (2021)

using LinearAlgebra

rng = StableRNG(42)
nch, ntime, ntr = 4, 50, 100
times = range(-0.2, 1.5, length = ntime)
bsl_ix = findall(t -> t <= 0.0, collect(times))

bsl_amp = [1.0, 2.0, 3.0, 4.0] # per-channel baseline amplitude
effect_amp = 2.5 # common effect

evts = DataFrame(:effect => randn(rng, ntr))
evts.bsl = [bsl_amp .* (1 .+ 0.1 .* randn(rng, nch)) for _ = 1:ntr] # trial-wise noisy baseline per channel

data = zeros(nch, ntime, ntr)
for (i, r) in enumerate(eachrow(evts))
    for ch = 1:nch
        data[ch, :, i] .= r.bsl[ch] .+ effect_amp .* r.effect .+ 0.2 .* randn(rng, ntime)
    end
end

f = @formula 0 ~ 1 + effect + bychannel(bsl)
m = fit(UnfoldLinearModel, [Any => (f, times)], evts, data; show_progress = false)

# the per-channel designmatrices, as used during fitting
X3d = Unfold.modelcols_bychannel(designmatrix(m)[1].formula.rhs, designmatrix(m)[1].events)

@testset "bychannel: fit" begin
    @test Unfold.get_coefnames(m) == ["(Intercept)", "effect", "bsl"]
    @test size(Unfold.coef(m)) == (nch, ntime, 3)
    X = modelmatrix(m)
    @test length(X) == 1
    @test size(X[1]) == (ntr, 3)

    # the stored designmatrix corresponds to the first channel
    @test X[1] == X3d[1, :, :]

    # per-channel designs: shared predictors identical across channels, bychannel ones not
    @test all(X3d[1, :, 1:2] .== X3d[2, :, 1:2]) # intercept + effect
    @test !all(X3d[1, :, 3] .== X3d[2, :, 3]) # bsl differs

    # estimates are identical to fitting each channel separately with its own designmatrix
    est_manual = similar(Unfold.coef(m))
    for ch = 1:nch
        Xch = X3d[ch, :, :]
        Y = permutedims(reshape(@view(data[ch:ch, :, :]), ntime, ntr))
        est_manual[ch, :, :] .= permutedims(pinv(Xch) * Y, (2, 1))
    end
    @test maximum(abs.(Unfold.coef(m) .- est_manual)) < 1e-10

    # recovers the ground truth
    beta_post = Unfold.coef(m)[:, 30, :]
    @test all(isapprox.(beta_post[:, 2], effect_amp; atol = 0.1)) # effect
    @test all(isapprox.(beta_post[:, 3], 1.0; atol = 0.1)) # bsl slope ≈ 1
end

@testset "bychannel: interaction" begin
    f2 = @formula 0 ~ 1 + effect * bychannel(bsl)
    m2 = fit(UnfoldLinearModel, [Any => (f2, times)], evts, data; show_progress = false)
    @test Unfold.get_coefnames(m2) == ["(Intercept)", "effect", "bsl", "effect & bsl"]
    @test size(Unfold.coef(m2)) == (nch, ntime, 4)

    # the interaction is 0 in the simulated data (noise-scale tolerance)
    @test all(isapprox.(Unfold.coef(m2)[:, 30, 4], 0.0; atol = 0.3))
end

@testset "bychannel: effects" begin
    # bychannel predictors not in the reference-grid are fixed to 0,
    # i.e. effects are conditional on them
    eff = Unfold.effects(Dict(:effect => [0.0, 1.0]), m)
    @test size(eff, 1) == nch * 2 * ntime
    @test :channel ∈ propertynames(eff)
    # conditional on baseline (bsl = 0): only intercept + effect
    yhat_eff0 = eff[(eff.effect .== 0.0) .& (eff.time .== times[30]), :]
    @test yhat_eff0.yhat ≈ Unfold.coef(m)[:, 30, 1] # intercept per channel

    # a bychannel predictor in the reference-grid: its values are ignored (with a warning),
    # effects are computed conditional on it (fixed to 0)
    eff = @test_logs (:warn, r"bychannel") Unfold.effects(
        Dict(:effect => [0.0], :bsl => [1.0]),
        m,
    )
    yhat_bsl = eff[(eff.effect .== 0.0) .& (eff.time .== times[30]), :]
    @test yhat_bsl.yhat ≈ Unfold.coef(m)[:, 30, 1] # intercept per channel (bsl fixed to 0)
end

@testset "bychannel: downstream" begin
    @test size(Unfold.stderror(m)) == (nch, ntime, 3)
    ct = Unfold.coeftable(m)
    @test size(ct, 1) == nch * ntime * 3
    @test unique(ct.channel) == 1:nch
    # vif is computed on the stored (first-channel) designmatrix
    v = Unfold.vif(m)
    @test names(v) == ["eventname", "coefname", "VIF"]
    @test size(v, 1) == 3
end

@testset "bychannel: errors" begin
    # bychannel() outside a formula
    @test_throws ErrorException bychannel(1)

    # channel-varying predictors are not supported for time-expanded models
    @test_throws ErrorException fit(
        UnfoldLinearModelContinuousTime,
        [Any => (f, firbasis(τ = (-0.1, 0.5), sfreq = 10))],
        evts,
        rand(rng, nch, 600),
    )

    # channel-vectors must have the same length in every row
    evts_bad = copy(evts)
    evts_bad.bsl = vcat(evts_bad.bsl[1:(end-1)], [randn(rng, nch + 1)])
    @test_throws ArgumentError fit(
        UnfoldLinearModel,
        [Any => (f, times)],
        evts_bad,
        data;
        show_progress = false,
    )

    # predict is not implemented for bychannel models (and residuals, which build on it)
    @test_throws ErrorException Unfold.predict(m)
    @test_throws ErrorException Unfold.residuals(m, data)
end

@testset "bychannel: save/load" begin
    f_tmp = tempname() * ".jld2"
    Unfold.FileIO.save(f_tmp, m)
    m_loaded = Unfold.FileIO.load(f_tmp, UnfoldLinearModel)
    @test maximum(abs.(Unfold.coef(m_loaded) .- Unfold.coef(m))) ≈ 0.0
    @test size(modelmatrix(m_loaded)[1]) == (ntr, 3)
    @test modelmatrix(m_loaded)[1] ≈ modelmatrix(m)[1]
    eff_l = Unfold.effects(Dict(:effect => [0.0]), m_loaded)
    eff_o = Unfold.effects(Dict(:effect => [0.0]), m)
    @test eff_l.yhat ≈ eff_o.yhat
end
