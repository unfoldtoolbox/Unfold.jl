using BenchmarkTools
using Random
using StableRNGs
using UnfoldSim
using Unfold
using UnfoldMixedModels
using DataFrames
using CategoricalArrays
using MixedModels
using BSplineKit
#const SUITE = BenchmarkGroup()
#Random.seed!(3)


#---


sfreq = 100
data_multsub, evts_multsub =
    UnfoldSim.predef_2x2(StableRNG(1); n_subjects = 20, signalsize = sfreq)
data_multsub_epochs, evts_multsub_epochs = UnfoldSim.predef_2x2(
    StableRNG(1);
    n_subjects = 20,
    signalsize = sfreq,
    return_epoched = true,
)
data, evts = UnfoldSim.predef_2x2(StableRNG(1); n_subjects = 1, signalsize = sfreq)
data_epochs, evts_epochs = UnfoldSim.predef_2x2(
    StableRNG(1);
    n_subjects = 1,
    signalsize = sfreq,
    return_epoched = true,
)


data_multsub = data_multsub[:]
data_multsub_epochs = reshape(data_multsub_epochs, size(data_multsub_epochs, 1), :)
transform!(evts_multsub, :subject => categorical => :subject)
transform!(evts_multsub, :item => categorical => :item)
evts.type = rand(StableRNG(1), ["A", "B"], nrow(evts))
evts_multsub.type = rand(StableRNG(1), ["A", "B"], nrow(evts_multsub))

for (ix, k) in
    enumerate([:continuousA, :continuousB, :continuousC, :continuousD, :continuousE])
    evts[!, k] = rand(StableRNG(ix), size(evts, 1))
    evts_multsub[!, k] = rand(StableRNG(ix), size(evts_multsub, 1))
    evts_epochs[!, k] = rand(StableRNG(ix), size(evts_epochs, 1))
    evts_multsub_epochs[!, k] = rand(StableRNG(ix), size(evts_multsub_epochs, 1))
end

ba1 = firbasis(τ = (-0.1, 1), sfreq = sfreq, name = "A") # names explicitly given due to benchmark?
ba2 = firbasis(τ = (-0.2, 1), sfreq = sfreq, name = "B")

f1 = @formula 0 ~ 1 + A
f1_spl = @formula 0 ~
         1 +
         A +
         spl(continuousA, 5) +
         spl(continuousB, 5) +
         spl(continuousC, 5) +
         spl(continuousD, 5) +
         spl(continuousE, 5)
f2 = @formula 0 ~ 1 + B

f1_lmm = @formula 0 ~ 1 + A + (1 + A | subject)
f2_lmm = @formula 0 ~ 1 + A + (1 + A | item)
dict_lin = Dict("A" => (f1, ba1), "B" => (f2, ba2))
dict_spl = Dict("A" => (f1_spl, ba1), "B" => (f1_spl, ba2))
dict_lmm = Dict("A" => (f1_lmm, ba1), "B" => (f2_lmm, ba2))

times = 1:size(data_epochs, 1)



#---
m_epoch_lin_f1 = fit(UnfoldModel, f1, evts_epochs, data_epochs, times)
m_epoch_lin_f1_spl = fit(UnfoldModel, f1_spl, evts_epochs, data_epochs, times)

m_lin_f1 = fit(UnfoldModel, dict_lin, evts, data, eventcolumn = "type")


m_lin_f1_spl = fit(UnfoldModel, dict_spl, evts, data, eventcolumn = "type")

if 1 == 0
    m_lin_f1_spl_ch =
        fit(UnfoldModel, dict_spl, evts, repeat(data, 1, 100)', eventcolumn = "type")
end
#---


SUITE = BenchmarkGroup()
SUITE["designmat"] = BenchmarkGroup(["designmat"])
SUITE["fit"] = BenchmarkGroup(["fit"])
SUITE["effects"] = BenchmarkGroup(["effects"])

# designmatrix generation
SUITE["designmat"]["lin"] =
    @benchmarkable designmatrix(UnfoldLinearModelContinuousTime, $f1, $evts, $ba1)
SUITE["designmat"]["lmm"] = @benchmarkable designmatrix(
    UnfoldMixedModels.UnfoldLinearMixedModelContinuousTime,
    $f1_lmm,
    $evts_multsub,
    $ba1,
)

# Model Fit
SUITE["fit"]["lin"] =
    @benchmarkable fit(UnfoldModel, $f1, $evts_epochs, $data_epochs, $times)
SUITE["fit"]["lmm"] = @benchmarkable fit(
    UnfoldModel,
    $f1_lmm,
    $evts_multsub_epochs,
    $data_multsub_epochs,
    $times,
)
SUITE["fit"]["lin_deconv"] =
    @benchmarkable fit(UnfoldModel, $dict_lin, $evts, $data, eventcolumn = "type");
#SUITE["fit"]["lmm_deconv"] =
#    @benchmarkable fit(UnfoldModel, $dict_lmm, $evts, $data, eventcolumn = "type");


SUITE["effects"]["lin"] =
    @benchmarkable effects($(Dict("A" => ["a_small", "a_big"])), $m_lin_f1)
SUITE["effects"]["lin_spl"] = @benchmarkable effects(
    $(Dict(:continuousA => collect(range(0.1, 0.9, length = 15)))),
    $m_lin_f1_spl,
)

#SUITE["fit"]["debugging"] = @benchmarkable read(run(`git diff Project.toml`))

if 1 == 0
    @benchmark designmatrix(UnfoldLinearModelContinuousTime, f1_spl, evts, ba1)

    ba1 = firbasis(τ = (-0.1, 1), sfreq = sfreq)
    ba1_int = firbasis(τ = (-0.1, 1), sfreq = sfreq, interpolate = true)
    using BenchmarkTools
    @benchmark X = designmatrix(UnfoldLinearModelContinuousTime, f1_spl, evts, ba1)
    @benchmark X2 = designmatrix(UnfoldLinearModelContinuousTime, f1_spl, evts, ba1_int)



end
