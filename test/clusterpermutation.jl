using Test
using MixedModelsSim
using Random
using DSP
include("test_utilities.jl")
τ = 1.5
fs = 18
evt, epoch_dat = simulate_lmm(
    MersenneTwister(1),
    τ,
    fs,
    β = [0.0, -1.0],
    σs = [1, 1, 1],
    σ = 1,
    n_sub = 20,
    n_item = 30,
    noise_type = "AR-exponential",
);

#mres,res = fit(UnfoldLinearMixedModel,f,evt,dat[chIx:chIx,:,:].*1e6 ,times,contrasts=Dict(:category => EffectsCoding(), :condition => EffectsCoding()));

times = range(0, stop = τ - 1 / fs, step = 1 / fs)

f2 = @formula 0 ~ 1 + stimType + (1 + stimType | subj) + (1 | item)
mres, res = fit(UnfoldLinearMixedModel, f2, evt, epoch_dat, times)

p_df = cluster_permutation_test(mres, epoch_dat, times, 2)

@test p_df[1, :pval] == 0.036

#--- Testing cluster_permutation call
tRange = 1:length(times)
nPerm = 10
permDat = cluster_permutation(mres, epoch_dat, tRange, 2, nPerm)
@test size(permDat, 2) == nPerm
@test size(permDat, 1) == length(tRange)

tRange = 5:length(times)-5
permDat = cluster_permutation(mres, epoch_dat, tRange, 2, nPerm)
@test size(permDat, 1) == length(tRange)


#------ testing pymne_cluster
z_obs = [
    m.z for m in coefpvalues(mres.modelinfo) if
    String(m.coefname) == coefnames(mres.X.formulas)[2][coeffOfInterest]
]

#-- clusterpermutation
obs_cluster = pymne_cluster(z_obs, 4)
@test all(abs.(obs_cluster[2]) .> 4) # clusters have to be larger 4
@test obs_cluster[1][1][1].start == 13 # start of cluster
@test obs_cluster[1][1][1].stop == 14 # end of cluster

obs_cluster = pymne_cluster(z_obs, 0.1)
@test all(abs.(obs_cluster[2]) .> 0.1) # all clusters have to be larger 0.1
@test length(obs_cluster[2]) == 12 # if nothing changed, 12 clusters should be found

#-- tfce
clusterFormingThreshold = Dict(:start => 0, :step => 0.2)
obs_cluster = pymne_cluster(z_obs, clusterFormingThreshold)
@test all(obs_cluster[2] .>= 0) # TFCE is always positive
@test findmax(obs_cluster[2])[2] == findmax(abs.(z_obs))[2] # tfce should not move maximum

#gen_han(τ,fs,3)
if 1 == 0
    # plot some design bits
    map = mapping(:subj, :dv, color = :stimType)#,linestyle=:channel)
    AlgebraOfGraphics.data(evt) * visual(Scatter) * map |> draw
end
