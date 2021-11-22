using Test
using Unfold
using StatsModels
using DataFrames
using Statistics
include("test_utilities.jl")

data, evts = loadtestdata("test_case_3a") #

data_r = reshape(data, (1, :))
data_e, times = Unfold.epoch(data = data_r, tbl = evts, τ = (0, 0.05), sfreq = 10) # cut the data into epochs

##
f = @formula 0 ~ 1 + conditionA + continuousA # 1
m_mul = fit(Unfold.UnfoldModel, Dict(Any=>(f,times)), evts, data_e)

# test simple case
eff = Unfold.effects(Dict(:conditionA => [0,1],:continuousA =>[0]),m_mul)
@test size(eff,1) == 2 # we specified 2 levels @ 1 time point
@test eff.conditionA ≈ [0.,1.] # we want to different levels
@test eff.yhat ≈ [2.0,5.0] # these are the perfect predicted values

# combination 2 levels /  6 values
eff = Unfold.effects(Dict(:conditionA => [0,1],:continuousA =>[-2,0,2]),m_mul)
@test size(eff,1) == 6 # we want 6 values
@test eff.conditionA ≈ [0.,0.,0.,1.,1.,1.] 
@test eff.continuousA ≈ [-2,0,2,-2,0,2.] 

# testing typical value
eff_man = Unfold.effects(Dict(:conditionA => [0,1],:continuousA =>[mean(evts.continuousA)]),m_mul)
eff_typ = Unfold.effects(Dict(:conditionA => [0,1]),m_mul)
@test eff_man.yhat ≈ eff_typ.yhat


## Testing Splines
f_spl = @formula 0 ~ 1 + conditionA + spl(continuousA, 4) # 1
m_mul_spl = fit(UnfoldModel, f_spl, evts, data_e, times)

eff = Unfold.effects(Dict(:conditionA => [0,1],:continuousA =>[0]),m_mul_spl)
@test size(eff,1) == 2 # we specified 2 levels @ 1 time point
@test eff.conditionA ≈ [0.,1.] # we want to different levels
@test eff.yhat ≈ [2.0,5.0] # these are the perfect predicted values

# combination 2 levels /  6 values
eff = Unfold.effects(Dict(:conditionA => [0,1],:continuousA =>[-0.5,0,0.5]),m_mul_spl)
@test size(eff,1) == 6 # we want 6 values
@test eff.conditionA ≈ [0.,0.,0.,1.,1.,1.] 
@test eff.continuousA ≈ [-0.5,0,0.5] 
@test eff.yhat ≈ [0,2,4,3,5,7]

# testing for safe predictions
eff = Unfold.effects(Dict(:conditionA => [0,1],:continuousA =>[2]),m_mul_spl)
@test all(ismissing.(eff.yhat ))


## Timeexpansion

f = @formula 0 ~ 1 + conditionA + continuousA # 1
uf = fit(Unfold.UnfoldModel, Dict(Any=>(f,firbasis([0,0.05],10))), evts, data_e)
eff = Unfold.effects(Dict(:conditionA => [0,1],:continuousA =>[0]),uf)
