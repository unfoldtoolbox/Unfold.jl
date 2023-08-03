using BenchmarkTools
using Random
using UnfoldSim
using Unfold
using DataFrames
using CategoricalArrays
const SUITE = BenchmarkGroup()
Random.seed!(3)





## Preparatory code
#include("../test/test_utilities.jl"); # to load the simulated data

#data,evts = loadtestdata("testCase6",dataPath=(@__DIR__)*"/../test/data") #

sfreq = 100
data,evts = UnfoldSim.predef_2x2(;n_subjects=20,signalsize=sfreq)
data_epochs,evts_epochs = UnfoldSim.predef_2x2(;n_subjects=20,signalsize=sfreq,return_epoched=true)
transform!(evts, :subject => categorical => :subject)
transform!(evts, :item => categorical => :item)
evts.type = rand(MersenneTwister(1),[0,1],nrow(evts))

ba1 = firbasis(τ=(0,1),sfreq = sfreq,name="evts1")
ba2 = firbasis(τ=(0,1),sfreq = sfreq,name="evts2")

f1_lmm  = @formula 0~1+A+(1+A|subject)
f2_lmm  = @formula 0~1+A+(1+A|item)

f1  = @formula 0~1+A
f2  = @formula 0~1+B

dict_lin = Dict(0 => (f1,ba1),1 => (f2,ba2))
dict_lmm = Dict(0 => (f1_lmm,ba1),1 => (f2_lmm,ba2))


X1_lmm  = designmatrix(UnfoldLinearMixedModel,f1_lmm,evts,ba1)
X2_lmm  = designmatrix(UnfoldLinearMixedModel,f2_lmm,evts,ba2)

X1   = designmatrix(UnfoldLinearModelContinuousTime,f1,evts,ba1)
X2   = designmatrix(UnfoldLinearModelContinuousTime,f2,evts,ba2)


times = 1:size(data_epochs,1)

SUITE = BenchmarkGroup()
SUITE["nodc"] = BenchmarkGroup(
    ["nodc"])
SUITE["dc"] = BenchmarkGroup(["dc"])

# designmatrix generation
SUITE["dc"]["X_gen_lin"] = @benchmarkable designmatrix(UnfoldLinearModelContinuousTime,$f1,$evts,$ba1)
SUITE["dc"]["X_gen_lmm"] = @benchmarkable designmatrix(UnfoldLinearMixedModelContinuousTime,$f1_lmm,$evts,$ba1)

# designmat concatination
SUITE["dc"]["X_cat_lin"] = @benchmarkable $X1+$X2
SUITE["dc"]["X_cat_lmm"] = @benchmarkable $X1_lmm+$X2_lmm

# Model Fit
SUITE["nodc"]["fit_lin"] = @benchmarkable fit(UnfoldModel,$f1,$evts_epochs,$data_epochs,$times)
SUITE["nodc"]["fit_lmm"] = @benchmarkable fit(UnfoldModel,$f1_lmm,$evts_epochs,$data_epochs,$times)
SUITE["dc"]["fit_lin"]   = @benchmarkable fit(UnfoldModel,$dict_lin,$evts,$data,eventcolumn="condA");
SUITE["dc"]["fit_lmm"] = @benchmarkable   fit(UnfoldModel,$dict_lmm,$evts,$data,eventcolumn="condA");