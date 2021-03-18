using BenchmarkTools
using unfold,StatsModels,MixedModels,DataFrames
using Random: randn,seed!


const SUITE = BenchmarkGroup()
seed!(3)





# Preparatory code
include("../test/test_utilities.jl"); # to load the simulated data

data,evts = loadtestdata("testCase6",dataPath=(@__DIR__)*"/../test/data") #

categorical!(evts,:subject);
categorical!(evts,:stimulus);
evts.subjectB = evts.subject;
evts1 = evts[evts.condA.==0,:]
evts2 = evts[evts.condA.==1,:]
data = data.+ 0.5*randn(size(data)) # we have to add minimal noise, else mixed models crashes.


f1_lmm  = @formula 0~1+condB+(1|subject)
f2_lmm  = @formula 0~1+condB+(1|subjectB)

f1  = @formula 0~1+condB
f2  = @formula 0~1+condB

data_r = reshape(data,(1,:)) #make ch x time
data_epochs,times = unfold.epoch(data=data_r,tbl=evts,τ=(-0.4,0.8),sfreq=100); #epoch
evts,data_epochs = unfold.dropMissingEpochs(evts,data_epochs) # remove missing

ba1 = firbasis(τ=(0,1),sfreq = 10,name="evts1")
ba2 = firbasis(τ=(0,1),sfreq = 10,name="evts2")

X1_lmm  = designmatrix(UnfoldLinearMixedModel,f1_lmm,evts1,ba1)
X2_lmm  = designmatrix(UnfoldLinearMixedModel,f2_lmm,evts2,ba2)

X1   = designmatrix(UnfoldLinearModel,f1,evts1,ba1)
X2   = designmatrix(UnfoldLinearModel,f2,evts2,ba2)




SUITE = BenchmarkGroup()
SUITE["nodc"] = BenchmarkGroup(
    ["nodc"])
SUITE["dc"] = BenchmarkGroup(["dc"])
# designmatrix generation
SUITE["dc"]["X_gen_lin"] = @benchmarkable designmatrix(UnfoldLinearModel,$f1,$evts1,$ba1)
SUITE["dc"]["X_gen_mix"] = @benchmarkable designmatrix(UnfoldLinearMixedModel,$f1_lmm,$evts1,$ba1)

# designmat concatination
SUITE["dc"]["X_cat_lin"] = @benchmarkable $X1+$X2
SUITE["dc"]["X_cat_mix"] = @benchmarkable $X1_lmm+$X2_lmm

# Model Fit
SUITE["nodc"]["fit_lin"] = @benchmarkable unfold.fit(UnfoldLinearModel,$f1,$evts1,$data_epochs,$times)
#SUITE["nodc"]["fit_mix"] = @benchmarkable unfold.fit(UnfoldLinearMixedModel,$f1_lmm,$evts1,$data_epochs,$times)
SUITE["dc"]["fit_lin"] = @benchmarkable unfold.unfoldfit(UnfoldLinearModel,$X1+$X2,$data_r);
SUITE["dc"]["fit_mix"] = @benchmarkable unfold.unfoldfit(UnfoldLinearMixedModel,$X1_lmm+$X2_lmm,$data_r);
