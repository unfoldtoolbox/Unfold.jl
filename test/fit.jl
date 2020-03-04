using Revise
using Test, StatsModels
using DataFrames

import unfold
include("test_utilities.jl")

data,evts = loadtestdata("test/testCase4") #
f  = @formula 0~1+conditionA+conditionB # 4
f = @formula 0~1+continuousA+continuousB # 1

data_r = reshape(data,(1,:))
# cut the data into epochs
data_e,times = unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1.9),sfreq=10)
# Additional step because the epoching function might return missing values
evts_e,data_e = unfold.dropMissingEpochs(evts,data_e)

# Mass Univariate Linear
m_mul = unfold.fit(unfold.UnfoldLinearModel,f,evts,data_e,times)
@test all(m_mul.results[(m_mul.results.time.==0.1),:estimate] .≈ [3.0 2.5 -1.5]')

# Timexpanded Univariate Linear
basisfunction = unfold.firbasis(τ=(-1,1.9),sfreq=100)
m_tul = unfold.fit(unfold.UnfoldLinearModel,f,evts,data,basisfunction)
@test all(m_tul.results[(m_tul.results.time.==0.1),:estimate] .≈ [3.0 2.5 -1.5]')


@time unfold.generateDesignmatrix(unfold.UnfoldLinearModel,f,evts,basisfunction)
# new version 2.7s

###############################
##  Mixed Model tests
###############################
data,evts = loadtestdata("test/testCase3") #
data = data.+ 1*randn(size(data)) # we have to add minimal noise, else mixed models crashes.
categorical!(evts,:subject)
f  = @formula 0~1+condA+condB + (1+condA+condB|subject)
f  = @formula 0~1 + (1|subject)


data_r = reshape(data,(1,:))
# cut the data into epochs
# TODO This ignores subject bounds
data_e,times = unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1.9),sfreq=10)
# Additional step because the epoching function might return missing values
evts_e,data_e = unfold.dropMissingEpochs(evts,data_e)


# Mass Univariate Mixed
@time m_mum = unfold.fit(unfold.UnfoldLinearMixedModel,f,evts,data_e    ,times,contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()))
#@test all(m_mul.results[(m_mul.results.time.==0.1),:estimate] .≈ [3.0 2.5 -1.5]')
plot(m_mum.results.time,m_mum.results.estimate,group=(m_mum.results.term,m_mum.results.group),layout=2,legend=:outerbottom)

# Timexpanded Univariate Mixed
basisfunction = unfold.firbasis(τ=(-0.2,0.3),sfreq=10)
@time m_tum = unfold.fit(unfold.UnfoldLinearMixedModel,f,evts,data,basisfunction, contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()) )
#2@test all(m_tul.results[(m_tul.results.time.==0.1),:estimate] .≈ [3.0 2.5 -1.5]')
plot(m_tum.results.time,m_tum.results.estimate,group=(m_tum.results.term,m_tum.results.group),layout=2,legend=:outerbottom)

# simulation fixef: [10 -2.5 0]
# simulation ranef: [5, 2 3]
##
Xs = unfold.unfoldDesignmatrix(unfold.UnfoldLinearMixedModel,f,evts,basisfunction)

# Fit the model
m = unfold.unfoldFit(unfold.UnfoldLinearMixedModel,Xs,dropdims(data_f,dims=1))

ufresult = unfold.condense(m,evts,times)




##########
