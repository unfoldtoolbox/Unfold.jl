using Test, StatsModels
using DataFrames

import unfold
include("test_utilities.jl")

data,evts = loadtestdata("testCase1") #
f  = @formula 0~1+continuousA+continuousB # 1

data_r = reshape(data,(1,:))
# cut the data into epochs
data_e,times = unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1.9),sfreq=10)
# Additional step because the epoching function might return missing values
evts_e,data_e = unfold.dropMissingEpochs(evts,data_e)

# Mass Univariate Linear
m_mul = unfold.fit(unfold.UnfoldLinearModel,f,evts,data_e,times)
@test all(m_mul.results[(m_mul.results.time.==0.1),:estimate] .≈ [3.0 2.5 -1.5]')

# Timexpanded Univariate Linear
basisfunction = unfold.firbasis(τ=(-1,1),sfreq=10,eventname="A")
m_tul = unfold.fit(unfold.UnfoldLinearModel,f,evts,data,basisfunction)
@test all(m_tul.results[(m_tul.results.time.==0.1),:estimate] .≈ [3.0 2.5 -1.5]')


data4,evts4 = loadtestdata("testCase4") #
f4  = @formula 0~1+conditionA+conditionB # 4
basisfunction4 = unfold.firbasis(τ=(-1,1),sfreq=1000,eventname="A")

@time unfold.generateDesignmatrix(unfold.UnfoldLinearModel,f4,evts4,basisfunction4)
# new version 7s-10s, dataset4, sfreq=1000, 1200stim,

###############################
##  Mixed Model tests
###############################
data,evts = loadtestdata("testCase3") #
data = data.+ 1*randn(size(data)) # we have to add minimal noise, else mixed models crashes.
categorical!(evts,:subject)
f  = @formula 0~1+condA+condB + (1+condA+condB|subject)
#f  = @formula 0~1 + (1|subject)


data_r = reshape(data,(1,:))
# cut the data into epochs
# TODO This ignores subject bounds
data_e,times = unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1.9),sfreq=10)
# Additional step because the epoching function might return missing values
evts_e,data_e = unfold.dropMissingEpochs(evts,data_e)


# Mass Univariate Mixed
@time m_mum = unfold.fit(unfold.UnfoldLinearMixedModel,f,evts,data_e    ,times,contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()))
#@test all(m_mul.results[(m_mul.results.time.==0.1),:estimate] .≈ [3.0 2.5 -1.5]')
#plot(m_mum)

# Timexpanded Univariate Mixed
basisfunction = unfold.firbasis(τ=(-0.2,0.3),sfreq=10,eventname="")
@time m_tum = unfold.fit(unfold.UnfoldLinearMixedModel,f,evts,data,basisfunction, contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()) )

f  = @formula 0~1+(1|subject)
@test_broken  m_tum = unfold.fit(unfold.UnfoldLinearMixedModel,f,evts,data,basisfunction, contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()) )
#2@test all(m_tul.results[(m_tul.results.time.==0.1),:estimate] .≈ [3.0 2.5 -1.5]')
#plot(m_tum.results.time,m_tum.results.estimate,group=(m_tum.results.term,m_tum.results.group),layout=2,legend=:outerbottom)

# simulation fixef: [10 5 10]
# simulation ranef: [3, 3 3]


##########
basisfunction = unfold.firbasis(τ=(-.1,.4),sfreq=10,eventname="A")
resAll = DataFrame()
f  = @formula 0~1+condA+condB
for s in unique(evts.subject)
    ##
    from = minimum(evts[evts.subject.==s,:latency])-10
    to = maximum(evts[evts.subject.==s,:latency])+40
    to = min(to,size(data,1))
    evts_s = evts[evts.subject.==s,:]
    evts_s.latency .-= from
    m = unfold.fit(unfold.UnfoldLinearModel,f,evts_s,data[from:to],basisfunction,contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()))
    m.results.subject = s
    append!(resAll,m.results)
end

results = resAll[resAll.term.=="(Intercept)",:]
#plot(results.time,results.estimate,
#        group=(results.subject),
#        layout=1,legend=true)

results[results.time .== .1,:]
