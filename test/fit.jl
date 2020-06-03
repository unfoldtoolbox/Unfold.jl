##
using Test, StatsModels
using DataFrames

import unfold
include("test_utilities.jl")

data,evts = loadtestdata("testCase1") #
f  = @formula 0~1+continuousA+continuousB # 1

# prepare data
data_r = reshape(data,(1,:))
data_r = vcat(data_r,data_r)#add second channel
data_e,times = unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1.9),sfreq=10) # cut the data into epochs
data_e_missing = data_e
data_e_missing[1,25,end-5:end] .= missing
data_missing = Array{Union{Missing,Number}}(undef,size(data_r))
data_missing .= data_r
data_missing[4500:4600] .= missing

## Mass Univariate Linear

m_mul,m_mul_results = unfold.fit(unfold.UnfoldLinearModel,f,evts,data_e,times)

@test all(m_mul_results[(m_mul_results.channel.==1).&(m_mul_results.colnames_basis .==0.1),:estimate] .≈ [3.0 2.5 -1.5]')
# Timexpanded Univariate Linear
basisfunction = unfold.firbasis(τ=(-1,1),sfreq=10,name="A")
m_tul,m_tul_results = unfold.fit(unfold.UnfoldLinearModel,f,evts,data_r,basisfunction)
@test all(m_tul_results[(m_tul_results.channel.==1).&(m_tul_results.colnames_basis .==0.1),:estimate] .≈ [3.0 2.5 -1.5]')


# Add Missing in Data
m_mul_missing,m_mul_missing_results = unfold.fit(unfold.UnfoldLinearModel,f,evts,data_e_missing,times)
@test m_mul_missing_results.estimate ≈ m_mul_results.estimate
# Timexpanded Univariate Linear
m_tul_missing,m_tul_missing_results = unfold.fit(unfold.UnfoldLinearModel,f,evts,data_missing,basisfunction)
@test  isapprox(m_tul_missing_results.estimate , m_tul_results.estimate,atol=1e-5)  # higher tol because we remove stuff

# runntime tests - does something explode?
for k in 1:3
    if k == 1
        f  = @formula 0~1
    elseif k == 2
        f  = @formula 0~1+continuousA
    elseif k == 3
        f  = @formula 0~1+continuousB
    end
    println("Testing Runtime $k with Formula:$f")
    Xs = unfold.unfoldDesignmatrix(unfold.UnfoldLinearModel,f,evts,basisfunction)

    # Fit the model
    df = unfold.unfoldFit(unfold.UnfoldLinearModel,Xs,data_e)
    c = unfold.condense_long(df,times)
    unfold.fit(unfold.UnfoldLinearModel,f,evts,data_e,times)
    unfold.fit(unfold.UnfoldLinearModel,f,evts,data,basisfunction)
end


##
data4,evts4 = loadtestdata("testCase4") #
data4 = reshape(data4,(1,:))

data4 = vcat(data4,data4)
f4  = @formula 0~1+conditionA+conditionB # 4
basisfunction4 = unfold.firbasis(τ=(-1,1),sfreq=1000,name="A")

@time unfold.generateDesignmatrix(unfold.UnfoldLinearModel,f4,evts4,basisfunction4)
# new version 7s-10s, dataset4, sfreq=1000, 1200stim,

###############################
##  Mixed Model tests
###############################
data,evts = loadtestdata("testCase3") #
append!(data,zeros(1000))
data = reshape(data,1,:)
#data = vcat(data,data)
data = data.+ 1*randn(size(data)) # we have to add minimal noise, else mixed models crashes.
data_missing = Array{Union{Missing,Number}}(undef,size(data))
data_missing .= data
data_missing[4500:4600] .= missing

categorical!(evts,:subject)
f  = @formula 0~1+condA+condB + (1+condA+condB|subject)
#f  = @formula 0~1 + (1|subject)



# cut the data into epochs
# TODO This ignores subject bounds
data_e,times = unfold.epoch(data=data,tbl=evts,τ=(-1.,1.9),sfreq=10)
evts_e,data_e = unfold.dropMissingEpochs(evts,data_e)

######################
##  Mass Univariate Mixed
@time m_mum = unfold.fit(unfold.UnfoldLinearMixedModel,f,evts_e,data_e    ,times,contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()))
#@test all(m_mul_results[(m_mul_results.time.==0.1),:estimate] .≈ [3.0 2.5 -1.5]')
#plot(m_mum)

# Timexpanded Univariate Mixed
basisfunction = unfold.firbasis(τ=(-0.2,0.3),sfreq=10)
@time m_tum = unfold.fit(unfold.UnfoldLinearMixedModel,f,evts,data,basisfunction, contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()) )





import Logging
Logging.global_logger(Logging.SimpleLogger(stdout, Logging.Debug))

basisfunction = unfold.firbasis(τ=(-0.1,.3),sfreq=10)

f  = @formula 0~1+condA+condB # 1
Xs = unfold.unfoldDesignmatrix(unfold.UnfoldLinearModel,f,evts_e)
ufModel_A = unfold.unfoldFit(unfold.UnfoldLinearModel,Xs,data_e)

basisfunction = unfold.firbasis(τ=(-0.1,.3),sfreq=10)
Xs = unfold.unfoldDesignmatrix(unfold.UnfoldLinearModel,f,evts,basisfunction)
basisfunction = unfold.firbasis(τ=(-0.1,.5),sfreq=10)
Xs2 = unfold.unfoldDesignmatrix(unfold.UnfoldLinearModel,f,evts,basisfunction)
#Xs = Xs+Xs2

ufModel_B = unfold.unfoldFit(unfold.UnfoldLinearModel,Xs,data)


f  = @formula 0~1+condA+condB + (1+condA|subject)
Xs = unfold.unfoldDesignmatrix(unfold.UnfoldLinearMixedModel,f,evts_e)

ufModel_C = unfold.unfoldFit(unfold.UnfoldLinearMixedModel,Xs,data_e)

basisfunction = unfold.firbasis(τ=(-0.1,.3),sfreq=10)
Xs = unfold.unfoldDesignmatrix(unfold.UnfoldLinearMixedModel,f,evts,basisfunction)

ufModel_D = unfold.unfoldFit(unfold.UnfoldLinearMixedModel,Xs,data)


ufA = unfold.condense_long(ufModel_A,times)
ufB = unfold.condense_long(ufModel_B)
ufC = unfold.condense_long(ufModel_C,times)
ufD = unfold.condense_long(ufModel_D)

plot(ufB.colnames_basis,ufB.estimate)
plot(ufB.colnames_basis,ufB.estimate)
plot(ufC.colnames_basis,ufC.estimate)
plot(ufD.colnames_basis,ufD.estimate)

##########
if 1 == 0
    # Fit mass-univariate 1st level for all subjects
    basisfunction = unfold.firbasis(τ=(-.1,.4),sfreq=10,name="A")
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

    results[results.time .== .1,:]
end
