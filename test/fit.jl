# ---
# Test
# ---
using Test, StatsModels
using DataFrames

using Unfold
include("test_utilities.jl")

data,evts = loadtestdata("test_case_3a") #
f  = @formula 0~1+conditionA+continuousA # 1

# prepare data
data_r = reshape(data,(1,:))
data_r = vcat(data_r,data_r)#add second channel

#--------------------------#
## Mass Univariate Linear ##
#--------------------------#
data_e,times = Unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1.9),sfreq=20) # cut the data into epochs
m_mul,m_mul_results = fit(UnfoldLinearModel,f,evts,data_e,times)
@test m_mul_results[(m_mul_results.channel.==1).&(m_mul_results.colname_basis .==0.1),:estimate] ≈ [2,3,4]


data_e_noreshape,times = Unfold.epoch(data=data,tbl=evts,τ=(-1.,1.9),sfreq=20) # cut the data into epochs
m_mul_noreshape,m_mul_results_noreshape = fit(UnfoldLinearModel,f,evts,data_e_noreshape,times)
@test m_mul_results_noreshape[(m_mul_results_noreshape.channel.==1).&(m_mul_results_noreshape.colname_basis .==0.1),:estimate] ≈ [2,3,4]
@test size(m_mul_results_noreshape)[1] ==size(m_mul_results)[1]/2

# Add Missing in Data
data_e_missing = data_e
data_e_missing[1,25,end-5:end] .= missing
m_mul_missing,m_mul_missing_results = Unfold.fit(UnfoldLinearModel,f,evts,data_e_missing,times)
@test m_mul_missing_results.estimate ≈ m_mul_results.estimate

# Special solver solver_lsmr_se with Standard Error
se_solver = solver=(x,y)->Unfold.solver_default(x,y,stderror=true)
m_mul_se,m_mul_results_se = Unfold.fit(UnfoldLinearModel,f,evts,data_e,times,solver=se_solver)
@test all(m_mul_results_se.estimate .≈ m_mul_results.estimate)
@test !all(isnothing.(m_mul_results_se.stderror ))


#---------------------------------#
## Timexpanded Univariate Linear ##
#---------------------------------#
basisfunction = firbasis(τ=(-1,1),sfreq=20,name="basisA")
m_tul,m_tul_results = fit(UnfoldLinearModel,f,evts,data_r,basisfunction)
@test isapprox(m_tul_results[(m_tul_results.channel.==1).&(m_tul_results.colname_basis .==0.1),:estimate],[2,3,4],atol=0.0001)

# test without reshape, i.e. 1 channel vector e.g. size(data) = (1200,)
m_tul_noreshape,m_tul_results_noreshape = fit(UnfoldLinearModel,f,evts,data,basisfunction);
@test size(m_tul_results_noreshape)[1] ==size(m_tul_results)[1]/2

# Test under missing data
data_missing = Array{Union{Missing,Number}}(undef,size(data_r))
data_missing .= data_r
data_missing[4500:4600] .= missing

m_tul_missing,m_tul_missing_results = fit(UnfoldLinearModel,f,evts,data_missing,basisfunction)
@test  isapprox(m_tul_missing_results.estimate , m_tul_results.estimate,atol=1e-4)  # higher tol because we remove stuff


## Test multiple basisfunctions
b1 = firbasis(τ=(-1,1),sfreq=20,name="basisA")
b2 = firbasis(τ=(-1,1),sfreq=20,name="basisB")

f1  = @formula 0~1+continuousA # 1
f2  = @formula 0~1+continuousA # 1

# Fast-lane new implementation
m,res = fit(UnfoldLinearModel,Dict(0=>(f1,b1),1=>(f1,b1)),evts,data,eventcolumn="conditionA")

# slow manual
X1 = designmatrix(UnfoldLinearModel,f1,filter(x->(x.conditionA==0),evts),b1)
X2 = designmatrix(UnfoldLinearModel,f2,filter(x->(x.conditionA==1),evts),b2)

@time m = unfoldfit(UnfoldLinearModel,X1+X2,data')
tmp = condense_long(m)

# test fast way & slow way to be identical
@test all(tmp.estimate .== res.estimate)



# runntime tests - does something explode?
for k in 1:4
    local f
    if k == 1
        f  = @formula 0~1
    elseif k == 2
        f  = @formula 0~1+conditionA
    elseif k == 3
        f  = @formula 0~0+conditionA
    elseif k == 4
        f  = @formula 0~1+continuousA
    end
    println("Testing Runtime $k with Formula:$f")
    Xs = designmatrix(UnfoldLinearModel,f,evts,basisfunction)

    # Fit the model
    df = unfoldfit(UnfoldLinearModel,Xs,data_e)
    c = condense_long(df,times)
    fit(UnfoldLinearModel,f,evts,data_e,times)
    fit(UnfoldLinearModel,f,evts,data,basisfunction)
end

# Special solver solver_lsmr_se with Standard Error
se_solver = solver=(x,y)->Unfold.solver_default(x,y,stderror=true)
m_tul_se,m_tul_results_se = fit(UnfoldLinearModel,f,evts,data_r,basisfunction,solver=se_solver)
@test all(m_tul_results_se.estimate .== m_tul_results.estimate)
@test !all(isnothing.(m_tul_results_se.stderror ))

#m_mul_se,m_mul_results_se = fit(UnfoldLinearModel,f,evts,data_e.+randn(size(data_e)).*5,times,solver=se_solver)
#plot_results(m_mul_results_se[m_mul_results_se.channel.==1,:],se=true)
#m_tul_se,m_tul_results_se = fit(UnfoldLinearModel,f,evts,data_r.+randn(size(data_r)).*5,basisfunction,solver=se_solver)
#plot_results(m_tul_results_se[m_tul_results_se.channel.==1,:],se=true)

##
data_long,evts_long = loadtestdata("test_case_1c") #
data_long = reshape(data_long,(1,:))

data_long = vcat(data_long,data_long)
f_long  = @formula 0~1
basisfunction_long = firbasis(τ=(-1,1),sfreq=1000,name="basisA")

@time designmatrix(UnfoldLinearModel,f_long,evts_long,basisfunction_long)
@time a,m_tul_long = fit(UnfoldLinearModel,f_long,evts_long,data_long,basisfunction_long);

@test isapprox(m_tul_long[(m_tul_long.channel.==1).&(m_tul_long.colname_basis .==0.1),:estimate],[2],atol=0.0001)
# ~21s, test_case_1c, sfreq = 1000, 6000 events
#
## Older numbers for "designmatrix" only. But see "benchmark/benchmarkjl" for better benchmarks
# new version 7s-10s, dataset4, sfreq=1000, 1200stim,
# ~13s, test_case_1c, sfreq = 1000, 6000 events


###############################
##  Mixed Model tests
###############################
data,evts = loadtestdata("testCase3",dataPath=(@__DIR__)*"/data") #
append!(data,zeros(1000))
data = reshape(data,1,:)
data = vcat(data,data)
data = data.+ 1*randn(size(data)) # we have to add minimal noise, else mixed models crashes.
data_missing = Array{Union{Missing,Number}}(undef,size(data))
data_missing .= data
data_missing[4500:4600] .= missing

categorical!(evts,:subject)
f  = @formula 0~1+condA+condB + (1+condA+condB|subject)
#f  = @formula 0~1 + (1|subject)



# cut the data into epochs
# TODO This ignores subject bounds
data_e,times = Unfold.epoch(data=data,tbl=evts,τ=(-1.,1.9),sfreq=10)
data_missing_e,times = Unfold.epoch(data=data_missing,tbl=evts,τ=(-1.,1.9),sfreq=10)
evts_e,data_e = Unfold.dropMissingEpochs(copy(evts),data_e)
evts_missing_e,data_missing_e = Unfold.dropMissingEpochs(copy(evts),data_missing_e)

######################
##  Mass Univariate Mixed
@time m_mum = fit(UnfoldLinearMixedModel,f,evts_e,data_e    ,times,contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()))
df = m_mum[2]
@test isapprox(df[(df.channel .== 1).&(df.term.=="condA: 1").&(df.colname_basis.==0.0),:estimate], [5.618,9.175],rtol=0.1)
@time m_mum = fit(UnfoldLinearMixedModel,f,evts_missing_e,data_missing_e    ,times,contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()))
df = m_mum[2]
@test isapprox(df[(df.channel .== 1).&(df.term.=="condA: 1").&(df.colname_basis.==0.0),:estimate], [5.618,9.175],rtol=0.1)


# Timexpanded Univariate Mixed
f  = @formula 0~1+condA+condB + (1+condA+condB|subject)
basisfunction = firbasis(τ=(-0.2,0.3),sfreq=10,name="ABC")
@time m_tum = fit(UnfoldLinearMixedModel,f,evts,data,basisfunction, contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()) )
df = m_tum[2]
@test isapprox(df[(df.channel .== 1).&(df.term.=="condA: 1").&(df.colname_basis.==0.0),:estimate], [5.618,9.175],rtol=0.1)


# missing data in LMMs
# not yet implemented
Test.@test_broken  m_tum = fit(UnfoldLinearMixedModel,f,evts,data_missing,basisfunction, contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()) )

evts.subjectB = evts.subject;
evts1 = evts[evts.condA.==0,:]
evts2 = evts[evts.condA.==1,:]

f0_lmm  = @formula 0~1+condB+(1|subject) + (1|subjectB)
@time m_tum = fit(UnfoldLinearMixedModel,f0_lmm,evts,data,basisfunction)


f1_lmm  = @formula 0~1+condB+(1|subject)
f2_lmm  = @formula 0~1+condB+(1|subjectB)

b1 = firbasis(τ=(-0.2,0.3),sfreq=10,name="A")
b2 = firbasis(τ=(-0.1,0.3),sfreq=10,name="B")


X1_lmm  = designmatrix(UnfoldLinearMixedModel,f1_lmm,evts1,b1)
X2_lmm  = designmatrix(UnfoldLinearMixedModel,f2_lmm,evts2,b2)

r = Unfold.unfoldfit(UnfoldLinearMixedModel,X1_lmm+X2_lmm,data);
df = condense_long(r)

@test isapprox(df[(df.channel .== 1).&(df.term.=="condB").&(df.colname_basis.==0.0),:estimate],[18.18,10.4],rtol=0.1)

# Fast-lane new implementation
m,res = fit(UnfoldLinearMixedModel,Dict(0=>(f1_lmm,b1),1=>(f2_lmm,b2)),evts,data,eventcolumn="condA")


if 1 == 0
    using WGLMakie,AlgebraOfGraphics
    m = mapping(:colname_basis,:estimate,color=:term,layout_x=:term,layout_y=:basisname)
    
    AlgebraOfGraphics.data(df[df.channel.==2,:]) * visual(Lines) * m  |> draw

end


##------



if 1 == 0
    # Fit mass-univariate 1st level for all subjects
    basisfunction = firbasis(τ=(-.1,.4),sfreq=10,name="A")
    resAll = DataFrame()
    f  = @formula 0~1+condA+condB
    for s in unique(evts.subject)
        ##
        from = minimum(evts[evts.subject.==s,:latency])-10
        to = maximum(evts[evts.subject.==s,:latency])+40
        to = min(to,size(data,1))
        evts_s = evts[evts.subject.==s,:]
        evts_s.latency .-= from
        m = fit(UnfoldLinearModel,f,evts_s,data[from:to],basisfunction,contrasts=Dict(:condA => EffectsCoding(), :condB => EffectsCoding()))
        m.results.subject = s
        append!(resAll,m.results)
    end

    results = resAll[resAll.term.=="(Intercept)",:]

    results[results.time .== .1,:]
end
