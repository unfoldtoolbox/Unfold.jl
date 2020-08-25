using Test, StatsModels
using DataFrames

using unfold
include("test_utilities.jl")

data,evts = loadtestdata("testCase1") #
##

f_spl  = @formula 0~1+spl(continuousA,4) # 1
f  = @formula 0~1+continuousA # 1
data_r = reshape(data,(1,:))
data_e,times = unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1.),sfreq=10) # cut the data into epochs

m_mul,m_mul_results = fit(UnfoldLinearModel,f,evts,data_e,times)
m_mul_spl,m_mul_results_spl = fit(UnfoldLinearModel,f_spl,evts,data_e,times)

# asking for 4 splines should generate 4 splines 
@test length(unique(m_mul_results_spl.term)) == 5 # XXX check back with unfold whether this is the same! could be n-1 splines in unfold. We should keep that comparable I guess


basisfunction = firbasis(τ=(-1,1),sfreq=10,name="A")
m_tul,m_tul_results = fit(UnfoldLinearModel,f,evts,data_r,basisfunction)
m_tul_spl,m_tul_results_spl = fit(UnfoldLinearModel,f_spl,evts,data_r,basisfunction)
newdf = DataFrame([:continuousA=>collect(1:100.)])

# results from timeexpanded and non should be equal
@test all(m_mul_results_spl.estimate .≈ m_tul_results_spl.estimate)

if 1 == 0
    using Makie
    using StatsMakie
a = lines(Data(unfold.predict(m_tul_spl,newdf)),Group(:continuousA),:times,:yhat)
b = lines(Data(unfold.predict(m_tul,newdf)),Group(:continuousA),:times,:yhat)
# visual check whether linear & spl look similar
vbox(a,b)
end

# test much higher number of splines
f_spl  = @formula 0~1+spl(continuousA,131) # 1
m_mul_spl,m_mul_results_spl = fit(UnfoldLinearModel,f_spl,evts,data_e,times)
@test length(unique(m_mul_results_spl.term)) == 132


