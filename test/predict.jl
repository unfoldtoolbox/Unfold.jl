using Test, StatsModels
using DataFrames
using StatsBase

using unfold
include("test_utilities.jl")


data,evts = loadtestdata("testCase1") #
data_e,times = unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1.1),sfreq=10) # cut the data into epochs
basisfunction = firbasis(τ=(-1,1),sfreq=10,name="A")

f  = @formula 0~1 # 1
m_tul,m_tul_results = fit(UnfoldLinearModel,f,evts,data_r,basisfunction)
m_mul,m_mul_results = fit(UnfoldLinearModel,f,evts,data_e,times)

yhat = predict(m_tul,evts)
@test yhat[end-3,:yhat] ≈ mean(data[data.!=0])
yhat2 = predict(m_mul,evts)
@test yhat.yhat ≈ yhat2.yhat
f  = @formula 0~1+continuousA# 1
m_tul,m_tul_results = fit(UnfoldLinearModel,f,evts,data_r,basisfunction)
m_mul,m_mul_results = fit(UnfoldLinearModel,f,evts,data_e,times)
yhat = predict(m_tul,evts)
@test predict(m_mul,evts).yhat ≈ yhat.yhat
@test yhat.continuousA[end] .* m_tul.beta[end-1] + m_tul.beta[20] ≈ yhat.yhat[end-3]
