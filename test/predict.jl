using Test, StatsModels
using DataFrames
using StatsBase

using unfold
include("test_utilities.jl")


data,evts = loadtestdata("test_case_3a") #
data_r = reshape(data,(1,:))

data_e,times = unfold.epoch(data=data_r,tbl=evts,τ=(-1.,1),sfreq=20) # cut the data into epochs
basisfunction = firbasis(τ=(-1,1),sfreq=20,name="basisA")

# f  = @formula 0~ 1 * beta_0 + contA*beta_1
# fit => beta_0 = -1.5, beta_1 = 2
# predict => (contA=0.5) y_hat = 1 * -1.5 + 0.5 * 2 = -0.5

f  = @formula 0~1 # 1
m_mul,m_mul_results = fit(UnfoldLinearModel,f,evts,data_e,times)
m_tul,m_tul_results = fit(UnfoldLinearModel,f,evts,data_r,basisfunction)

conditionA  = [0,1.]
continuousA = [-1.,0,1.]

tmp = reshape([ [x,y]  for x=conditionA, y=continuousA ],length(conditionA)*length(continuousA))
evts_grid = DataFrame(hcat(tmp...)')
rename!(evts_grid,["conditionA","continuousA"])

yhat_tul = predict(m_tul,evts_grid)
yhat_mul = predict(m_mul,evts_grid)

# there seems to be a bug in yhat, that increases the times + 1 sample at the end. Investigate!
@test_broken size(m_tul_results)[1] * size(evts_grid)[1] == size(yhat_tul)[1]
@test size(m_mul_results)[1] * size(evts_grid)[1] == size(yhat_mul)[1]


@test yhat_mul[end-3,:yhat] ≈ mean(data[data.!=0])
@test_broken yhat_tul[end-3,:yhat] ≈ mean(data[data.!=0])


@test_broken yhat_mul.times[0] == -1.
@test_broken yhat_tul.yhat ≈ yhat_mul.yhat


f  = @formula 0~1+conditionA+continuousA# 1
m_tul,m_tul_results = fit(UnfoldLinearModel,f,evts,data_r,basisfunction)
m_mul,m_mul_results = fit(UnfoldLinearModel,f,evts,data_e,times)
yhat_tul = predict(m_tul,evts_grid)
yhat_mul = predict(m_mul,evts_grid)
@test predict(m_mul,evts).yhat ≈ yhat.yhat
@test yhat.continuousA[end] .* m_tul.beta[end-1] + m_tul.beta[20] ≈ yhat.yhat[end-3]


@test isapprox(yhat_tul[yhat_tul.times.==0.5,:yhat], yhat_mul[yhat_mul.times.==30,:yhat],atol=0.001) # fixme: Time fix in mul_results
@test isapprox(yhat_tul[yhat_tul.times.==0.5,:yhat], [-2,1,2,5,6,9.],atol=0.0001)

if 1 == 0
    #for visualization
    using AlgebraOfGraphics
    yhat_mul.conditionA = categorical(yhat_mul.conditionA)
    m = mapping(:times,:yhat,color=:continuousA,linestyle=:conditionA)
    df = yhat_mul
    AlgebraOfGraphics.data(df) * visual(Lines) * m  |> draw
end