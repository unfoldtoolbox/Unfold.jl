
include("test_utilities.jl")
include("setup.jl")
include("../src/circsplinepredictors.jl")

data, evts = loadtestdata("test_case_3a") #
#
##
f_circspl = @formula 0 ~ 1 + conditionA + circspl(continuousA, 4, 0.1, 10) # 1
f_spl = @formula 0 ~ 1 + conditionA + spl(continuousA, 4) # 1
f = @formula 0 ~ 1 + conditionA + continuousA # 1
data_r = reshape(data, (1, :))
data_e, times = Unfold.epoch(data = data_r, tbl = evts, τ = (-1.0, 1.0), sfreq = 10) # cut the data into epochs

m_mul = coeftable(fit(UnfoldModel, f, evts, data_e, times))
# to experiment with circular splines. TODO: remove them once test cases got implemented
#test_circspl = fit(UnfoldModel, f_circspl, evts, data_e, times)
#test_spl = fit(UnfoldModel, f_spl, evts, data_e, times)
m_mul_spl = coeftable(fit(UnfoldModel, f_spl, evts, data_e, times))
m_mul_circspl = coeftable(fit(UnfoldModel, f_circspl, evts, data_e, times))

# asking for 4 splines should generate 4 splines 
@test length(unique(m_mul_spl.coefname)) == 6 # XXX check back with Unfold whether this is the same! could be n-1 splines in Unfold. We should keep that comparable I guess

# test safe prediction
m = fit(UnfoldModel, f_spl, evts, data_e, times)
r = predict(m,DataFrame(conditionA=[0,0],continuousA=[0.9,1.9]))
@test all(ismissing.(r.yhat[r.continuousA.==1.9]))
@test !any(ismissing.(r.yhat[r.continuousA.==0.9]))


# test time expanded
basisfunction = firbasis(τ = (-1, 1), sfreq = 10, name = "A")
m_tul = coeftable(fit(UnfoldModel, f, evts, data_r, basisfunction))
m_tul_spl = coeftable(fit(UnfoldModel, f_spl, evts, data_r, basisfunction))

# test safe predict
m = fit(UnfoldModel, f_spl, evts, data_r, basisfunction)
@test_broken predict(m,DataFrame(conditionA=[0,0,0],continuousA=[0.9,0.9,1.9]))
#@test_broken all(ismissing.)

#evts_grid = gridexpand() 
# results from timeexpanded and non should be equal
#yhat_tul  = predict(m_tul_spl,evts_grid)
#yhat_mul  = predict(m_mul_spl,evts_grid)
if 1 == 0
    using AlgebraOfGraphics
    yhat_mul.conditionA = categorical(yhat_mul.conditionA)
    yhat_mul.continuousA = categorical(yhat_mul.continuousA)
    m = mapping(:times, :yhat, color = :continuousA, linestyle = :conditionA)
    df = yhat_mul
    AlgebraOfGraphics.data(df) * visual(Lines) * m |> draw
end

# test much higher number of splines
f_spl_many = @formula 0 ~ 1 + spl(continuousA, 131) # 1
m_mul_spl_many = coeftable(fit(UnfoldModel, f_spl_many, evts, data_e, times))
@test length(unique(m_mul_spl_many.coefname)) == 132


## some circular test_spl
m = fit(UnfoldModel, f_circspl, evts, data_e, times)
f_evaluated = Unfold.formula(m)
@test size(f_evaluated.rhs.term[3].fun([1])) == (1,5)
@test size(modelcols(form.rhs.terms[3],DataFrame(conditionA=[0],continuousA=[1]))) == (1,4)