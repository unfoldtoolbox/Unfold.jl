data, evts = loadtestdata("test_case_3a") #
##

f_spl = @formula 0 ~ 1 + conditionA + spl(continuousA, 4) # 1
f = @formula 0 ~ 1 + conditionA + continuousA # 1
data_r = reshape(data, (1, :))
data_e, times = Unfold.epoch(data = data_r, tbl = evts, τ = (-1.0, 1.0), sfreq = 10) # cut the data into epochs

m_mul = coeftable(fit(UnfoldModel, f, evts, data_e, times))
m_mul_spl = coeftable(fit(UnfoldModel, f_spl, evts, data_e, times))

# asking for 4 splines should generate 4 splines 
@test length(unique(m_mul_spl.coefname)) == 5

s = Unfold.formulas(fit(UnfoldModel, f_spl, evts, data_e, times))[1].rhs.terms[3]
@test Unfold.width(s) == 3
@test length(coefnames(s)) == 3
@test s.df == 4

@testset "outside bounds" begin
    # test safe prediction
    m = fit(UnfoldModel, f_spl, evts, data_e, times)
    r = predict(m, DataFrame(conditionA = [0, 0], continuousA = [0.9, 1.9]))
    @test all(ismissing.(r.yhat[r.continuousA.==1.9]))
    @test !any(ismissing.(r.yhat[r.continuousA.==0.9]))
end

basisfunction = firbasis(τ = (-1, 1), sfreq = 10, name = "A")
@testset "timeexpanded" begin
    # test time expanded
    m_tul = coeftable(fit(UnfoldModel, f, evts, data_r, basisfunction))
    m_tul_spl = coeftable(fit(UnfoldModel, f_spl, evts, data_r, basisfunction))
end
@testset "safe prediction outside bounds" begin
    # test safe predict
    m = fit(UnfoldModel, f_spl, evts, data_r, basisfunction)

    p = predict(m, DataFrame(conditionA = [0, 0, 0], continuousA = [0.9, 0.9, 1.9]))
    @test all(ismissing.(p[p.continuousA.==1.9, :yhat]))

end
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
@testset "many splines" begin
    # test much higher number of splines
    f_spl_many = @formula 0 ~ 1 + spl(continuousA, 131) # 1
    m_mul_spl_many = coeftable(fit(UnfoldModel, f_spl_many, evts, data_e, times))
    @test length(unique(m_mul_spl_many.coefname)) == 131
end

@testset "PeriodicSplines" begin
    f_circspl = @formula 0 ~ 1 + circspl(continuousA, 10, -1, 1) # 1
    m = fit(UnfoldModel, f_circspl, evts, data_e, times)
    f_evaluated = Unfold.formula(m)

    effValues = [-1, -0.99, 0, 0.99, 1]
    effValues = range(-1.1, 1.1, step = 0.1)
    effSingle = effects(Dict(:continuousA => effValues), m)
    tmp = subset(effSingle, :time => x -> x .== -1.0)
    @test tmp.yhat[tmp.continuousA.==-1.1] ≈ tmp.yhat[tmp.continuousA.==0.9]
    @test tmp.yhat[tmp.continuousA.==-1.0] ≈ tmp.yhat[tmp.continuousA.==1]
    @test tmp.yhat[tmp.continuousA.==-0.9] ≈ tmp.yhat[tmp.continuousA.==1.1]

end

@testset "minimal number of splines" begin
    f_spl = @formula 0 ~ 1 + conditionA + spl(continuousA, 3) # 1
    @test_throws AssertionError fit(UnfoldModel, f_spl, evts, data_e, times)

end
