data, evts = loadtestdata("test_case_3a")
##
f_spl = @formula 0 ~ 1 + conditionA + spl(continuousA, 4)
f = @formula 0 ~ 1 + conditionA + continuousA
data_r = reshape(data, (1, :))
data_e, times = Unfold.epoch(data = data_r, tbl = evts, τ = (-1.0, 1.0), sfreq = 10)

@testset "spline basis" begin
    m_mul_spl = coeftable(fit(UnfoldModel, f_spl, evts, data_e, times))

    # asking for 4 splines should generate 4 splines
    @test length(unique(m_mul_spl.coefname)) == 5

    s = formulas(fit(UnfoldModel, f_spl, evts, data_e, times))[1].rhs.terms[3]
    @test Unfold.width(s) == 3
    @test length(coefnames(s)) == 3
    @test s.df == 4
end

@testset "outside bounds" begin
    # test safe prediction
    m = fit(UnfoldModel, f_spl, evts, data_e, times)
    r = predicttable(m, DataFrame(conditionA = [0, 0], continuousA = [0.9, 1.9]))
    @test all(ismissing.(r.yhat[r.continuousA .== 1.9]))
    @test !any(ismissing.(r.yhat[r.continuousA .== 0.9]))
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

    p = predicttable(m, DataFrame(conditionA = [0, 0, 0], continuousA = [0.9, 0.9, 1.9]))
    @test all(ismissing.(p[p.continuousA .== 1.9, :yhat]))
end

@testset "many splines" begin
    # test much higher number of splines
    f_spl_many = @formula 0 ~ 1 + spl(continuousA, 131)
    m_mul_spl_many = coeftable(fit(UnfoldModel, f_spl_many, evts, data_e, times))
    @test length(unique(m_mul_spl_many.coefname)) == 131
end

@testset "PeriodicSplines" begin
    f_circspl = @formula 0 ~ 1 + circspl(continuousA, 10, -1, 1)
    m = fit(UnfoldModel, f_circspl, evts, data_e, times)

    effValues = range(-1.1, 1.1, step = 0.1)
    effSingle = effects(Dict(:continuousA => effValues), m)
    tmp = subset(effSingle, :time => x -> x .== -1.0)
    @test tmp.yhat[tmp.continuousA .== -1.1] ≈ tmp.yhat[tmp.continuousA .== 0.9]
    @test tmp.yhat[tmp.continuousA .== -1.0] ≈ tmp.yhat[tmp.continuousA .== 1]
    @test tmp.yhat[tmp.continuousA .== -0.9] ≈ tmp.yhat[tmp.continuousA .== 1.1]
end

@testset "minimal number of splines" begin
    f_spl = @formula 0 ~ 1 + conditionA + spl(continuousA, 3)
    @test_throws AssertionError fit(UnfoldModel, f_spl, evts, data_e, times)
end

@testset "SplineBasis shift_onset" begin
    bf = Unfold.splinebasis(τ = (-0.2, 1.0), sfreq = 100, nsplines = 10, name = "test")
    @test Unfold.shift_onset(bf) == Int64(floor(-0.2 * 100)) # should be -20
end
