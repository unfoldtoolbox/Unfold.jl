# Tests for UnfoldLinearModelContinuousTimeDenseDesign: a deconvolution model identical to
# UnfoldLinearModelContinuousTime, but storing the time-expanded designmatrix as a dense
# Matrix instead of a SparseMatrixCSC. Same user-facing interface, same default solver
# (solver_default). Cross-validation & regularized solvers are not implemented.

data, evts = loadtestdata("test_case_3a")
f = @formula 0 ~ 1 + conditionA + continuousA
data_r = reshape(data, (1, :))
data_r = vcat(data_r, deepcopy(data_r)) # add second channel

basisfunction = firbasis(τ = (-1, 1), sfreq = 20)

@testset "dense designmatrix construction" begin
    uf_s = fit(UnfoldLinearModelContinuousTime, f, evts, data_r, basisfunction)
    uf_d = fit(UnfoldLinearModelContinuousTimeDenseDesign, f, evts, data_r, basisfunction)

    @test typeof(uf_d) == UnfoldLinearModelContinuousTimeDenseDesign{Float64}
    @test typeof(designmatrix(uf_d)[1]) ==
          Unfold.DesignMatrixLinearModelContinuousTimeDense{Float64}
    @test modelmatrix(uf_d) isa Matrix
    # existing model unchanged
    @test typeof(uf_s) == UnfoldLinearModelContinuousTime{Float64}
    @test modelmatrix(uf_s) isa SparseMatrixCSC
    # identical numerics of the designmatrix
    @test modelmatrix(uf_d) == Matrix(modelmatrix(uf_s))
end

@testset "dense == sparse estimates (same solver)" begin
    uf_s = fit(UnfoldLinearModelContinuousTime, f, evts, data_r, basisfunction)
    uf_d = fit(UnfoldLinearModelContinuousTimeDenseDesign, f, evts, data_r, basisfunction)

    β_s = coef(uf_s)
    β_d = coef(uf_d)
    @test size(β_d) == size(β_s)
    # identical data + identical solver (lsmr), only the matrix storage differs
    @test β_d ≈ β_s rtol = 1e-5 atol = 1e-8

    # multi-event design with different basisfunction lengths: exercises the dense
    # extend_to_larger (vcat-0-rows + hcat) path
    evts_local = deepcopy(evts)
    evts_local.type .= repeat(["A", "B"], nrow(evts) ÷ 2)
    design = [
        "A" => (f, firbasis(τ = (-1, 1), sfreq = 20)),
        "B" => (f, firbasis(τ = (0, 0.5), sfreq = 20)),
    ]
    uf_s2 =
        fit(UnfoldLinearModelContinuousTime, design, evts_local, data_r; eventcolumn = :type)
    uf_d2 = fit(
        UnfoldLinearModelContinuousTimeDenseDesign,
        design,
        evts_local,
        data_r;
        eventcolumn = :type,
    )
    @test size(coef(uf_d2)) == size(coef(uf_s2))
    @test coef(uf_d2) ≈ coef(uf_s2) rtol = 1e-5 atol = 1e-8
end

@testset "dense model user-facing interface" begin
    uf_s = fit(UnfoldLinearModelContinuousTime, f, evts, data_r, basisfunction)
    uf_d = fit(UnfoldLinearModelContinuousTimeDenseDesign, f, evts, data_r, basisfunction)

    # coeftable
    tbl = coeftable(uf_d)
    @test nrow(tbl) == length(coef(uf_d))
    @test tbl.estimate ≈ coeftable(uf_s).estimate

    # predict (full overlap: X*b)
    @test predict(uf_d) ≈ predict(uf_s)

    # residuals & r2
    @test Unfold.r2(uf_d, data_r) ≈ Unfold.r2(uf_s, data_r)

    # effects: tolerance-based, estimates differ at solver-FP level between sparse/dense storage
    d = Dict(:conditionA => [0, 1])
    eff_d = effects(d, uf_d)
    eff_s = effects(d, uf_s)
    @test nrow(eff_d) == nrow(eff_s) && names(eff_d) == names(eff_s)
    @test eff_d[!, Not(:yhat)] == eff_s[!, Not(:yhat)]
    @test all(map(
        (a, b) -> isequal(a, b) || isapprox(a, b; rtol = 1e-5, atol = 1e-8),
        eff_d.yhat,
        eff_s.yhat,
    ))

    # modelmatrix without timeexpansion (widened from the sparse concrete type)
    @test modelmatrix(uf_d, false) == modelmatrix(uf_s, false)

    # times (Unfold.times: `times` is shadowed by a range in the shared test scope)
    @test Unfold.times(uf_d) == Unfold.times(uf_s)

    # manual pathway: construct + designmatrix! + fit! with default solver
    uf_manual = UnfoldLinearModelContinuousTimeDenseDesign{Float64}([Any => (f, basisfunction)])
    designmatrix!(uf_manual, evts; eventcolumn = "type")
    fit!(uf_manual, data_r)
    @test modelmatrix(uf_manual) isa Matrix
    @test !isempty(coef(uf_manual))
    @test coef(uf_manual) ≈ coef(uf_d) rtol = 1e-5 atol = 1e-8
end

@testset "dense end-to-end on known simulated effects" begin
    data, evts = UnfoldSim.predef_eeg(StableRNG(1); n_repeats = 1, noiselevel = 1)
    f = @formula 0 ~ 1 + condition
    m_d = fit(
        UnfoldLinearModelContinuousTimeDenseDesign,
        [Any => (f, firbasis((-0.1, 1), 100))],
        evts,
        data;
    )
    # known value measured with the identical sparse model (test/predict.jl "r2")
    _r2 = Unfold.r2(m_d, data)
    @test length(_r2) == 1
    @test isapprox(_r2[1], 0.74, atol = 0.01)

    # the simulated condition effect (face > car) is recovered: with the default
    # EffectsCoding there is a single "condition: face" coefficient = face - mean(face,car)
    tbl = coeftable(m_d)
    face = tbl[tbl.coefname.=="condition: face", :estimate]
    @test maximum(face) > 0
end

@testset "dense save/load roundtrip" begin
    uf_d = fit(UnfoldLinearModelContinuousTimeDenseDesign, f, evts, data_r, basisfunction)
    file = joinpath(tempdir(), "densedesign_roundtrip.jld2")
    try
        save(file, uf_d)
        uf_loaded = load(file, UnfoldModel, generate_Xs = true)
        @test typeof(uf_loaded) == typeof(uf_d)
        @test coef(uf_loaded) ≈ coef(uf_d)
        @test modelmatrix(uf_loaded) == modelmatrix(uf_d)
    finally
        rm(file, force = true)
    end
end
