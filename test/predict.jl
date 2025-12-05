
data, evts = loadtestdata("test_case_3a") #
data_r = reshape(data, (1, :))

data_e, times = Unfold.epoch(data = data_r, evts = evts, τ = (0.0, 0.95), sfreq = 20) # cut the data into epochs
basisfunction = firbasis(τ = (0.0, 0.95), sfreq = 20)

f = @formula 0 ~ 1 # 1
m_mul = fit(UnfoldModel, f, evts, data_e, times)
m_tul = fit(UnfoldModel, f, evts, data_r, basisfunction)
m_mul_results = coeftable(m_mul)
m_tul_results = coeftable(m_tul)

conditionA = [0, 1.0]
continuousA = [-1.0, 0, 1.0]

tmp = reshape(
    [[x, y] for x in conditionA, y in continuousA],
    length(conditionA) * length(continuousA),
)
evts_grid = DataFrame(collect(hcat(tmp...)'), ["conditionA", "continuousA"])


yhat_mul = predict(m_mul, evts_grid)
yhat_tul = predict(m_tul, evts_grid)

@test yhat_mul ≈ yhat_tul

Unfold.result_to_table(m_mul, yhat_mul, [evts_grid])

@test all(yhat_mul[1][:] .≈ mean(data[data .!= 0]))
@test all(yhat_tul[1][:] .≈ mean(data[data .!= 0]))

## Case with multiple formulas
f = @formula 0 ~ 1 + conditionA + continuousA# 1
m_tul = fit(UnfoldModel, f, evts, data_r, basisfunction)
m_mul = fit(UnfoldModel, f, evts, data_e, times)

m_mul_results = coeftable(m_mul)
m_tul_results = coeftable(m_tul)

yhat_tul = predict(m_tul, evts_grid)
yhat_mul = predict(m_mul, evts_grid)
#@test predict(m_mul,evts).yhat ≈ yhat.yhat


@test isapprox(sum(coef(m_tul)[[5, 25, 45]]), yhat_tul[1][end-3], atol = 0.001)


@test isapprox(yhat_tul[1], yhat_mul[1], atol = 0.001)
@test isapprox(yhat_tul[1][1, 5, :], [-2, 1, 2, 5, 6, 9.0], atol = 0.0001)

## two events
data, evts = loadtestdata("test_case_4a") #
data_e, times = Unfold.epoch(data = data, evts = evts, τ = (0.0, 0.95), sfreq = 20) # cut the data into epochs
b1 = firbasis(τ = (0.0, 0.95), sfreq = 20)
b2 = firbasis(τ = (0.0, 0.95), sfreq = 20)
f = @formula 0 ~ 1 # 1
m_tul = fit(
    UnfoldModel,
    ["eventA" => (f, b1), "eventB" => (f, b2)],
    evts,
    data,
    eventcolumn = "type",
)

p = predict(m_tul, DataFrame(:Cond => [1]))
@test length(p) == 2
@test size(p[1], 2) == size(p[2], 2) == 20

## two events mass-univariate
m_mul = fit(
    UnfoldModel,
    ["eventA" => (f, times), "eventB" => (f, times)],
    evts,
    data_e,
    eventcolumn = "type",
)

p = predict(m_mul, DataFrame(:Cond => [1, 2, 3]))
@test length(p) == 2
@test size(p[2]) == (1, 20, 3)

## result_to_table
data, evts = UnfoldSim.predef_eeg(StableRNG(1); n_repeats = 5, noiselevel = 0.8)
m = fit(
    UnfoldModel,
    [
        "car" => (@formula(0 ~ 1 + spl(continuous, 4)), firbasis((-0.1, 1), 100)),
        "face" => (@formula(0 ~ 1 + continuous^2), firbasis((-0.2, 0.4), 100)),
    ],
    evts,
    repeat(data, 1, 3)';
    eventcolumn = "condition",
    show_progress = false,
)
p = predict(m; overlap = false)
pt = Unfold.result_to_table(m, p, repeat([evts], 2))

@show pt[[1, 2, 3], :yhat]
#@test_broken all(isapprox.(pt[[1, 2, 3], :yhat], 0.24672; atol = 0.01)) # test broken until UnfoldSim.jl is updated!!
@test all(pt[[1, 2, 3], :channel] .== [1, 2, 3])
@test all(
    pt[[1, 6 * 112 + 1, 3 * 112 + 1], :continuous] .≈
    [-0.5555555555555556, 2.7777777777777777, 3.888888888888889],
)



@testset "residuals" begin
    data, evts = UnfoldSim.predef_eeg(StableRNG(1); n_repeats = 5, noiselevel = 0.8)

    # time expanded
    m = fit(UnfoldModel, [Any => (@formula(0 ~ 1), firbasis((-0.1, 1), 100))], evts, data;)
    @test size(Unfold.residuals(m, data)) == (1, length(data))

    # time expanded + multichannel
    m = fit(
        UnfoldModel,
        [Any => (@formula(0 ~ 1), firbasis((-0.1, 1), 100))],
        evts,
        repeat(data, 1, 3)';
    )
    @test size(Unfold.residuals(m, data)) == (3, length(data))


    # time expanded, data longer
    m = fit(
        UnfoldModel,
        [Any => (@formula(0 ~ 1), firbasis((-0.1, 0.1), 100))],
        evts,
        repeat(data, 1, 3)';
    )
    resids = Unfold.residuals(m, repeat(data, 1, 3)')
    @test size(resids) == (3, length(data))
    @test all(resids[1, (end-2):end] .≈ data[(end-2):end])

    #

    data_e, evts =
        UnfoldSim.predef_eeg(; n_repeats = 5, noiselevel = 0.8, return_epoched = true)



    times = 1:size(data_e, 1)
    m_mul = fit(UnfoldModel, @formula(0 ~ 1), evts, data_e, times)
    resids_e = Unfold.residuals(m_mul, data_e)

    @test size(resids_e)[2:3] == size(data_e)
    @test maximum(abs.(data_e .- (resids_e .+ predict(m_mul)[1])[1, :, :])) < 0.0000001


    ##


    @test all(Unfold._residuals(UnfoldModel, [1 2 3; 3 4 5], [1 2 3; 3 4 5]) .== 0)

    #  y longer
    res = Unfold._residuals(UnfoldModel, [1 2 3; 3 4 5], [1 2 3 4; 3 4 5 6])
    @test all(res[:, 1:3] .== 0)
    @test res[:, 4] == [4, 6]

    # yhat longer
    @test all(Unfold._residuals(UnfoldModel, [1 2 3 4; 3 4 5 6], [1 2 3; 3 4 5]) .== 0)

end


@testset "r2" begin
    data, evts = UnfoldSim.predef_eeg(StableRNG(1); n_repeats = 1, noiselevel = 1)
    data_e, _ = UnfoldSim.predef_eeg(
        StableRNG(1);
        n_repeats = 1,
        noiselevel = 1,
        return_epoched = true,
    )

    # time expanded
    m = fit(
        UnfoldModel,
        [Any => (@formula(0 ~ 1 + condition), firbasis((-0.1, 1), 100))],
        evts,
        data;
    )
    m_e = fit(
        UnfoldModel,
        [Any => (@formula(0 ~ 1 + condition), 1:size(data_e, 1))],
        evts,
        data_e;
    )
    _r2 = Unfold.r2(m, data)
    @test length(_r2) == 1
    @test isapprox(_r2[1], 0.74, atol = 0.01)
    _r2 = Unfold.r2(m_e, data_e)
    @test length(_r2) == size(data_e, 1)
    @test all(_r2 .< 1)
    @test isapprox(_r2[1], 0.001, atol = 0.01)
    @test isapprox(_r2[16], 0.82, atol = 0.01)

    data_reshape = reshape(data, 1, :)
    data_e_reshape = reshape(data_e, 1, size(data_e)...)
    m_reshape = fit(
        UnfoldModel,
        [Any => (@formula(0 ~ 1 + condition), firbasis((-0.1, 1), 100))],
        evts,
        data_reshape;
    )
    m_e_reshape = fit(
        UnfoldModel,
        [Any => (@formula(0 ~ 1 + condition), 1:size(data_e, 1))],
        evts,
        data_e_reshape;
    )
    _r2 = Unfold.r2(m_reshape, data_reshape)
    @test length(_r2) == 1
    _r2 = Unfold.r2(m_e_reshape, data_e_reshape)
    @test size(_r2) == (1, size(data_e_reshape, 2))

end
