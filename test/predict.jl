
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

@test all(yhat_mul[1][:] .≈ mean(data[data.!=0]))
@test all(yhat_tul[1][:] .≈ mean(data[data.!=0]))

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
pt = Unfold.result_to_table(
    m,
    p,
    [
        subset(evts, :condition => ByRow(==("car"))),
        subset(evts, :condition => ByRow(==("face"))),
    ],
) #repeat([evts], 2))

@show pt[[1, 2, 3], :yhat]
#@test_broken all(isapprox.(pt[[1, 2, 3], :yhat], 0.24672; atol = 0.01)) # test broken until UnfoldSim.jl is updated!!
@test all(pt[[1, 2, 3], :channel] .== [1, 2, 3])
# spot check to see if the order changed somehow
@test all(
    pt[[1, 5000, 25123], :yhat] .≈
    [0.23833130331025282, 0.07879460692911115, 0.016934637133599384],
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
    @test maximum(abs.(data_e .- (resids_e.+predict(m_mul)[1])[1, :, :])) < 0.0000001


    ##


    @test all(Unfold._residuals(UnfoldModel, [1 2 3; 3 4 5], [1 2 3; 3 4 5]) .== 0)

    #  y longer
    res = Unfold._residuals(UnfoldModel, [1 2 3; 3 4 5], [1 2 3 4; 3 4 5 6])
    @test all(res[:, 1:3] .== 0)
    @test res[:, 4] == [4, 6]

    # yhat longer
    @test all(Unfold._residuals(UnfoldModel, [1 2 3 4; 3 4 5 6], [1 2 3; 3 4 5]) .== 0)

end


@testset "non_zero_rows" begin
    dense = [0 1 0 3; 2 0 3 4; 0 0 0 0]
    @test Unfold.non_zero_rows(dense) == [true, true, false]
    @test Unfold.non_zero_rows(zeros(2, 3)) == [false, false]

    sparse_dense = sparse(dense)
    @test sort(Unfold.non_zero_rows(sparse_dense)) == [1, 2]
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
    m_e2 = fit(
        UnfoldModel,
        [
            "face" => (@formula(0 ~ 1), 1:size(data_e, 1)),
            "car" => (@formula(0 ~ 1), 1:size(data_e, 1)),
        ],
        evts,
        data_e;
        eventcolumn = :condition,
    )


    _r2 = Unfold.r2(m, data)
    @test length(_r2) == 1
    @test isapprox(_r2[1], 0.74, atol = 0.01)
    _r2 = Unfold.r2(m_e, data_e)
    _r2_e2 = Unfold.r2(m_e2, data_e)
    @test all(_r2 .≈ _r2_e2)
    @test length(_r2) == size(data_e, 1)
    @test all(_r2 .< 1)
    @test isapprox(_r2[1], 0.001, atol = 0.01)
    @test isapprox(_r2[16], 0.806, atol = 0.01)

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


    # Create data with sparse events (some rows with no modelling)
    data_sparse = allowmissing(data)
    data_sparse[1:2] .= missing
    evts_sparse = deepcopy(evts)
    evts_e_sparse = deepcopy(evts)

    data_e_sparse = allowmissing(data_e)
    data_e_sparse[1, :] .= missing
    data_e_sparse[:, 1] .= missing
    data_e_sparse[2, 2] = missing

    # Create models with sparse conditions to have non-modelled rows
    m_sparse = fit(
        UnfoldModel,
        [Any => (@formula(0 ~ 1 + condition), firbasis((-0.1, 1), 100))],
        evts_sparse,
        data_sparse;
    )
    m_e_sparse = fit(
        UnfoldModel,
        [Any => (@formula(0 ~ 1 + condition), 1:size(data_e_sparse, 1))],
        evts_sparse,
        data_e_sparse;
    )

    # check skipmissing false
    _r2 = Unfold.r2(m_sparse, data_sparse; skipmissing = false)
    @test ismissing(_r2[1])
    _r2 = Unfold.r2(m_e_sparse, data_e_sparse; skipmissing = false)
    @test all(ismissing.(_r2))


    # skipping 0 entries should increase r2
    # let's put it to the exterme, let's add a huge outlier
    data_sparse[3] = 1000
    _r2_skip = Unfold.r2(m_sparse, data_sparse; skip_notmodelled = true)
    _r2_notskip = Unfold.r2(m_sparse, data_sparse; skip_notmodelled = false)
    @test _r2_skip > _r2_notskip
    @test _r2_skip[1] > 0.75
    @test _r2_notskip[1] < 0.01


    # Test skip_notmodelled = true for epoched
    data_e_sparse[3, 3] = 1000
    _r2_skip = Unfold.r2(m_e_sparse, data_e_sparse; skip_notmodelled = true)
    _r2_notskip = Unfold.r2(m_e_sparse, data_e_sparse; skip_notmodelled = false)

    # Test skip_notmodelled with skipmissing=false
    _r2_no_skip_miss =
        Unfold.r2(m_sparse, data_sparse; skip_notmodelled = true, skipmissing = false)
    @test length(_r2_no_skip_miss) == 1

    # Test consistency between skip_notmodelled=false and default
    _r2_all = Unfold.r2(m_sparse, data_sparse; skip_notmodelled = false)
    _r2_default = Unfold.r2(m_sparse, data_sparse)
    @test isapprox(_r2_all[1], _r2_default[1])

end
