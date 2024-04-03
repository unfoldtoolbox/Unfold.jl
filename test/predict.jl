
data, evts = loadtestdata("test_case_3a") #
data_r = reshape(data, (1, :))

data_e, times = Unfold.epoch(data = data_r, tbl = evts, τ = (0.0, 0.95), sfreq = 20) # cut the data into epochs
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

@test yhat_mul[1][1, :, :] ≈ yhat_tul[1][1, 1:20, :]

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
b1 = firbasis(τ = (0.0, 0.95), sfreq = 20)
b2 = firbasis(τ = (0.0, 0.95), sfreq = 20)
f = @formula 0 ~ 1 # 1
m_tul = fit(
    UnfoldModel,
    Dict("eventA" => (f, b1), "eventB" => (f, b2)),
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
    Dict("eventA" => (f, times), "eventB" => (f, times)),
    evts,
    data_e,
    eventcolumn = "type",
)

p = predict(m_mul, DataFrame(:Cond => [1, 2, 3]))
@test length(p) == 2
@test size(p[2]) == (1, 20, 3)

## result_to_table
data, evts = UnfoldSim.predef_eeg(; n_repeats = 5, noiselevel = 0.8)
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

@test all(pt[[1, 2, 3], :yhat] .== 0.293292)
@test all(pt[[1, 2, 3], :channel] .== [1, 2, 3])
@test all(pt[[1, 2, 3], :channel] .== [1, 2, 3])
@test all(
    pt[[1, 6 * 112 + 1, 3 * 112 + 1], :continuous] .≈ [5, 1.6666666667, -2.7777777778],
)
