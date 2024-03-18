
data, evts = loadtestdata("test_case_3a") #
data_r = reshape(data, (1, :))

data_e, times = Unfold.epoch(data = data_r, tbl = evts, τ = (0.0, 0.95), sfreq = 20) # cut the data into epochs
basisfunction = firbasis(τ = (0.0, 0.95), sfreq = 20, name = "basisA")

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


yhat_tul = predict(m_tul, evts_grid)
yhat_mul = predict(m_mul, evts_grid)

@test unique(yhat_mul.time) == times

@test size(m_tul_results)[1] * size(evts_grid)[1] == size(yhat_tul)[1]
@test size(m_mul_results)[1] * size(evts_grid)[1] == size(yhat_mul)[1]


@test yhat_mul[end-3, :yhat] ≈ mean(data[data.!=0])
@test yhat_tul[end-3, :yhat] ≈ mean(data[data.!=0])


@test yhat_mul.time[1] == 0.0
@test length(unique(yhat_mul.time)) == 20

@test length(unique(yhat_mul.yhat)) == 1

## Case with multiple formulas

f = @formula 0 ~ 1 + conditionA + continuousA# 1
m_tul = fit(UnfoldModel, f, evts, data_r, basisfunction)
m_mul = fit(UnfoldModel, f, evts, data_e, times)

m_mul_results = coeftable(m_mul)
m_tul_results = coeftable(m_tul)

yhat_tul = predict(m_tul, evts_grid)
yhat_mul = predict(m_mul, evts_grid)
#@test predict(m_mul,evts).yhat ≈ yhat.yhat


@test isapprox(sum(coef(m_tul)[[5, 25, 45]]), yhat_tul.yhat[end-3], atol = 0.001)


@test isapprox(
    yhat_tul[yhat_tul.time.==0.5, :yhat],
    yhat_mul[yhat_mul.time.==0.5, :yhat],
    atol = 0.001,
)
@test isapprox(yhat_tul[yhat_tul.time.==0.5, :yhat], [-2, 1, 2, 5, 6, 9.0], atol = 0.0001)

## two events

data, evts = loadtestdata("test_case_4a") #
b1 = firbasis(τ = (0.0, 0.95), sfreq = 20, name = "basisA")
b2 = firbasis(τ = (0.0, 0.95), sfreq = 20, name = "basisB")
f = @formula 0 ~ 1 # 1
m_tul = fit(
    UnfoldModel,
    Dict("eventA" => (f, b1), "eventB" => (f, b2)),
    evts,
    data,
    eventcolumn = "type",
)

p = predict(m_tul, DataFrame(:Cond => [1]))

@test size(p, 1) == 40
@test length(unique(p.time)) == 20
@test unique(p.eventname) == ["basisA", "basisB"]




##
if 1 == 0
    #for visualization
    using AlgebraOfGraphics
    yhat_mul.conditionA = categorical(yhat_mul.conditionA)
    m = mapping(:time, :yhat, color = :continuousA, linestyle = :conditionA)
    df = yhat_mul
    AlgebraOfGraphics.data(df) * visual(Lines) * m |> draw
end
