##---
# VIF: variance inflation factor of the designmatrix, checks for collinearity

@testset "core matrix" begin
    # r = 0.5 -> VIF = 1/(1-r^2) = 4/3 for both columns
    @test all(Unfold._vif([1.0 2.0; 2.0 0.0; 3.0 4.0]) .≈ 4/3)
    # r^2 = 25/28 -> VIF = 28/3
    v = Unfold._vif([1.0 1.0 1.0; 1.0 2.0 2.0; 1.0 3.0 6.0])
    @test isinf(v[1]) # constant column (e.g. intercept) has no defined VIF
    @test all(v[2:3] .≈ 28/3)
    # orthogonal columns -> VIF = 1
    v = Unfold._vif([1.0 1.0 1.0; 1.0 2.0 -1.0; 1.0 3.0 -1.0; 1.0 4.0 1.0])
    @test isinf(v[1])
    @test all(v[2:3] .≈ 1.0)
    # perfectly collinear columns (col3 = col2): one of them must be Inf
    v = Unfold._vif([1.0 1.0 1.0; 1.0 2.0 2.0; 1.0 3.0 3.0; 1.0 4.0 4.0])
    @test isinf(v[1])
    @test count(isinf.(v[2:3])) == 1
    @test all(v[2:3][.!isinf.(v[2:3])] .≈ 1.0)
    # scaled duplicate (col3 = 2*col2) is detected as collinear as well
    v = Unfold._vif([1.0 1.0 2.0; 1.0 2.0 4.0; 1.0 3.0 6.0; 1.0 4.0 8.0])
    @test count(isinf.(v)) == 2
    # near-collinear -> large but finite (r^2 ~ 0.9931 -> VIF ~ 145)
    v = Unfold._vif([1.0 1.0 1.1; 1.0 2.0 1.9; 1.0 3.0 3.1; 1.0 4.0 3.9])
    @test all(isfinite.(v[2:3]))
    @test all(v[2:3] .≈ 145.0)
    # degenerate cases
    @test only(Unfold._vif(reshape([1.0, 2.0, 3.0], 3, 1))) ≈ 1.0
    @test all(isinf.(Unfold._vif(ones(4, 2))))
end

@testset "sparse == dense" begin
    onsets = [3, 10, 17, 24]
    w = 4
    rng = MersenneTwister(42)
    Xs = sparse(
        [o + t - 1 for o in onsets for t = 1:w for j = 1:2],
        [(j - 1) * w + t for o in onsets for t = 1:w for j = 1:2],
        rand(rng, length(onsets) * w * 2),
        30,
        2 * w,
    )
    @test norm(Unfold._vif(Xs) - Unfold._vif(Matrix(Xs)), Inf) < 1e-10
end

data, evts = loadtestdata("test_case_3a")
data_r = reshape(data, (1, :))
f = @formula 0 ~ 1 + conditionA + continuousA

@testset "time-expanded model" begin
    basisfunction = firbasis(τ = (-1.0, 0.9), sfreq = 20)
    m = fit(
        UnfoldModel,
        [Any => (f, basisfunction)],
        evts,
        data_r;
        fit = false,
        show_warnings = false,
    )
    X = modelmatrix(m)
    v = vif(m)
    @test size(v) == (size(X, 2), 2)
    @test v.coefname == Unfold.get_coefnames(m)
    @test all(v.VIF .≈ Unfold._vif(Matrix(X)))
    # cross-check against a direct (dense) correlation-matrix computation
    Xf = Matrix(X)
    Xc = Xf .- mean(Xf, dims = 1)
    ss = sqrt.(vec(sum(Xc .^ 2, dims = 1)))
    R = (Xc' * Xc) ./ (ss * ss')
    @test all(v.VIF .≈ diag(inv(R)))
end

@testset "mass-univariate model" begin
    data_e, times = Unfold.epoch(data = data_r, tbl = evts, τ = (-1.0, 1.9), sfreq = 20)
    # single event, unfitted: vif only needs the designmatrix
    m1 = fit(
        UnfoldModel,
        [Any => (f, times)],
        evts,
        data_e;
        fit = false,
        show_warnings = false,
    )
    v1 = vif(m1)
    @test size(v1) == (3, 3)
    @test v1.eventname == repeat(["Any"], 3)
    @test v1.coefname == ["(Intercept)", "conditionA", "continuousA"]
    @test isinf(v1.VIF[1])
    @test all(v1.VIF[2:3] .≈ Unfold._vif(modelmatrix(m1)[1])[2:3])

    # multiple events
    evts2 = deepcopy(evts)
    evts2.type .= repeat(["A", "B"], nrow(evts2) ÷ 2)
    m2 = fit(
        UnfoldModel,
        ["A" => (f, times), "B" => (f, times)],
        evts2,
        data_e;
        eventcolumn = "type",
        fit = false,
        show_warnings = false,
    )
    v2 = vif(m2)
    @test size(v2) == (6, 3)
    @test v2.eventname == repeat(["A", "B"], inner = 3)
    @test all(isinf.(v2.VIF[[1, 4]]))
    @test all(v2.VIF[2:3] .≈ Unfold._vif(modelmatrix(m2)[1])[2:3])
    @test all(v2.VIF[5:6] .≈ Unfold._vif(modelmatrix(m2)[2])[2:3])
end
