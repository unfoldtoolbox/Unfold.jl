using Logging

@testset "solver_cv" begin
    @test isa(solver_cv(; n_folds = 3), Function)

    X = reshape(collect(1.0:24.0), 12, 2)
    y = reshape(collect(1.0:12.0), 1, 1, 12)

    inner_solver =
        (X_, y_) -> begin
            train_n = size(X_, 1)
            train_sum = sum(y_)
            estimate = reshape([train_n, train_sum], 1, 1, 2)
            standarderror = reshape([train_n / 10, train_sum / 10], 1, 1, 2)
            return Unfold.LinearModelFit(
                estimate,
                (n = train_n, ysum = train_sum),
                standarderror,
            )
        end

    cv = solver_cv(;
        n_folds = 3,
        shuffle = false,
        inner_solver = inner_solver,
        fit_test = true,
    )
    fit_cv = cv(X, y)

    @test fit_cv isa LinearModelFitCV
    @test size(fit_cv.estimate) == (1, 1, 2, 3)
    @test size(fit_cv.standarderror) == (1, 1, 2, 3)
    @test length(fit_cv.folds) == 3
    @test length(fit_cv.info.inner_infos) == 3
    @test size(fit_cv.info.test_estimates) == (1, 1, 2, 3)

    # For 12 trials and 3 folds, each train split has 8 samples.
    @test all(fit_cv.estimate[1, 1, 1, :] .== 8)
    # With fit_test=true the held-out splits have 4 samples each.
    @test all(fit_cv.info.test_estimates[1, 1, 1, :] .== 4)

    all_test_indices = sort!(vcat([collect(f.test) for f in fit_cv.folds]...))
    @test all_test_indices == collect(1:12)

    # the per-fold test sets can be extracted with fold_assignment
    @test sort(vcat(fold_assignment(fit_cv)...)) == collect(1:12)
end

@testset "LinearModelFitCV" begin
    estimate = reshape(collect(1.0:24.0), 2, 3, 2, 2)
    standarderror = estimate ./ 10
    m = LinearModelFitCV(estimate, (;), standarderror, [])

    @test coef(m) == dropdims(mean(estimate, dims = 4), dims = 4)
    @test Unfold.stderror(m) == dropdims(mean(standarderror, dims = 4), dims = 4)
end

@testset "_cat" begin
    m1 = Unfold.LinearModelFit(
        reshape([1.0, 2.0], 1, 2, 1),
        "first",
        reshape([0.1, 0.2], 1, 2, 1),
    )
    m2 = Unfold.LinearModelFit(
        reshape([3.0, 4.0], 1, 2, 1),
        "second",
        reshape([0.3, 0.4], 1, 2, 1),
    )

    cat_lm = Unfold._cat([m1, m2])
    @test cat_lm isa Unfold.LinearModelFit
    @test size(cat_lm.estimate) == (1, 2, 2)
    @test size(cat_lm.standarderror) == (1, 2, 2)
    @test cat_lm.info == ["first", "second"]
    @test cat_lm.estimate[:, :, 1] == m1.estimate[:, :, 1]
    @test cat_lm.estimate[:, :, 2] == m2.estimate[:, :, 1]

    c1 = LinearModelFitCV(
        reshape([1.0, 2.0], 1, 1, 1, 2),
        (id = 1,),
        reshape([0.1, 0.2], 1, 1, 1, 2),
        [(train = [1, 2], test = [3])],
    )
    c2 = LinearModelFitCV(
        reshape([3.0, 4.0], 1, 1, 1, 2),
        (id = 2,),
        reshape([0.3, 0.4], 1, 1, 1, 2),
        [(train = [4, 5], test = [6])],
    )

    cat_cv = Unfold._cat([c1, c2])
    @test cat_cv isa LinearModelFitCV
    @test size(cat_cv.estimate) == (1, 1, 2, 2)
    @test size(cat_cv.standarderror) == (1, 1, 2, 2)
    @test cat_cv.info == [(id = 1,), (id = 2,)]
    @test cat_cv.folds == [c1.folds, c2.folds]
end

using Logging

@testset "solver_cv deconvolved" begin
    # a sparse design matrix with 3 active blocks separated by all-zero gaps
    n = 100
    X = spzeros(Float64, n, 2)
    X[1:30, 1] .= 1.0
    X[1:30, 2] .= 0.5
    X[40:70, 1] .= 1.0
    X[40:70, 2] .= -0.5
    X[80:100, 1] .= 1.0
    X[80:100, 2] .= 0.25
    y = rand(MersenneTwister(1), 2, n)

    # stub inner solver that records how many rows it was fitted on
    inner_solver = (X_, y_) -> begin
        n_rows = size(X_, 1)
        estimate = [n_rows n_rows; sum(y_) sum(y_)]
        standarderror = estimate ./ 10
        return Unfold.LinearModelFit(estimate, (n_rows = n_rows,), standarderror)
    end

    # --- auto-detected segments (maximal stretches of non-zero rows) ---
    m = solver_cv(n_folds = 2, shuffle = false, inner_solver = inner_solver)(X, y)

    @test m isa LinearModelFitCV{Float64, 3}
    @test size(m.estimate) == (2, 2, 2)
    @test size(m.standarderror) == (2, 2, 2)
    @test m.info.segments == UnitRange{Int}[1:30, 40:70, 80:100]
    @test m.info.test_estimates === nothing
    @test length(m.folds) == 2
    # every segment is held out exactly once
    @test sort(vcat([collect(f.test) for f in m.folds]...)) == [1, 2, 3]
    # no leakage: each fold is fitted exactly on the complement of its test segments
    for (f, info) in zip(m.folds, m.info.inner_infos)
        n_test_rows = sum(length(s) for s in m.info.segments[f.test])
        @test info.n_rows == n - n_test_rows
    end
    # fold-averaged coef / stderror over the (3rd) fold dimension
    @test coef(m) == dropdims(mean(m.estimate, dims = 3), dims = 3)
    @test Unfold.stderror(m) == dropdims(mean(m.standarderror, dims = 3), dims = 3)

    # reproducibility with a fixed seed
    m_a = solver_cv(
        n_folds = 2,
        shuffle = true,
        rng = MersenneTwister(3),
        inner_solver = inner_solver,
    )(X, y)
    m_b = solver_cv(
        n_folds = 2,
        shuffle = true,
        rng = MersenneTwister(3),
        inner_solver = inner_solver,
    )(X, y)
    @test m_a.estimate == m_b.estimate
    @test m_a.folds == m_b.folds

    # --- manual breaks in quiet intervals: no warning ---
    @test_logs min_level = Logging.Warn begin
        m_clean = solver_cv(
            n_folds = 3,
            shuffle = false,
            inner_solver = inner_solver,
            breaks = [40, 80],
        )(X, y)
        @test m_clean.info.segments == UnitRange{Int}[1:39, 40:79, 80:100]
    end

    # --- manual break inside an active response: warn about the boundary ---
    @test_logs (:warn, r"row\(s\) \[49\]") begin
        m_dirty = solver_cv(
            n_folds = 3,
            shuffle = false,
            inner_solver = inner_solver,
            breaks = [50, 80],
        )(X, y)
        @test m_dirty.info.segments == UnitRange{Int}[1:49, 50:79, 80:100]
    end

    # breaks can also be given as a DataFrame with a :latency column
    m_df = solver_cv(
        n_folds = 3,
        shuffle = false,
        inner_solver = inner_solver,
        breaks = DataFrame(latency = [40, 80]),
    )(X, y)
    @test m_df.info.segments == UnitRange{Int}[1:39, 40:79, 80:100]

    # n_folds cannot exceed the number of segments
    @test_throws ArgumentError solver_cv(n_folds = 5, shuffle = false, inner_solver = inner_solver)(
        X, y,
    )

    # fit_test=true also returns the estimates of the held-out segments
    m_ft = solver_cv(
        n_folds = 2,
        shuffle = false,
        inner_solver = inner_solver,
        fit_test = true,
    )(X, y)
    @test size(m_ft.info.test_estimates) == (2, 2, 2)
end

@testset "solver_cv deconvolved fold_assignment" begin
    n = 100
    # breaks at 20, 40, 60, 80 -> 5 segments (lengths 19, 20, 20, 20, 21);
    # the active rows are strictly inside the segments, so all boundaries
    # are clean (no leakage warnings)
    X = spzeros(Float64, n, 2)
    X[5:15, 1] .= 1.0
    X[24:35, 2] .= 1.0
    X[44:55, 1] .= 1.0
    X[64:75, 2] .= 1.0
    X[84:95, 1] .= 1.0
    y = rand(MersenneTwister(1), 2, n)

    inner_solver = (X_, y_) -> begin
        n_rows = size(X_, 1)
        estimate = [n_rows; sum(y_)]
        return Unfold.LinearModelFit(estimate, (n_rows = n_rows,), estimate ./ 10)
    end

    breaks = [20, 40, 60, 80]
    # with shuffle=false and n_folds = n_segments, fold i holds out segment i
    m = solver_cv(
        n_folds = 5, shuffle = false, inner_solver = inner_solver, breaks = breaks,
    )(X, y)
    @test fold_assignment(m) == [[1], [2], [3], [4], [5]]
    @test [info.n_rows for info in m.info.inner_infos] == [81, 80, 80, 80, 79]

    # reusing the same split (even with a different rng/shuffle) gives the
    # identical fit
    m2 = solver_cv(
        n_folds = 5,
        shuffle = true,
        rng = MersenneTwister(42),
        inner_solver = inner_solver,
        breaks = breaks,
        fold_assignment = fold_assignment(m),
    )(X, y)
    @test fold_assignment(m2) == fold_assignment(m)
    @test m2.estimate == m.estimate
    @test m2.folds == m.folds

    # a permuted assignment: fold 1 now holds out segment 2 (20 rows)
    fa = [[2], [1], [3], [4], [5]]
    m3 = solver_cv(
        n_folds = 5,
        shuffle = false,
        inner_solver = inner_solver,
        breaks = breaks,
        fold_assignment = fa,
    )(X, y)
    @test fold_assignment(m3) == fa
    @test m3.info.inner_infos[1].n_rows == 80
    @test m3.info.inner_infos[2].n_rows == 81

    # invalid assignments are rejected
    @test_throws ArgumentError solver_cv(
        n_folds = 5, inner_solver = inner_solver, breaks = breaks,
        fold_assignment = [[1], [2], [3]],
    )(X, y) # too few entries
    @test_throws ArgumentError solver_cv(
        n_folds = 5, inner_solver = inner_solver, breaks = breaks,
        fold_assignment = [[1], [2], [3], [5], [6]],
    )(X, y) # segment 4 missing, 6 out of range
    @test_throws ArgumentError solver_cv(
        n_folds = 5, inner_solver = inner_solver, breaks = breaks,
        fold_assignment = [[1, 1], [2], [3], [4], [5]],
    )(X, y) # duplicated segment
end

@testset "solver_cv deconvolved end-to-end" begin
    # box-shaped responses to 8 block onsets: the CV estimates should
    # predict the held-out blocks well
    sfreq = 100
    onsets = [50, 175, 300, 425, 550, 675, 800, 925]
    evts = DataFrame(latency = onsets, condition = repeat(["a", "b"], outer = 4))
    n = 1000
    data = 0.5 * randn(MersenneTwister(2), 2, n)
    τ = (-0.2, 0.6)
    for o in onsets
        w = (o - 20):(o + 60)
        w = w[(w .>= 1) .& (w .<= n)]
        data[:, w] .+= 8.0
    end
    bf = firbasis(τ = τ, sfreq = sfreq)
    f = @formula 0 ~ 1 + condition

    m0 = fit(UnfoldModel, [Any => (f, bf)], evts, data; fit = false)
    Xd = modelmatrix(m0)
    n_used = min(size(Xd, 1), n)

    m_cv = fit(
        UnfoldModel,
        [Any => (f, bf)],
        evts,
        data;
        solver = solver_cv(n_folds = 4, rng = MersenneTwister(1)),
    )
    mf = modelfit(m_cv)

    @test mf isa LinearModelFitCV
    @test size(mf.estimate) == (2, size(Xd, 2), 4)
    # coeftable works on the CV fit and reports the fold-averaged estimates
    @test size(coeftable(m_cv), 1) == 2 * size(Xd, 2)

    # out-of-sample R^2 against the grand mean should be high for all folds
    yall = data[:, 1:n_used]
    ȳ = mean(yall, dims = 2)
    for (f_, k) in zip(mf.folds, 1:4)
        test_rows = vcat(mf.info.segments[f_.test]...)
        test_rows = test_rows[test_rows .<= n_used]
        ytest = yall[:, test_rows]
        yhat = mf.estimate[:, :, k] * Xd[test_rows, :]'
        ss_res = vec(sum((ytest .- yhat).^2, dims = 2))
        ss_tot = vec(sum((ytest .- ȳ).^2, dims = 2))
        @test all((1 .- ss_res ./ ss_tot) .> 0.5)
    end
end

@testset "solver_cv deconvolved fold_assignment alignment" begin
    # delta-R² use case: a full model and a reduced model are evaluated on
    # the identical held-out segments by reusing the full model's split
    sfreq = 100
    onsets = [50, 175, 300, 425, 550, 675, 800, 925]
    evts = DataFrame(latency = onsets, condition = repeat(["a", "b"], outer = 4))
    n = 1000
    bf = firbasis(τ = (-0.2, 0.6), sfreq = sfreq)
    f_full = @formula 0 ~ 1 + condition
    f_red = @formula 0 ~ 1
    # breaks at the start of each response window: 9 segments with clean boundaries
    breaks = onsets .- 20

    # condition modulates the response amplitude, so the full model is clearly
    # better than the reduced (common-response-only) model
    data = 0.5 * randn(MersenneTwister(2), 2, n)
    for (o, c) in zip(onsets, evts[!, :condition])
        data[:, (o - 20):(o + 60)] .+= (c == "a" ? 6.0 : 10.0)
    end

    m_full = fit(
        UnfoldModel,
        [Any => (f_full, bf)],
        evts,
        data;
        solver = solver_cv(n_folds = 4, rng = MersenneTwister(1), breaks = breaks),
    )
    fa = fold_assignment(modelfit(m_full))
    @test length(fa) == 4
    @test sort(vcat(fa...)) == collect(1:9)

    # the reduced model's design has the same all-zero gaps (event-level terms
    # only), but the same manual breaks plus the extracted assignment keep the
    # split aligned with the full model's
    m_red = fit(
        UnfoldModel,
        [Any => (f_red, bf)],
        evts,
        data;
        solver = solver_cv(n_folds = 4, breaks = breaks, fold_assignment = fa),
    )
    mf_full = modelfit(m_full)
    mf_red = modelfit(m_red)

    @test fold_assignment(mf_red) == fa
    # identical held-out segments (and hence rows) per fold
    for (f_full, f_red) in zip(mf_full.folds, mf_red.folds)
        @test vcat(mf_full.info.segments[f_full.test]...) ==
              vcat(mf_red.info.segments[f_red.test]...)
    end

    # the full model predicts the held-out segments at least as well as the
    # reduced model, and on average clearly better (delta-R² > 0)
    Xd_full = modelmatrix(fit(UnfoldModel, [Any => (f_full, bf)], evts, data; fit = false))
    Xd_red = modelmatrix(fit(UnfoldModel, [Any => (f_red, bf)], evts, data; fit = false))
    n_used = min(size(Xd_full, 1), n)
    yall = data[:, 1:n_used]
    ȳ = mean(yall, dims = 2)
    function oos_r2(mf, Xd)
        r2s = Vector{Vector{Float64}}()
        for (f_, k) in zip(mf.folds, 1:4)
            test_rows = vcat(mf.info.segments[f_.test]...)
            test_rows = test_rows[test_rows .<= n_used]
            ytest = yall[:, test_rows]
            yhat = mf.estimate[:, :, k] * Xd[test_rows, :]'
            ss_res = vec(sum((ytest .- yhat) .^ 2, dims = 2))
            ss_tot = vec(sum((ytest .- ȳ) .^ 2, dims = 2))
            push!(r2s, 1 .- ss_res ./ ss_tot)
        end
        return hcat(r2s...)
    end
    r2_full = oos_r2(mf_full, Xd_full)
    r2_red = oos_r2(mf_red, Xd_red)
    @test all(r2_full .>= r2_red)
    @test mean(r2_full) > mean(r2_red)
    @test mean(r2_full) > 0.5
end
