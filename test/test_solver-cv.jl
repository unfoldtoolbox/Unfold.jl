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
