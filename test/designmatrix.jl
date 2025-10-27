##---
using DataFrames, StatsModels
tbl = DataFrame([1 4]', [:latency])
X = ones(size(tbl))
shouldBeNeg = zeros(4, 4)
shouldBeNeg[1, :] = [1, 0, 0, 1]
shouldBeNeg[2, :] = [0, 1, 0, 0]
shouldBeNeg[3, :] = [0, 0, 1, 0]
shouldBeNeg[4, :] = [0, 0, 0, 1]

shouldBePos = zeros(4, 4)
shouldBePos[1, :] = [0, 0, 0, 0]
shouldBePos[2, :] = [1, 0, 0, 0]
shouldBePos[3, :] = [0, 1, 0, 0]
shouldBePos[4, :] = [0, 0, 1, 0]

#---
@testset "basic designmat" begin
    ## test negative
    basisfunction = firbasis(τ = (-3, 0), sfreq = 1, name = "testing")
    timeexpandterm =
        Unfold.TimeExpandedTerm(FormulaTerm(Term, Term), basisfunction, :latency)
    Xdc = Unfold.time_expand(X, timeexpandterm, tbl)
    @test all(isapprox.(Matrix(Xdc)[1:4, 1:4], shouldBeNeg, atol = 1e-15))

    ## Test Positive only
    basisfunction = firbasis(τ = (1, 4), sfreq = 1, name = "testing")
    timeexpandterm =
        Unfold.TimeExpandedTerm(FormulaTerm(Term, Term), basisfunction, :latency)
    Xdc = Unfold.time_expand(X, timeexpandterm, tbl)


    @test all(isapprox.(Matrix(Xdc)[1:4, 1:4], shouldBePos, atol = 1e-15))
end


@testset "basic designmat interpolated yes/no" begin
    ## test negative
    basisfunction = firbasis(τ = (-3, 0), sfreq = 1, name = "testing")
    timeexpandterm =
        Unfold.TimeExpandedTerm(FormulaTerm(Term, Term), basisfunction, :latency)
    Xdc = Unfold.time_expand(X, timeexpandterm, tbl)
    @test length(Unfold.times(timeexpandterm)) == 4
    @test size(Xdc) == (4, 4)
    @test Unfold.width(basisfunction) == 4
    @test Unfold.height(basisfunction) == 4
    basisfunction = firbasis(τ = (-3, 0), sfreq = 1, name = "testing", interpolate = true)
    timeexpandterm =
        Unfold.TimeExpandedTerm(FormulaTerm(Term, Term), basisfunction, :latency)
    Xdc = Unfold.time_expand(X, timeexpandterm, tbl)
    @test length(Unfold.times(timeexpandterm)) == 5
    @test Unfold.width(basisfunction) == 4
    @test Unfold.height(basisfunction) == 5

    @test size(Xdc) == (5, 4)
end
@testset "customized eventfields" begin
    tbl2 = tbl = DataFrame([1 4]', [:onset])

    timeexpandterm_latency = Unfold.TimeExpandedTerm(FormulaTerm(Term, Term), basisfunction)
    timeexpandterm_onset = Unfold.TimeExpandedTerm(
        FormulaTerm(Term, Term),
        basisfunction,
        eventfields = [:onset],
    )
    Xdc = Unfold.time_expand(X, timeexpandterm_onset, tbl)
    @test_throws ArgumentError Unfold.time_expand(X, timeexpandterm_latency, tbl)
end
@testset "combining designmatrices" begin
    tbl = DataFrame([1 4]', [:latency])
    X = ones(size(tbl))
    basisfunction1 = firbasis(τ = (0, 1), sfreq = 10, name = "basis1")
    basisfunction2 = firbasis(τ = (0, 0.5), sfreq = 10, name = "basis2")
    f = @formula 0 ~ 1
    Xdc1 = designmatrix(UnfoldLinearModelContinuousTime, f, tbl, basisfunction1)
    Xdc2 = designmatrix(UnfoldLinearModelContinuousTime, f, tbl .+ 1, basisfunction2)

    Xdc = Xdc1 + Xdc2
    @test size(modelmatrix(Xdc), 2) ==
          size(modelmatrix(Xdc1), 2) + size(modelmatrix(Xdc2), 2)
    @test length(Unfold.events(Xdc)) == 2

    Xdc_3 = Xdc1 + Xdc2 + Xdc2

    @test size(modelmatrix(Xdc_3), 2) ==
          size(modelmatrix(Xdc1), 2) + 2 * size(modelmatrix(Xdc2), 2)
    @test length(Unfold.events(Xdc_3)) == 3
end

#basisfunction2 = firbasis(τ = (0, 0.5), sfreq = 10, name = "basis2")
@testset "Missings in Events" begin
    tbl = DataFrame(
        :a => [1.0, 2, 3, 4, 5, 6, 7, 8],
        :b => [1.0, 1, 1, 2, 2, 2, 3, missing],
        :c => [1.0, 2, 3, 4, 5, 6, 7, missing],
        :d => ["1", "2", "3", "4", "5", "6", "7", "8"],
        :e => ["1", "2", "3", "4", "5", "6", "7", missing],
        :event => [1, 1, 1, 1, 2, 2, 2, 2],
        :latency => [10, 20, 30, 40, 50, 60, 70, 80],
    )
    tbl.event = string.(tbl.event)
    designmatrix(UnfoldLinearModel, @formula(0 ~ a), tbl)
    @test_throws ErrorException designmatrix(UnfoldLinearModel, @formula(0 ~ a + b), tbl)
    @test_throws ErrorException designmatrix(UnfoldLinearModel, @formula(0 ~ e), tbl)

    # including an actual missing doesnt work
    design = [
        "1" => (@formula(0 ~ a + b + c + d + e), firbasis((0, 1), 1)),
        "2" => (@formula(0 ~ a + b + c + d + e), firbasis((0, 1), 1)),
    ]
    uf = UnfoldLinearModelContinuousTime(design)
    @test_throws ErrorException designmatrix(uf, tbl)

    # but if the missing is in another event, no problem
    design = [
        "1" => (@formula(0 ~ a + b + c + d + e), firbasis((0, 1), 1)),
        "2" => (@formula(0 ~ a + d), firbasis((0, 1), 1)),
    ]
    uf = UnfoldLinearModelContinuousTime(design)
    designmatrix(uf, tbl)

    # prior to the Missing disallow sanity check, this gave an error
    design = [
        "1" => (@formula(0 ~ spl(a, 4) + spl(b, 4) + d + e), firbasis((0, 1), 1)),
        "2" => (@formula(0 ~ a + d), firbasis((0, 1), 1)),
    ]
    uf = UnfoldLinearModelContinuousTime(design)
    designmatrix(uf, tbl)
end
@testset "Integer Covariate Splines" begin
    tbl = DataFrame(
        :a => [1.0, 2, 3, 4, 5, 6, 7, 8],
        :event => [1, 1, 1, 1, 2, 2, 2, 2],
        :latency => [10, 20, 30, 40, 50, 60, 70, 80],
    )
    tbl.event = string.(tbl.event)
    designmatrix(UnfoldLinearModel, @formula(0 ~ a), tbl)

    # prior to the Missing disallow sanity check, this gave an error
    design = ["1" => (@formula(0 ~ spl(a, 4)), firbasis((0, 1), 1))]
    uf = UnfoldLinearModelContinuousTime(design)
    designmatrix(uf, tbl)
end
