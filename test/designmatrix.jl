##---

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
    println(Matrix(Xdc))

    @test all(isapprox.(Matrix(Xdc)[1:4, 1:4], shouldBePos, atol = 1e-15))
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
    @test length(events(Xdc)) == 2

    Xdc_3 = Xdc1 + Xdc2 + Xdc2

    @test size(modelmatrix(Xdc_3), 2) ==
          size(modelmatrix(Xdc1), 2) + 2 * size(modelmatrix(Xdc2), 2)
    @test length(events(Xdc_3)) == 3
end
@testset "combining MixedModel Designmatrices" begin

    basisfunction1 = firbasis(τ = (0, 1), sfreq = 10, name = "basis1")
    basisfunction2 = firbasis(τ = (0, 0.5), sfreq = 10, name = "basis2")

    tbl = DataFrame(
        [1 4 10 15 20 22 31 37; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
        [:latency, :subject, :item],
    )
    tbl2 = DataFrame(
        [2 3 12 18 19 25 40 43; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
        [:latency, :subject, :itemB],
    )
    y = Float64.([collect(range(1, stop = 100))...])'
    transform!(tbl, :subject => categorical => :subject)
    transform!(tbl2, :itemB => categorical => :itemB)
    transform!(tbl, :item => categorical => :item)
    #tbl.itemB = tbl.item
    f3 = @formula 0 ~ 1 + (1 | item) + (1 | subject)
    f4 = @formula 0 ~ 1 + (1 | itemB)
    f4_wrong = @formula 0 ~ 1 + (1 | item)
    Xdc3 = designmatrix(UnfoldLinearMixedModel, f3, tbl, basisfunction1)
    Xdc4 = designmatrix(UnfoldLinearMixedModel, f4, tbl2, basisfunction2)
    Xdc4_wrong = designmatrix(UnfoldLinearMixedModel, f4_wrong, tbl, basisfunction2)

    Xdc = Xdc3 + Xdc4
    @test typeof(Xdc.modelmatrix[1][1]) <: SparseArrays.SparseMatrixCSC
    @test length(Xdc.modelmatrix) == 4 # one FeMat  + 3 ReMat
    @test_throws String Xdc3 + Xdc4_wrong
    uf = UnfoldLinearMixedModelContinuousTime(Dict(), Xdc, [])
    Unfold.fit!(uf, y)

    Xdc = Xdc3 + Xdc4
    df = fit(
        UnfoldLinearMixedModelContinuousTime,
        Xdc,
        rand(1, size(Unfold.modelmatrix(Xdc)[1], 1)),
    )
    @test size(Unfold.coef(df), 2) == 17

    ## Speedtest
    if 1 == 0
        x = collect(range(10, stop = 10000, step = 10))
        tbl = DataFrame(reshape(x, 1000, 1), [:latency])
        X = ones(size(tbl, 1), 3) .* [1, 2, 3]'


        basisfunction = firbasis(τ = (0, 1), sfreq = 100, name = "test")
        term = Unfold.TimeExpandedTerm(Term, basisfunction, :latency)
        @time Xdc = Matrix(Unfold.time_expand(X, term, tbl))

    end

end

@testset "equalizeReMatLengths" begin
    bf1 = firbasis(τ = (0, 1), sfreq = 10, name = "basis1")
    bf2 = firbasis(τ = (0, 0.5), sfreq = 10, name = "basis2")

    tbl1 = DataFrame(
        [1 4 10 15 20 22 31 37; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
        [:latency, :subject, :item],
    )
    tbl2 = DataFrame(
        [2 3 12 18 19 25 40 43; 1 1 1 2 2 2 3 3; 1 2 3 1 2 3 1 2]',
        [:latency, :subject, :itemB],
    )

    transform!(tbl1, :subject => categorical => :subject)
    transform!(tbl1, :item => categorical => :item)
    transform!(tbl2, :itemB => categorical => :itemB)
    #tbl.itemB = tbl.item
    f1 = @formula 0 ~ 1 + (1 | item) + (1 | subject)
    f2 = @formula 0 ~ 1 + (1 | itemB)

    form = apply_schema(f1, schema(f1, tbl1), MixedModels.LinearMixedModel)
    form = Unfold.apply_basisfunction(form, bf1, nothing)
    X1 = modelcols.(form.rhs, Ref(tbl1))

    form = apply_schema(f2, schema(f2, tbl2), MixedModels.LinearMixedModel)
    form = Unfold.apply_basisfunction(form, bf2, nothing)
    X2 = modelcols.(form.rhs, Ref(tbl2))

    # no missmatch, shouldnt change anything then
    X = deepcopy(X1[2:end])
    if !isdefined(Base, :get_extension)
        include("../ext/UnfoldMixedModelsExt/UnfoldMixedModelsExt.jl")
        ext = UnfoldMixedModelsExt
    else
        ext = Base.get_extension(Unfold, :UnfoldMixedModelsExt)
    end
    ext.equalize_ReMat_lengths!(X)
    @test all([x[1] for x in size.(X)] .== 48)

    X = (deepcopy(X1[2:end])..., deepcopy(X2[2:end])...)
    @test !all([x[1] for x in size.(X)] .== 49) # not alllenghts the same
    ext.equalize_ReMat_lengths!(X)
    @test all([x[1] for x in size.(X)] .== 49) # now all lengths the same :-)


    X = deepcopy(X2[2])

    @test size(X)[1] == 49
    ext.change_ReMat_size!(X, 52)
    @test size(X)[1] == 52

    X = deepcopy(X2[2])
    @test size(X)[1] == 49
    ext.change_ReMat_size!(X, 40)
    @test size(X)[1] == 40


    X = (deepcopy(X1)..., deepcopy(X2[2:end])...)
    @test size(X[1])[1] == 48
    @test size(X[2])[1] == 48
    @test size(X[3])[1] == 48
    @test size(X[4])[1] == 49
    XA, XB = ext.change_modelmatrix_size!(52, X[1], X[2:end])
    @test size(XA)[1] == 52
    @test size(XB)[1] == 52

    XA, XB = ext.change_modelmatrix_size!(40, X[1], X[2:end])
    @test size(XA)[1] == 40
    @test size(XB)[1] == 40

    XA, XB = ext.changeMatSize!(30, Matrix(X[1]), X[2:end])
    @test size(XA)[1] == 30
    @test size(XB)[1] == 30
end

@testset "Some LinearMixedModel tests" begin

    data, evts = loadtestdata("testCase3", dataPath = (@__DIR__) * "/data") #
    evts.subject = categorical(evts.subject)


    f_zc = @formula 0 ~ 1 + condA + condB + zerocorr(1 + condA + condB | subject)
    basisfunction = firbasis(τ = (-0.1, 0.1), sfreq = 10, name = "ABC")
    Xdc_zc = designmatrix(UnfoldLinearMixedModel, f_zc, evts, basisfunction)

    @test length(Xdc_zc.modelmatrix[2].inds) == 9
    f = @formula 0 ~ 1 + condA + condB + (1 + condA + condB | subject)
    Xdc = designmatrix(UnfoldLinearMixedModel, f, evts, basisfunction)
    @test length(Xdc.modelmatrix[2].inds) == (9 * 9 + 9) / 2

    # test bug with not sequential subjects
    evts_nonseq = copy(evts)
    evts_nonseq = evts_nonseq[.!(evts_nonseq.subject .== 2), :]
    Xdc_nonseq = designmatrix(UnfoldLinearMixedModel, f_zc, evts_nonseq, basisfunction)
    # This used to lead to problems here:
    fit(UnfoldLinearMixedModel, Xdc_nonseq, data')

end

#basisfunction2 = firbasis(τ = (0, 0.5), sfreq = 10, name = "basis2")
@testset "Missings in Events" begin
    tbl = DataFrame(
        :a => [1, 2, 3, 4, 5, 6, 7, 8],
        :b => [1, 1, 1, 2, 2, 2, 3, missing],
        :c => [1, 2, 3, 4, 5, 6, 7, missing],
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
