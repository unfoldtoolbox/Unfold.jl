
@testset "lmm tests" begin
    ###############################
    ##  Mixed Model tests
    ###############################
    data, evts = loadtestdata("testCase3", dataPath = (@__DIR__) * "/data") #
    append!(data, zeros(1000))
    data = reshape(data, 1, :)
    data = vcat(data, data)
    data = data .+ 1 * randn(size(data)) # we have to add minimal noise, else mixed models crashes.
    data_missing = Array{Union{Missing,Number}}(undef, size(data))
    data_missing .= deepcopy(data)

    data_missing[4500:4600] .= missing

    transform!(evts, :subject => categorical => :subject)

    f = @formula 0 ~ 1 + condA + condB + (1 + condA + condB | subject)
    #f  = @formula 0~1 + (1|subject)



    # cut the data into epochs
    # TODO This ignores subject bounds
    data_e, times = Unfold.epoch(data = data, tbl = evts, τ = (-1.0, 1.9), sfreq = 10)
    data_missing_e, times =
        Unfold.epoch(data = data_missing, tbl = evts, τ = (-1.0, 1.9), sfreq = 10)
    evts_e, data_e = Unfold.drop_missing_epochs(copy(evts), data_e)
    evts_missing_e, data_missing_e = Unfold.drop_missing_epochs(copy(evts), data_missing_e)

    ######################
    ##  Mass Univariate Mixed
    @time m_mum = fit(
        UnfoldModel,
        f,
        evts_e,
        data_e,
        times,
        contrasts = Dict(:condA => EffectsCoding(), :condB => EffectsCoding()),
        show_progress = false,
    )
    df = Unfold.coeftable(m_mum)
    @test isapprox(
        df[(df.channel.==1).&(df.coefname.=="condA: 1").&(df.time.==0.0), :estimate],
        [5.618, 9.175],
        rtol = 0.1,
    )



    # with missing
    @time m_mum = fit(
        UnfoldModel,
        f,
        evts_missing_e,
        data_missing_e,
        times,
        contrasts = Dict(:condA => EffectsCoding(), :condB => EffectsCoding()),
        show_progress = false,
    )
    df = coeftable(m_mum)
    @test isapprox(
        df[(df.channel.==1).&(df.coefname.=="condA: 1").&(df.time.==0.0), :estimate],
        [5.618, 9.175],
        rtol = 0.1,
    )


    # Timexpanded Univariate Mixed
    f = @formula 0 ~ 1 + condA + condB + (1 + condA | subject)
    basisfunction = firbasis(τ = (-0.2, 0.3), sfreq = 10)
    @time m_tum = fit(
        UnfoldModel,
        f,
        evts,
        data,
        basisfunction,
        contrasts = Dict(:condA => EffectsCoding(), :condB => EffectsCoding()),
        show_progress = false,
    )
    df = coeftable(m_tum)
    @test isapprox(
        df[(df.channel.==1).&(df.coefname.=="condA: 1").&(df.time.==0.0), :estimate],
        [5.618, 9.175],
        rtol = 0.1,
    )


    # missing data in LMMs
    # not yet implemented
    Test.@test_broken m_tum = fit(
        UnfoldModel,
        f,
        evts,
        data_missing,
        basisfunction,
        contrasts = Dict(:condA => EffectsCoding(), :condB => EffectsCoding()),
    )


    evts.subjectB = evts.subject
    evts1 = evts[evts.condA.==0, :]
    evts2 = evts[evts.condA.==1, :]

    f0_lmm = @formula 0 ~ 1 + condB + (1 | subject) + (1 | subjectB)
    @time m_tum = coeftable(
        fit(UnfoldModel, f0_lmm, evts, data, basisfunction; show_progress = false),
    )


    f1_lmm = @formula 0 ~ 1 + condB + (1 | subject)
    f2_lmm = @formula 0 ~ 1 + condB + (1 | subjectB)

    b1 = firbasis(τ = (-0.2, 0.3), sfreq = 10, name = 0)
    b2 = firbasis(τ = (-0.1, 0.3), sfreq = 10, name = 1)

    ext = Base.get_extension(Unfold, :UnfoldMixedModelsExt)
    X1_lmm = designmatrix(ext.UnfoldLinearMixedModelContinuousTime, f1_lmm, evts1, b1)
    X2_lmm = designmatrix(ext.UnfoldLinearMixedModelContinuousTime, f2_lmm, evts2, b2)

    r = fit(
        ext.UnfoldLinearMixedModelContinuousTime,
        X1_lmm + X2_lmm,
        data;
        show_progress = false,
    )
    df = coeftable(r)

    @test isapprox(
        df[(df.channel.==1).&(df.coefname.=="condB").&(df.time.==0.0), :estimate],
        [18.21, 17.69],
        rtol = 0.1,
    )

    # Fast-lane new implementation
    m = coeftable(
        fit(
            UnfoldModel,
            [0 => (f1_lmm, b1), 1 => (f2_lmm, b2)],
            evts,
            data,
            eventcolumn = "condA",
        ),
    )

end
## Condense check for multi channel, multi
@testset "LMM multi channel, multi basisfunction" begin
    data, evts = loadtestdata("testCase3", dataPath = (@__DIR__) * "/data")
    transform!(evts, :subject => categorical => :subject)
    data = vcat(data', data')

    bA0 = firbasis(τ = (-0.0, 0.1), sfreq = 10, name = 0)
    bA1 = firbasis(τ = (0.1, 0.2), sfreq = 10, name = 1)
    evts.subject2 = evts.subject
    fA0 = @formula 0 ~ 1 + condB + zerocorr(1 | subject)
    fA1 = @formula 0 ~ 1 + condB + zerocorr(1 | subject2)
    m = fit(
        UnfoldModel,
        Dict(0 => (fA0, bA0), 1 => (fA1, bA1)),
        evts,
        data;
        eventcolumn = "condA",
        show_progress = false,
    )

    res = coeftable(m)

    @test all(last(.!isnothing.(res.group), 8))
    @test all(last(res.coefname, 8) .== "(Intercept)")
end


@testset "LMM bug reorder #115" begin

    data, evts = UnfoldSim.predef_2x2(;
        return_epoched = true,
        n_subjects = 10,
        noiselevel = 1,
        onset = NoOnset(),
    )

    data = reshape(data, size(data, 1), :)

    designList = [
        [
            Any => (
                @formula(
                    0 ~
                        1 + A + B + zerocorr(1 + B + A | subject) + zerocorr(1 + B | item)
                ),
                range(0, 1, length = size(data, 1)),
            ),
        ],
        [
            Any => (
                @formula(
                    0 ~
                        1 + A + B + zerocorr(1 + A + B | subject) + zerocorr(1 + B | item)
                ),
                range(0, 1, length = size(data, 1)),
            ),
        ],
        [
            Any => (
                @formula(0 ~ 1 + zerocorr(1 + A + B | subject) + zerocorr(1 | item)),
                range(0, 1, length = size(data, 1)),
            ),
        ],
    ]
    #des = designList[1]
    #des = designList[2]
    for des in designList
        @test_throws AssertionError fit(UnfoldModel, des, evts, data)
        #
    end

    #counter check

    des = [
        Any => (
            @formula(0 ~ 1 + zerocorr(1 | item) + zerocorr(1 + A + B | subject)),
            range(0, 1, length = size(data, 1)),
        ),
    ]

    #= fails but not in the repl...?
    uf = fit(UnfoldModel, des, evts, data; show_progress = false)
    @test 3 ==
          unique(
        @subset(
            coeftable(uf),
            @byrow(:group == Symbol("subject")),
            @byrow :time == 0.0
        ).coefname,
    ) |> length
    =#
end


@testset "LMM bug reshape #110" begin
    data, evts =
        UnfoldSim.predef_2x2(; return_epoched = true, n_subjects = 10, noiselevel = 1)
    data = reshape(data, size(data, 1), :)

    des = [
        Any => (
            @formula(
                0 ~ 1 + A + B + zerocorr(1 + B + A | item) + zerocorr(1 + B | subject)
            ),
            range(0, 1, length = size(data, 1)),
        ),
    ]
    uf = fit(UnfoldModel, des, evts, data; show_progress = false)
    @test size(coef(uf)) == (1, 100, 3)
end
