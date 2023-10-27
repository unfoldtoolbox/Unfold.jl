data, evts = loadtestdata("test_case_3a") #
f = @formula 0 ~ 1 + conditionA + continuousA # 1

# prepare data
data_r = reshape(data, (1, :))
data_r = vcat(data_r, data_r)#add second channel

#--------------------------#
## Mass Univariate Linear ##
#--------------------------#
data_e, times = Unfold.epoch(data = data_r, tbl = evts, τ = (-1.0, 1.9), sfreq = 20) # cut the data into epochs

@testset "test manual pathway" begin
    uf = UnfoldLinearModel(Dict(Any => (f, times)))
    designmatrix!(uf, evts; eventcolumn = "type")
    fit!(uf, data_e)


    @test typeof(uf.modelfit) == Unfold.LinearModelFit
    @test !isempty(coef(uf.modelfit))

end

@testset "epoched auto multi-event" begin
    
    evts_local = deepcopy(evts)
    evts_local.type .= repeat(["A","B"],nrow(evts)÷2)

    uf = fit(UnfoldModel, Dict("A" => (f, times)),evts_local,data_e;eventcolumn="type")
    @test size(coef(uf)) == (2,59,3)
    uf_2events = fit(UnfoldModel, Dict("A" => (f, times),"B"=>(@formula(0~1),times)),evts_local,data_e;eventcolumn="type")
    @test size(coef(uf_2events)) == (2,59,4)

    c = coeftable(uf)
    c2 = coeftable(uf_2events)
    @test c2[c2.basisname .== "event: A",:] == c

    e_uf2 = effects(Dict(:condtionA=>[0,1]),uf_2events)
    e_uf = effects(Dict(:condtionA=>[0,1]),uf)
    
end

@testset "test Autodetection" begin
    @test Unfold.designToModeltype(Dict(Any => (@formula(0 ~ 1), 0:10))) == UnfoldLinearModel
    @test Unfold.designToModeltype(Dict(Any => (@formula(0 ~ 1 + A), 0:10))) ==
        UnfoldLinearModel
    @test Unfold.designToModeltype(
        Dict(Any => (@formula(0 ~ 1 + A), firbasis(τ = (-1, 1), sfreq = 20, name = "basisA"))),
    ) == UnfoldLinearModelContinuousTime
    @test Unfold.designToModeltype(Dict(Any => (@formula(0 ~ 1 + (1 | test)), 0:10))) ==
        UnfoldLinearMixedModel
    @test Unfold.designToModeltype(
        Dict(
            Any => (
                @formula(0 ~ 1 + (1 | test)),
                firbasis(τ = (-1, 1), sfreq = 20, name = "basisA"),
            ),
        ),
    ) == UnfoldLinearMixedModelContinuousTime


    @testset "Bad Input" begin
        # check that if UnfoldLinearModel or UnfoldLinearModelContinuousTime is defined, that the design is appropriate
        basisfunction = firbasis(τ = (-1, 1), sfreq = 20, name = "basisA")
        @test_throws "InputError" fit(UnfoldLinearModel,Dict(Any => (@formula(0 ~ 1),basisfunction)),evts,data_r)
        @test_throws "InputError" fit(UnfoldLinearModel,Dict(Any => (@formula(0 ~ 1),basisfunction)),evts,data_e)
        @test_throws "InputError" fit(UnfoldLinearModelContinuousTime,Dict(Any => (@formula(0 ~ 1),0:0.1:1)),evts,data_r)

    end
end

@testset "automatic, non-dictionary call" begin
    times = -1.0:0.05:1.9
    m_mul = coeftable(fit(UnfoldLinearModel, f, evts, data_e, times))
    
    @test m_mul[(m_mul.channel.==1).&(m_mul.time.==0.1), :estimate] ≈ [2, 3, 4]
    


    data_e_noreshape, times = Unfold.epoch(data = data, tbl = evts, τ = (-1.0, 1.9), sfreq = 20) # cut the data into epochs
    m_mul_noreshape = coeftable(fit(UnfoldLinearModel, f, evts, data_e_noreshape, times))

    @test m_mul_noreshape[
        (m_mul_noreshape.channel.==1).&(m_mul_noreshape.time.==0.1),
        :estimate,
    ] ≈ [2, 3, 4]
    @test size(m_mul_noreshape)[1] == size(m_mul)[1] / 2
    # test 2D call for UnfoldLinearModel
    m_mul_autoreshape = coeftable(fit(UnfoldLinearModel, f, evts, data_e_noreshape[1,:,:], times))
    m_mul_autoreshape == m_mul_noreshape
    # Add Missing in Data
    data_e_missing = deepcopy(data_e)
    data_e_missing[1, 25, end-5:end] .= missing
    m_mul_missing = coeftable(Unfold.fit(UnfoldLinearModel, f, evts, data_e_missing, times))

    @test m_mul_missing.estimate ≈ m_mul.estimate

    # Special solver solver_lsmr_se with Standard Error
    se_solver = solver = (x, y) -> Unfold.solver_default(x, y, stderror = true)
    m_mul_se = coeftable(Unfold.fit(UnfoldModel, f, evts, data_e, times, solver = se_solver))
    @test all(m_mul_se.estimate .≈ m_mul.estimate)
    @test !all(isnothing.(m_mul_se.stderror))

    # robust solver
    rob_solver = (x, y) -> Unfold.solver_robust(x, y)#,rlmOptions=(initial_coef=zeros(3 *length(times)),))
    data_outlier = copy(data_e)
    data_outlier[:,31,1] .= 1000
    m_mul_rob = coeftable(Unfold.fit(UnfoldModel, f, evts, data_outlier, times, solver = rob_solver))
    ix = findall(m_mul_rob.time .≈ 0.5)
    @test all(m_mul_rob.estimate[ix] .≈ m_mul.estimate[ix])

    m_mul_outlier = coeftable(Unfold.fit(UnfoldModel, f, evts, data_outlier, times))
end

#---------------------------------#
## Timexpanded Univariate Linear ##
#---------------------------------#
basisfunction = firbasis(τ = (-1, 1), sfreq = 20, name = "basisA")

@testset "timeexpanded univariate linear+missings" begin
        
    m_tul = coeftable(fit(UnfoldModel, f, evts, data_r, basisfunction))

    @test isapprox(
        m_tul[(m_tul.channel.==1).&(m_tul.time.==0.1), :estimate],
        [2, 3, 4],
        atol = 0.01,
    )

    # test without reshape, i.e. 1 channel vector e.g. size(data) = (1200,)
    m_tul_noreshape = coeftable(fit(UnfoldModel, f, evts, data, basisfunction));
    @test size(m_tul_noreshape)[1] == size(m_tul)[1] / 2

    # Test under missing data
    data_missing = Array{Union{Missing,Number}}(undef, size(data_r))
    data_missing .= deepcopy(data_r)
    data_missing[4500:4600] .= missing

    m_tul_missing = coeftable(fit(UnfoldModel, f, evts, data_missing, basisfunction))

    @test isapprox(m_tul_missing.estimate, m_tul.estimate, atol = 1e-4)  # higher tol because we remove stuff


    ## Test multiple basisfunctions
    b1 = firbasis(τ = (-1, 1), sfreq = 20, name = "basisA")
    b2 = firbasis(τ = (-1, 1), sfreq = 20, name = "basisB")

    f1 = @formula 0 ~ 1 + continuousA # 1
    f2 = @formula 0 ~ 1 + continuousA # 1

    # Fast-lane new implementation
    res = coeftable(
        fit(
            UnfoldModel,
            Dict(0 => (f1, b1), 1 => (f2, b2)),
            evts,
            data_r,
            eventcolumn = "conditionA",
        ),
    )

    # slow manual
    X1 = designmatrix(
        UnfoldLinearModelContinuousTime,
        f1,
        filter(x -> (x.conditionA == 0), evts),
        b1,
    )
    X2 = designmatrix(
        UnfoldLinearModelContinuousTime,
        f2,
        filter(x -> (x.conditionA == 1), evts),
        b2,
    )
    uf = UnfoldLinearModelContinuousTime(Dict(0 => (f1, b1), 1 => (f2, b2)), X1 + X2, [])
    @time fit!(uf, data_r)
    tmp = coeftable(uf)

    # test fast way & slow way to be identical
    @test all(tmp.estimate .== res.estimate)

end

@testset "runtime tests" begin
        
    # runntime tests - does something explode?
    for k = 1:4
        local f
        if k == 1
            f = @formula 0 ~ 1
        elseif k == 2
            f = @formula 0 ~ 1 + conditionA
        elseif k == 3
            f = @formula 0 ~ 0 + conditionA
        elseif k == 4
            f = @formula 0 ~ 1 + continuousA
        end

        fit(UnfoldModel, f, evts, data_e, times)
        fit(UnfoldModel, f, evts, data, basisfunction)
    end
end

@testset "automatic, non-dictionary call" begin
    m_mul = coeftable(fit(UnfoldLinearModel, f, evts, data_e, times))

    @test m_mul[(m_mul.channel.==1).&(m_mul.time.==0.1), :estimate] ≈ [2, 3, 4]
end


@testset "Special solver solver_lsmr_se with Standard Error" begin
    se_solver = solver = (x, y) -> Unfold.solver_default(x, y, stderror = true)
    m_tul_se = coeftable(fit(UnfoldModel, f, evts, data_r, basisfunction, solver = se_solver)) #with
    m_tul = coeftable(fit(UnfoldModel, f, evts, data_r, basisfunction)) #without
    @test all(m_tul_se.estimate .== m_tul.estimate)
    @test !all(isnothing.(m_tul_se.stderror))

    #m_mul_se,m_mul_se = fit(UnfoldLinearModel,f,evts,data_e.+randn(size(data_e)).*5,times,solver=se_solver)
    #plot(m_mul_se[m_mul_se.channel.==1,:],se=true)
    #m_tul_se,m_tul_se = fit(UnfoldLinearModel,f,evts,data_r.+randn(size(data_r)).*5,basisfunction,solver=se_solver)
    #plot(m_tul_se[m_tul_se.channel.==1,:],se=true)

end

@testset "long data" begin
    data_long, evts_long = loadtestdata("test_case_1c") #
    data_long = reshape(data_long, (1, :))

    data_long = vcat(data_long, data_long)
    f_long = @formula 0 ~ 1
    basisfunction_long = firbasis(τ = (-1, 1), sfreq = 1000, name = "basisA")

    @time designmatrix(UnfoldLinearModel, f_long, evts_long, basisfunction_long)
    @time m_tul_long =
        coeftable(fit(UnfoldModel, f_long, evts_long, data_long, basisfunction_long))

    @test isapprox(
        m_tul_long[(m_tul_long.channel.==1).&(m_tul_long.time.==0.1), :estimate],
        [2],
        atol = 0.0001,
    )
    # ~21s, test_case_1c, sfreq = 1000, 6000 events
    #
    ## Older numbers for "designmatrix" only. But see "benchmark/benchmarkjl" for better benchmarks
    # new version 7s-10s, dataset4, sfreq=1000, 1200stim,
    # ~13s, test_case_1c, sfreq = 1000, 6000 events

end

if 1 == 0
    using WGLMakie, AlgebraOfGraphics
    m = mapping(
        :time,
        :estimate,
        color = :colname,
        layout_x = :colname,
        layout_y = :basisname,
    )

    AlgebraOfGraphics.data(df[df.channel.==2, :]) * visual(Lines) * m |> draw

end


##------



if 1 == 0
    # Fit mass-univariate 1st level for all subjects
    basisfunction = firbasis(τ = (-.1, 0.4), sfreq = 10, name = "A")
    resAll = DataFrame()
    f = @formula 0 ~ 1 + condA + condB
    for s in unique(evts.subject)
        ##
        from = minimum(evts[evts.subject.==s, :latency]) - 10
        to = maximum(evts[evts.subject.==s, :latency]) + 40
        to = min(to, size(data, 1))
        evts_s = evts[evts.subject.==s, :]
        evts_s.latency .-= from
        m = fit(
            UnfoldLinearModel,
            f,
            evts_s,
            data[from:to],
            basisfunction,
            contrasts = Dict(:condA => EffectsCoding(), :condB => EffectsCoding()),
        )
        m.results.subject = s
        append!(resAll, m.results)
    end

    results = resAll[resAll.coefname.=="(Intercept)", :]

    results[results.time.==0.1, :]
end
