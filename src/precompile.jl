using PrecompileTools: @setup_workload, @compile_workload




@setup_workload begin
    data = rand(Float64, 2, 100)
    evts = DataFrame(
        :latency => (range(10, 70, length = 3)),
        :condition => [:a, :a, :b],
        :continuous => 1:3,
    )


    @compile_workload begin

        # cut the data into epochs
        data_epochs, t = Unfold.epoch(data = data, tbl = evts, τ = (-0.1, 0.1), sfreq = 10) # channel x timesteps x trials
        f = @formula 0 ~ 1 + condition + continuous # note the formulas left side is `0 ~ ` for technical reasons`


        basisfunction = firbasis(τ = (-0.1, 0.1), sfreq = 10)
        bf_vec = [Any => (f, basisfunction)]


        for m in [
            @suppress(fit(UnfoldModel, f, evts, data_epochs, t)),
            @suppress(fit(UnfoldModel, [Any => (f, t)], evts, data_epochs)),
            @suppress(
                fit(
                    UnfoldModel,
                    [Any => (f, t)],
                    Unfold.drop_missing_epochs(evts, data_epochs)[1],
                    Unfold.drop_missing_epochs(evts, data_epochs)[2],
                )
            ),
            @suppress(fit(UnfoldModel, bf_vec, evts, data)),
        ]
            coef(m)
            coeftable(m)
            effects(Dict(:condition => [:a, :b]), m)
            effects(Dict(:continuous => 1:2), m)
            effects(Dict(:continuous => collect(1:2)), m)

        end

    end
end
