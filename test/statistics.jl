@testset "LMM LRT" begin
    data, evts = UnfoldSim.predef_eeg(10; n_items = 20, sfreq = 10, return_epoched = true)
    data = reshape(data, 1, size(data, 1), size(data, 2) * size(data, 3)) # add second channel
    data = vcat(data, data)

    times = range(0, 1, size(data, 2))

    f0 = @formula 0 ~ 1 + condition + (1 | subject)
    f1 = @formula 0 ~ 1 + condition + continuous + (1 | subject)

    m0 = fit(UnfoldModel, [Any => (f0, times)], evts, data)
    m1 = fit(UnfoldModel, [Any => (f1, times)], evts, data)

    uf_lrt = likelihoodratiotest(data, m0, m1)
    uf_pval = reshape(vcat(pvalues(uf_lrt)...), size(data, 2), 2)'
    @test uf_pval[1, :] == uf_pval[2, :]

    #---
    for tix = 1:3
    evts[!, :y] = data[1, tix, :]

    f0 = @formula y ~ 1 + condition + (1 | subject)
    f1 = @formula y ~ 1 + condition + continuous + (1 | subject)

    lmm0 = fit(MixedModel, f0, evts)
    lmm1 = fit(MixedModel, f1, evts)




        mm_lrt = MixedModels.likelihoodratiotest(lmm0, lmm1)
        #---
        @test isapprox(mm_lrt.pvalues[1], uf_pval[1, tix]; atol = 1e-10)
        #---
    end

end
