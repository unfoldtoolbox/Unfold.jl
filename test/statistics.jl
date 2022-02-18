@testset "LMM LRT" begin
    data,evts = loadtestdata("testCase3",dataPath="data");
    
        evts[!,:subject] .= string.(evts.subject);
        data_epoch,times = Unfold.epoch(data=data.+1.1.*randn(size(data)),tbl=evts,τ=(-0.1,0.3),sfreq=10);
        evts_epoch,data_epoch = Unfold.dropMissingEpochs(evts,data_epoch)
   
        
        f0  = @formula 0~1+condA + (1|subject);
        f1  = @formula 0~1+condA+condB + (1|subject);
            
        m1 = fit(UnfoldModel,Dict(Any=>(f0,times)),evts_epoch,data_epoch);
        m2 = fit(UnfoldModel,Dict(Any=>(f1,times)),evts_epoch,data_epoch);
        
        evts[!,:y] = data_epoch[1,1,:]

        f0  = @formula y~1+condA + (1|subject);
        f1  = @formula y~1+condA+condB + (1|subject);
        
        lmm0 = fit(MixedModel,f0,evts)
        lmm1 = fit(MixedModel,f1,evts)

        uf_lrt = likelihoodratiotest(m1,m2)
        uf_lrt[1] == MixedModels.likelihoodratiotest(lmm0,lmm1)


        @test MixedModels.likelihoodratiotest(lmm0,lmm1).pvalues ≈ uf_lrt[1].pvalues


        
        
end