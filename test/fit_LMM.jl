
@testset "LMM bug reorder #115" begin
        
    data,evts = UnfoldSim.predef_2x2(;return_epoched=true,n_subjects=10,noiselevel=1)

designList = [Dict(Any=>(@formula(0~1+A+B+zerocorr(1+B+A|subject)+zerocorr(1+B|item)),range(0,1,length=size(data,1)))),
              Dict(Any=>(@formula(0~1+A+B+zerocorr(1+A+B|subject)+zerocorr(1+B|item)),range(0,1,length=size(data,1)))),
              Dict(Any=>(@formula(0~1+zerocorr(1+A+B|subject)+zerocorr(1|item)),range(0,1,length=size(data,1))))]
#des = designList[1]
#des = designList[2]
    for des = designList
        @test_throws AssertionError fit(UnfoldModel,des,evts,data)
        #
    end

    #counter check

    des = Dict(Any=>(@formula(0~1+zerocorr(1|item)+zerocorr(1+A+B|subject)),range(0,1,length=size(data,1))))
    uf = fit(UnfoldModel,des,evts,data)
    @test 3 == unique(@subset(coeftable(uf),@byrow(:group == Symbol("subject")),@byrow :time == 0.0).coefname) |> length
end

@testset "LMM bug reshape #110" begin
    des =    Dict(Any=>(@formula(0~1+A+B+zerocorr(1+B+A|item)+zerocorr(1+B|subject)),range(0,1,length=size(data,1))))
    data,evts = UnfoldSim.predef_2x2(;return_epoched=true,n_subjects=10,noiselevel=1)
    uf = fit(UnfoldModel,des,evts,data)
    @test size(coef(uf)) ==(1,100,3)
end