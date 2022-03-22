@testset "epoch" begin
    d = collect(1:100)
    evt = DataFrame(:latency=>(50))
    
    ep = τ -> Unfold.epoch(d,evt,τ,1)[1][1,:,1]
    # check time delays
    @test ep((0,10.)) ≈ collect(50:60.)
    @test ep((-10,10.)) ≈ collect(40:60.)
    @test ep((-10,0.)) ≈ collect(40:50.)
    @test ep((5,15)) ≈ collect(55:65.)
    @test ep((-15,-5)) ≈ collect(35:45.)

    # check corner cases (sample doesnt end on sampling rate)
    @test ep((0.6,2)) ≈ collect(51:52.)
    @test ep((0.2,2)) ≈ collect(50:52.)
end


#evt = DataFrame(:latency=>(512.7))
#ep = τ -> Unfold.epoch(d,evt,τ,256)[1][1,:,1]

#ep((-0.4, 0.8))
#tau = (-0.2, 0.8); srate = 256 -> 256 samples per epoch -> no error