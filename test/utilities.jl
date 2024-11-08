@testset "epoch" begin
    d = collect(1:100)
    evt = DataFrame(:latency => (50))

    ep = τ -> Unfold.epoch(d, evt, τ, 1)[1][1, :, 1]
    # check time delays
    @test ep((0, 10.0)) ≈ collect(50:60.0)
    @test ep((-10, 10.0)) ≈ collect(40:60.0)
    @test ep((-10, 0.0)) ≈ collect(40:50.0)
    @test ep((5, 15)) ≈ collect(55:65.0)
    @test ep((-15, -5)) ≈ collect(35:45.0)

    # check corner cases (sample doesnt end on sampling rate)
    @test ep((0.6, 2)) ≈ collect(51:52.0)
    @test ep((0.2, 2)) ≈ collect(50:52.0)


    # test sampling frequencies

    ep = τ -> Unfold.epoch(d, evt, τ, 2)[1][1, :, 1]
    @test ep((-1.0, 2)) ≈ collect(48:54.0)

    ep = τ -> Unfold.epoch(d, evt, τ, 0.5)[1][1, :, 1]
    @test ep((-4.0, 8)) ≈ collect(48:54.0)

    # rounding bug when latency was .5 -> bug #78
    d = zeros((1, 1270528))
    evt = DataFrame(:latency => (181603.5))
    ep = τ -> Unfold.epoch(d, evt, τ, 256.0)[1][1, :, 1]
    ep((-0.1, 0.8))
end
