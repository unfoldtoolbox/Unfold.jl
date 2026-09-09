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


@testset "detectbad_peak_to_peak" begin
    # The following tests were created using LLM (but checked)

    # ---------- clean signals -------------------------------------------------
    @testset "clean signal -> nothing flagged" begin
        sig = 10.0 .* sin.(2π .* (1:1000) ./ 100)        # |sig| ≤ 10  =>  ptp ≤ 20
        mask = detectbad_peak_to_peak(
            sig;
            threshold = 50.0,
            sfreq = 100,
            window = 0.05,
            stepsize = 0.01,
        )  # w = 5
        @test mask isa BitVector
        @test size(mask) == size(sig)
        @test count(mask) == 0
    end

    @testset "return type and length" begin
        mask = detectbad_peak_to_peak(
            Float64[0.0, 10.0, 20.0];
            sfreq = 100,
            window = 0.02,
            stepsize = 0.01,
        )  # w = 1
        @test length(mask) == 3
        @test count(mask) == 0
    end

    # ---------- single artifact -----------------------------------------------
    @testset "a spike flags exactly the union of the windows covering it" begin
        sig = zeros(100)
        sig[50] = 100.0
        mask = detectbad_peak_to_peak(
            sig;
            threshold = 50.0,
            sfreq = 100,
            window = 0.05,
            stepsize = 0.01,
        )     # w = 5
        # windows covering index 50: i = 46..50  ->  union of samples 46..54
        @test findall(mask) == collect(46:54)
        @test count(mask) == 9              # 9 samples flagged (not just the peak)
        @test count(mask[1:45]) == 0
        @test count(mask[55:end]) == 0
    end

    # ---------- strict '>' semantics ------------------------------------------
    @testset "ptp == threshold is NOT flagged" begin
        sig = zeros(10)
        sig[2] = 50.0                        # ptp = 50.0 == threshold
        mask = detectbad_peak_to_peak(
            sig;
            threshold = 50.0,
            sfreq = 100,
            window = 0.02,
            stepsize = 0.01,
        )     # w = 2
        @test count(mask) == 0
    end

    @testset "ptp just above threshold IS flagged" begin
        sig = zeros(10)
        sig[2] = 50.1
        mask = detectbad_peak_to_peak(
            sig;
            threshold = 50.0,
            sfreq = 100,
            window = 0.02,
            stepsize = 0.01,
        )     # w = 2
        @test findall(mask) == collect(1:3)   # windows covering index 2: i = 1..2
    end

    # ---------- default threshold ---------------------------------------------
    @testset "threshold defaults to 150.0" begin
        sig = zeros(20)
        sig[5] = 151.0
        mask = detectbad_peak_to_peak(sig; sfreq = 100, window = 0.05, stepsize = 0.01)  # w = 5, threshold = 50 (default)
        @test findall(mask) == collect(1:9)
    end

    # ---------- window size ---------------------------------------------------
    @testset "larger window flags a wider span" begin
        sig = zeros(100)
        sig[50] = 100.0
        m5 = detectbad_peak_to_peak(
            sig;
            threshold = 50.0,
            sfreq = 100,
            window = 0.05,
            stepsize = 0.01,
        )      # w = 5
        m10 = detectbad_peak_to_peak(
            sig;
            threshold = 50.0,
            sfreq = 100,
            window = 0.10,
            stepsize = 0.01,
        )      # w = 10
        @test findall(m5) == collect(46:54)
        @test findall(m10) == collect(41:59)
        @test count(m10) > count(m5)
    end

    @testset "window length = round(sfreq * window)" begin
        sig = zeros(10)
        sig[5] = 100.0
        mask = detectbad_peak_to_peak(
            sig;
            threshold = 50.0,
            sfreq = 100,
            window = 0.02,
            stepsize = 0.01,
        )     # 2.0 -> w = 2
        @test findall(mask) == collect(4:6)
    end

    @testset "window clamped to 1 sample never flags" begin
        sig = [0.0, 1000.0, 0.0]
        mask = detectbad_peak_to_peak(
            sig;
            threshold = 50.0,
            sfreq = 100,
            window = 0.002,
            stepsize = 0.001,
        )    # round(0.2)=0 -> w = 1
        @test count(mask) == 0
    end

    # ---------- multiple artifacts --------------------------------------------
    @testset "two separate spikes -> two flagged regions" begin
        sig = zeros(100)
        sig[20] = 100.0
        sig[80] = 100.0
        mask = detectbad_peak_to_peak(
            sig;
            threshold = 50.0,
            sfreq = 100,
            window = 0.05,
            stepsize = 0.01,
        )     # w = 5
        @test findall(mask) == vcat(16:24, 76:84)
    end

    # ---------- edges / robustness --------------------------------------------
    @testset "signal shorter than window -> all false, no error" begin
        sig = [0.0, 100.0, 0.0]
        mask = detectbad_peak_to_peak(
            sig;
            threshold = 50.0,
            sfreq = 100,
            window = 0.1,
            stepsize = 0.01,
        )      # w = 10 > n = 3
        @test count(mask) == 0
        @test size(mask) == (3,)
    end

    @testset "accepts non-Float64 (Integer) input" begin
        sig = Int[0, 0, 100, 0, 0]
        mask = detectbad_peak_to_peak(
            sig;
            threshold = 50,
            sfreq = 100,
            window = 0.05,
            stepsize = 0.01,
        )       # w = 5, n = 5
        @test findall(mask) == collect(1:5)
    end

    @testset "sfreq is typed ::Int (a Float throws)" begin
        @test_throws TypeError detectbad_peak_to_peak(
            Float64[0.0, 100.0];
            sfreq = 100.0,
            window = 0.05,
            stepsize = 0.01,
        )
    end

    # ---------- differential / property check ---------------------------------
    @testset "matches a brute-force reference" begin
        rng = MersenneTwister(1234)
        sig = 5.0 .* randn(rng, 400)
        sig[100] = 500.0
        sig[250] = -400.0              # guarantee artifacts
        sfreq, window, threshold = 100, 0.04, 30.0
        w = round(Int, sfreq * window)                     # = 4

        got = detectbad_peak_to_peak(
            sig;
            threshold = threshold,
            sfreq = sfreq,
            window = window,
            stepsize = 0.01,
        )

        ref = falses(length(sig))
        for i = 1:(length(sig)-w+1)
            seg = sig[i:(i+w-1)]
            if maximum(seg) - minimum(seg) > threshold
                ref[i:(i+w-1)] .= true
            end
        end
        @test got == ref
    end

    # ---------- matrix overload (rows = channels) -----------------------------
    @testset "matrix applies the per-channel mask to each row" begin
        sig = zeros(2, 20)                                 # 2 channels × 20 samples
        sig[1, 5] = 100.0                                  # spike in channel 1 @ sample 5
        mask = detectbad_peak_to_peak(
            sig;
            threshold = 50.0,
            sfreq = 100,
            window = 0.05,
            stepsize = 0.01,
        )     # w = 5
        @test size(mask) == (2, 20)
        @test findall(mask[1, :]) == collect(1:9)
        @test count(mask[2, :]) == 0
    end

end
