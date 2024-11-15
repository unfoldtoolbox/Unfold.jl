
@testset "FIR" begin
    firbase = firbasis((-1, 1), 10)

    # test optional call
    @test Unfold.kernel(firbase, 1) == Unfold.kernel(firbasis((-1, 1), 10), 1)

    @test typeof(Unfold.name(firbase)) <: String
    # test basics of basisfunction
    @test length(collect(Unfold.colnames(firbase))) == 21
    @test unique(Unfold.kernel(firbase, 1)) == [1.0, 0.0]

    # test length consistency
    @test length(Unfold.colnames(firbase)) == size(Unfold.kernel(firbase, 3.1))[2]
    @test length(Unfold.times(firbase)) == size(Unfold.kernel(firbase, 3.1))[1]

    # testing the non-sampling rate samples

    # test the interpolate / true false
    firbase_off = firbasis((-1, 1), 10; interpolate = false)
    firbase_on = firbasis((-1, 1), 10; interpolate = true)

    @test Unfold.kernel(firbase_off, 0.5)[1:3, 1:3] == [1 0 0; 0 1 0; 0 0 1]
    @test Unfold.kernel(firbase_on, 0.5)[1:3, 1:3] ==
          [0.5 0.0 0.0; 0.5 0.5 0.0; 0.0 0.5 0.5]

    # defining a firkernel with a duration = 1 sample, should be equal to the non-duration one
    f_dur = Unfold.firkernel([103.3; 1], range(-0.1, step = 0.01, stop = 0.31))
    f_fir = Unfold.firkernel(103.3, range(-0.1, step = 0.01, stop = 0.31))
    @test f_dur == f_fir

    # test duration for samples = 4
    f_dur = Unfold.kernel(firbase, [1, 4])
    @test all(sum(f_dur, dims = 1) .== 4)
end
@testset "BOLD" begin

    boldbase = hrfbasis(2.0, name = "test")

    @test Unfold.name(boldbase) == "test"
    @test Unfold.kernel(boldbase, 0) == Unfold.kernel(boldbase, 1)
    @test Unfold.kernel(boldbase, 0.1) != Unfold.kernel(boldbase, 1) # testing fractional kernel generation
    @test findmax(Unfold.kernel(boldbase, 0.3))[2][1] == 4


end

@testset "timespline" begin
    splinebase = Unfold.splinebasis(τ = (-1, 1), sfreq = 20, nsplines = 10, name = "basisA")

    @test length(Unfold.colnames(splinebase)) == size(Unfold.kernel(splinebase, 3.1))[2]
    @test length(Unfold.times(splinebase)) == size(Unfold.kernel(splinebase, 3.1))[1]

end
