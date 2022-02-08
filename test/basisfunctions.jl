

firbase = firbasis((-1, 1), 10)

# test optional call
@test Unfold.kernel(firbase) == Unfold.kernel(firbasis((-1, 1), 10))

@test typeof(Unfold.name(firbase)) <: String
# test basics of basisfunction
@test length(collect(Unfold.colnames(firbase))) == 21
@test unique(Unfold.kernel(firbase)(1)) == [1.0, 0.0]

# test length consistency
@test length(Unfold.colnames(firbase)) ==
      length(Unfold.times(firbase)) ==
      size(Unfold.kernel(firbase)(3.1))[2]

# testing the non-sampling rate samples
@test Unfold.kernel(firbase)(0.5)[1:3, 1:3] == [0.5 0.0 0.0; 0.5 0.5 0.0; 0.0 0.5 0.5]


boldbase = hrfbasis(2.0, name = "test")

@test Unfold.name(boldbase) == "test"
@test Unfold.kernel(boldbase)(0) == Unfold.kernel(boldbase)(1)
@test Unfold.kernel(boldbase)(0.1) != Unfold.kernel(boldbase)(1) # testing fractional kernel generation
@test findmax(Unfold.kernel(boldbase)(0.3))[2][1] == 4
