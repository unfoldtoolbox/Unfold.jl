using Test
import unfold

firbase = firbasis(τ=(-1,1),sfreq=10)

# test optional call
@test firbase.kernel == firbasis((-1,1),10).kernel

@test typeof(firbase.name) <: String
# test basics of basisfunction
@test length(collect(firbase.colnames)) == 21
@test unique(firbase.kernel(1))==[1.0,0.0]

# testing the non-sampling rate samples
@test firbase.kernel(0.5)[1:3,1:3] == [0.5 0.0 0.0; 0.5 0.5 0.0; 0.0 0.5 0.5]


boldbase = hrfbasis(2.,name="test")

@test boldbase.name == "test"
@test boldbase.kernel(0) == boldbase.kernel(1)
@test boldbase.kernel(0.1) != boldbase.kernel(1) # testing fractional kernel generation
@test findmax(boldbase.kernel(0.3))[2][1] == 4
