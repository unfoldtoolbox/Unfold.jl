using Test
import unfold

firbase = unfold.firbasis(τ=(-1,1),sfreq=10)

# test optional call
@test firbase == unfold.firbasis((-1,1),10)

@test firbase.name == ""
# test basics of basisfunction
@test length(collect(firbase.colnames)) == 21
@test unique(firbase.kernel(1))==[1.0,0.0]

# testing the non-sampling rate samples
@test firbase.kernel(0.5)[1:3,1:3] == [0.5 0.0 0.0; 0.5 0.5 0.0; 0.0 0.5 0.5]
