using Test,unfold

firbase = unfold.firbasis(τ=(-1,1),sfreq=10)

# test optional call
@test firbase == unfold.firbasis((-1,1),10)

# test basics of basisfunction
@test length(collect(firbase.times)) == 21
@test unique(firbase.kernel(1))==[1.0,0.0]

# testing the non-sampling rate samples
@test firbase.kernel(0.5)[1:3,1:3] == [0.5 0.0 0.0; 0.5 0.5 0.0; 0.0 0.5 0.5]

@test firbase.exact == true

firbase_notexact =  unfold.firbasis((-1,1),10,false)

# thing is for this test, I'm not sure yet whether I want to have this functionality. Exact could just be a flag.
@test_broken   firbase_notexact.kernel(1) == firbase_notexact.kernel(1.5)
