using Unfold

include("setup.jl")


# sanity checks, auto-quality control
#using Aqua
#Aqua.test_all(Unfold)

@testset "BasisFunctions" begin
    include("basisfunctions.jl")
end

@testset "Fitting" begin
    include("fit.jl")
end

@testset "Designmatrix" begin
    include("designmatrix.jl")
end

@testset "Splines" begin
    include("splines.jl")
end

@testset "Predict" begin
    include("predict.jl")
end

@testset "Effects" begin
    include("effects.jl")
end


@testset "Utilities" begin
    include("utilities.jl")
end

@testset "IO" begin
    include("io.jl")
end
#@testset "ClusterPermutation" begin
#include("clusterpermutation.jl")
#end
