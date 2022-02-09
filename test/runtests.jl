include("setup.jl")


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

#@testset "ClusterPermutation" begin
#include("clusterpermutation.jl")
#end
