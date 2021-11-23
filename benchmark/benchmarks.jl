using BenchmarkTools
using Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.add(["BenchmarkTools", "PkgBenchmark","Random"])
Pkg.resolve()

using Random
const SUITE = BenchmarkGroup()
Random.seed!(3)


SUITE = BenchmarkGroup()
SUITE["nodc"] = BenchmarkGroup(
    ["nodc"])
SUITE["dc"] = BenchmarkGroup(["dc"])
# designmatrix generation
SUITE["dc"]["X_gen_lin"] = @benchmarkable 1==1
SUITE["nodc"]["fit_lin"] = @benchmarkable 1==2
