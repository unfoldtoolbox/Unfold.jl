using BenchmarkTools
using Random

const SUITE = BenchmarkGroup()
Random.seed!(3)



SUITE = BenchmarkGroup()
SUITE["nodc"] = BenchmarkGroup(
    ["nodc"])
SUITE["dc"] = BenchmarkGroup(["dc"])
# designmatrix generation
SUITE["dc"]["X_gen_lin"] = @benchmarkable read(run(`git diff Project.toml`)) 
SUITE["nodc"]["fit_lin"] = @benchmarkable 1==2
read(run(`git status`))
