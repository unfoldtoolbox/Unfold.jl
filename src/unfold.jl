__precompile__(false)

module unfold

using SparseArrays
using StatsModels
using IterativeSolvers
using DataFrames
using MixedModels
using StatsBase
using LinearAlgebra
using Tables
import MixedModels.FeMat
include("linearmodels.jl")
include("basisfunctions.jl")
include("designmatrix.jl")
include("unfoldfit.jl")
include("utilities.jl")
include("condense.jl")
end # module
