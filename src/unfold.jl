__precompile__(false)

module unfold

using SparseArrays
using StatsModels
using StatsBase
using IterativeSolvers
using DataFrames
using MixedModels
using StatsBase
using LinearAlgebra
using Tables
using GLM
import MixedModels.FeMat
using TimerOutputs
using DSP
using StatsModels
import Base.(+)
using Distributions: Gamma, pdf # TODO replace this with direct implementation (used in basisfunction.jl)
include("linearmodels.jl")
include("basisfunctions.jl")
include("designmatrix.jl")
include("fit.jl")
include("utilities.jl")
include("condense.jl")
#include("plot.jl") # don't include for now
end # module
