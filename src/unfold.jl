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
using StaticArrays
using DocStringExtensions
import Base.(+)
using Distributions: Gamma, pdf # TODO replace this with direct implementation (used in basisfunction.jl)

include("basisfunctions.jl")
include("designmatrix.jl")
include("linearmodels.jl")
include("fit.jl")
include("utilities.jl")
include("condense.jl")
#include("plot.jl") # don't include for now
export fit, designmatrix, firbasis,hrfbasis,condense_long,UnfoldLinearModel,UnfoldLinearMixedModel
export unfoldfit # might be renamend to fit! in the future
end # module
