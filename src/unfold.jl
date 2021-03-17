module unfold

using PyMNE
using SparseArrays
using StatsModels
using StatsBase
using IterativeSolvers
using DataFrames
using MixedModels
using Missings
using StatsBase
using LinearAlgebra
using Tables # not sure we need it
using GLM # not sure we need it
import MixedModels.FeMat # extended for sparse femats, type piracy => issue on MixedModels.jl github
using TimerOutputs # debugging fitting times etc. not strictly needed
using DSP
using StatsModels
using StaticArrays # for MixedModels extraction of parametrs (inherited from MixedModels.jl, not strictly needed )
using ProgressMeter
using DocStringExtensions # for Docu
using MLBase # for crossVal
using Splines2 # for spline predictors


#using IncompleteLU
import Base.(+) # overwrite for DesignMatrices
using Distributions: Gamma, pdf # TODO replace this with direct implementation (used in basisfunction.jl)

include("basisfunctions.jl")
include("designmatrix.jl")
include("linearmodels.jl")
include("fit.jl")
include("utilities.jl")
include("condense.jl")
include("solver.jl")
include("predict.jl")
include("splinepredictors.jl")
include("clusterpermutation.jl")
#include("plot.jl") # don't include for now
export fit, designmatrix, firbasis,hrfbasis,condense_long,UnfoldLinearModel,UnfoldLinearMixedModel
export unfoldfit # might be renamend to fit! in the future
export predict,spl
export cluster_permutation_test 
end # module
