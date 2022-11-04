module Unfold

using PyMNE # for clusterpermutation
using MixedModelsPermutations

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
using BSplines # for spline predictors

using RobustModels # for robust modelling
#using CategoricalArrays
import StatsBase: fit
import StatsBase: coef
import StatsBase: fit!
import StatsBase: coefnames
import StatsBase: modelmatrix
import StatsModels: width
import StatsModels: terms

import MixedModels.likelihoodratiotest

import Base.(+) # overwrite for DesignMatrices
using Distributions: Gamma, pdf # TODO replace this with direct implementation (used in basisfunction.jl)

include("linearmodels.jl")
include("basisfunctions.jl")
include("designmatrix.jl")
include("fit.jl")
include("utilities.jl")
include("condense.jl")
include("solver.jl")
include("predict.jl")
include("splinepredictors.jl")
include("effects.jl")
include("statistics.jl")
#include("plot.jl") # don't include for now
export fit, fit!, designmatrix!
export firbasis, hrfbasis, condense_long
export UnfoldLinearModel,
    UnfoldLinearMixedModel,
    UnfoldModel,
    UnfoldLinearMixedModelContinuousTime,
    UnfoldLinearModelContinuousTime
export FIRBasis, HRFBasis, SplineBasis
export modelmatrix
export formula, design, designmatrix, coef
export coeftable
export modelfit
export predict
export spl
export cluster_permutation_test

export likelihoodratiotest # statistics.jl
export pvalues # statistics.jl

export effects # effects.jl
import StatsModels.@formula # for exporting
export @formula
end # module
