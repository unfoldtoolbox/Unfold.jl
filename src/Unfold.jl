module Unfold

using SparseArrays
using StatsModels
using StatsBase
using IterativeSolvers
using DataFrames
#using MixedModels
using Missings
using StatsBase
using LinearAlgebra
using Tables # not sure we need it
using GLM # not sure we need it

using TimerOutputs # debugging fitting times etc. not strictly needed
using DSP
using StatsModels
using ProgressMeter
using DocStringExtensions # for Docu
using MLBase # for crossVal
#using BSplineKit # for spline predictors

#using RobustModels # for robust modelling
#using CategoricalArrays
import StatsBase: fit
import StatsBase: coef
import StatsBase: fit!
import StatsBase: coefnames
import StatsBase: modelmatrix
import StatsModels: width
import StatsModels: terms



import StatsBase.quantile

import Base.show
import Base.(+) # overwrite for DesignMatrices
using Distributions: Gamma, pdf # TODO replace this with direct implementation (used in basisfunction.jl)

include("typedefinitions.jl")
include("basisfunctions.jl")
include("timeexpandedterm.jl")
include("designmatrix.jl")
include("fit.jl")
include("utilities.jl")
include("condense.jl")
include("solver.jl")
include("predict.jl")
include("splinepredictors.jl")
include("effects.jl")
include("statistics.jl")
include("io.jl")


## Extension Compatabality with julia  pre 1.9

if !isdefined(Base, :get_extension)
    include("../ext/UnfoldRobustModelsExt.jl")
    solver_robust = UnfoldRobustModelsExt.solver_robust

    include("../ext/UnfoldMixedModelsExt.jl")
    pvalues = UnfoldMixedModelsExt.pvalues
    likelihoodratiotest = UnfoldMixedModelsExt.likelihoodratiotest
else

    function solver_robust(args...;kwargs...)
        ext = Base.get_extension(@__MODULE__(),:UnfoldRobustModelsExt)
        if ext !== nothing
            return ext.solver_robust(args...;kwargs...)
        else
            throw_error("RobustModels not loaded. Please use ]add RobustModels, using RobustModels to install it prior to using")
        end
    end
    function pvalues(args...;kwargs...)
        ext = Base.get_extension(@__MODULE__(),:UnfoldMixedModelsExt)
        if ext !== nothing
            return ext.pvalues(args...;kwargs...)
        else
            throw_error("MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using")
        end
    end

    function likelihoodratiotest(args...;kwargs...)
        ext = Base.get_extension(@__MODULE__(),:UnfoldMixedModelsExt)
        if ext !== nothing
            return ext.likelihoodratiotest(args...;kwargs...)
        else
            throw_error("MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using")
        end
    end
end


#include("plot.jl") # don't include for now
export fit, fit!, designmatrix!
export firbasis, hrfbasis
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
export spl,circspl

export likelihoodratiotest # statistics.jl
export pvalues # statistics.jl

export effects # effects.jl
import StatsModels.@formula # for exporting
export @formula

export save
export load
end # module
