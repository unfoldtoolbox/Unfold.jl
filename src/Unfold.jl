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
#using DSP
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
 checkFun(sym) = Base.get_extension(@__MODULE__(),sym)
    function solver_robust(args...;kwargs...)
        ext = checkFun(:UnfoldRobustModelsExt) 
        msg = "RobustModels not loaded. Please use ]add RobustModels, using RobustModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.solver_robust(args...;kwargs...)
    end
    function pvalues(args...;kwargs...)
        ext = checkFun(:UnfoldMixedModelsExt) 
        msg = "MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.pvalues(args...;kwargs...)
    end
    function likelihoodratiotest(args...;kwargs...)
        ext = checkFun(:UnfoldMixedModelsExt) 
        msg = "MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.likelihoodratiotest(args...;kwargs...)
    end

    function rePCA(args...;kwargs...)
        ext = checkFun(:UnfoldMixedModelsExt) 
        msg = "MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.rePCA(args...;kwargs...)
    end

    function check_groupsorting(args...;kwargs...)
        ext = checkFun(:UnfoldMixedModelsExt) 
        msg = "MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.check_groupsorting(args...;kwargs...)
    end
    function lmm_combineMats!(args...;kwargs...)
        ext = checkFun(:UnfoldMixedModelsExt) 
        msg = "MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.lmm_combineMats!(args...;kwargs...)
    end
    function splinebasis(args...;kwargs...)
        ext = checkFun(:UnfoldBSplineKitExt) 
        msg = "BSplineKit not loaded. Please use `]add BSplineKit, using BSplineKit` to install/load it, if you want to use splines"
        isnothing(ext) ? throw(msg) : ext.splinebasis(args...;kwargs...)
    end
    


    
    
    function spl(args...;kwargs...)
        ext = checkFun(:UnfoldBSplineKitExt) 
        msg = "BSplineKit not loaded. Please use `]add BSplineKit, using BSplineKit` to install/load it, if you want to use splines"
        isnothing(ext) ? throw(msg) : ext.spl(args...;kwargs...)
    end
    function circspl(args...;kwargs...)
        ext = checkFun(:UnfoldBSplineKitExt) 
        msg = "BSplineKit not loaded. Please use `]add BSplineKit, using BSplineKit` to install/load it, if you want to use splines"
        isnothing(ext) ? throw(msg) : ext.circspl(args...;kwargs...)
    end
  

    function solver_krylov(args...;kwargs...)
        ext = checkFun(:UnfoldKrylovExt) 
        msg = "Krylov or CUDA not loaded. Please use `]add Krylov,CUDA, using Krylov,CUDA` to install/load it, if you want to use GPU-fitting"
        isnothing(ext) ? throw(msg) : ext.solver_krylov(args...;kwargs...)
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
