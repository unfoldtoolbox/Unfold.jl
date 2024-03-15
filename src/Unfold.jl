module Unfold

using SimpleTraits

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

#using PooledArrays

#using Tullio
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



#include("plot.jl") # don't include for now
export fit, fit!, designmatrix!
export firbasis, hrfbasis
export AbstractDesignMatrix, AbstractModelFit, UnfoldModel
export UnfoldLinearModel,
    UnfoldLinearMixedModel,
    UnfoldLinearMixedModelContinuousTime,
    UnfoldLinearModelContinuousTime
export FIRBasis, HRFBasis, SplineBasis
export modelmatrix
export formula, design, designmatrix, coef
export coeftable
export modelfit
export predict


if !isdefined(Base, :get_extension)
    ## Extension Compatabality with julia  pre 1.9
    include("../ext/UnfoldRobustModelsExt.jl")
    solver_robust = UnfoldRobustModelsExt.solver_robust

    include("../ext/UnfoldMixedModelsExt/UnfoldMixedModelsExt.jl")
    pvalues = UnfoldMixedModelsExt.pvalues
    using MixedModels
    rePCA = MixedModels.rePCA
    lmm_combine_modelmatrix! = UnfoldMixedModelsExt.lmm_combine_modelmatrix!
    likelihoodratiotest = UnfoldMixedModelsExt.likelihoodratiotest
    check_groupsorting = UnfoldMixedModelsExt.check_groupsorting

    spl() = error("dummy / undefined")
    circspl() = error("dummy / undefined")

    include("../ext/UnfoldBSplineKitExt/UnfoldBSplineKitExt.jl")
    splinebasis = UnfoldBSplineKitExt.splinebasis
    # need this definition unfortunatly
    #circspl = UnfoldBSplineKitExt.circspl
    #spl = UnfoldBSplineKitExt.spl

    include("../ext/UnfoldKrylovExt.jl")
    solver_krylov = UnfoldKrylovExt.solver_krylov
else
    # Jl 1.9+: we need dummy functions, in case the extension was not activated to warn the user if they try to use a functionality that is not yet defined
    checkFun(sym) = Base.get_extension(@__MODULE__(), sym)
    function solver_robust(args...; kwargs...)
        ext = checkFun(:UnfoldRobustModelsExt)
        msg = "RobustModels not loaded. Please use ]add RobustModels, using RobustModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.solver_robust(args...; kwargs...)
    end
    function pvalues(args...; kwargs...)
        ext = checkFun(:UnfoldMixedModelsExt)
        msg = "MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.pvalues(args...; kwargs...)
    end
    function likelihoodratiotest(args...; kwargs...)
        ext = checkFun(:UnfoldMixedModelsExt)
        msg = "MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.likelihoodratiotest(args...; kwargs...)
    end

    function rePCA(args...; kwargs...)
        ext = checkFun(:UnfoldMixedModelsExt)
        msg = "MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.rePCA(args...; kwargs...)
    end

    function check_groupsorting(args...; kwargs...)
        ext = checkFun(:UnfoldMixedModelsExt)
        msg = "MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.check_groupsorting(args...; kwargs...)
    end
    function lmm_combine_modelmatrix!(args...; kwargs...)
        ext = checkFun(:UnfoldMixedModelsExt)
        msg = "MixedModels not loaded. Please use ]add MixedModels, using MixedModels to install it prior to using"
        isnothing(ext) ? throw(msg) : ext.lmm_combine_modelmatrix!(args...; kwargs...)
    end
    function splinebasis(args...; kwargs...)
        ext = checkFun(:UnfoldBSplineKitExt)
        msg = "BSplineKit not loaded. Please use `]add BSplineKit, using BSplineKit` to install/load it, if you want to use splines"
        isnothing(ext) ? throw(msg) : ext.splinebasis(args...; kwargs...)
    end





    function spl(args...; kwargs...)
        ext = checkFun(:UnfoldBSplineKitExt)
        msg = "BSplineKit not loaded. Please use `]add BSplineKit, using BSplineKit` to install/load it, if you want to use splines"
        isnothing(ext) ? throw(msg) : ext.spl(args...; kwargs...)
    end
    function circspl(args...; kwargs...)
        ext = checkFun(:UnfoldBSplineKitExt)
        msg = "BSplineKit not loaded. Please use `]add BSplineKit, using BSplineKit` to install/load it, if you want to use splines"
        isnothing(ext) ? throw(msg) : ext.circspl(args...; kwargs...)
    end


    function solver_krylov(args...; kwargs...)
        ext = checkFun(:UnfoldKrylovExt)
        msg = "Krylov or CUDA not loaded. Please use `]add Krylov,CUDA, using Krylov,CUDA` to install/load it, if you want to use GPU-fitting"
        isnothing(ext) ? throw(msg) : ext.solver_krylov(args...; kwargs...)
    end

end




export spl, circspl

export likelihoodratiotest # statistics.jl
export pvalues # statistics.jl

export effects # effects.jl
import StatsModels.@formula # for exporting
export @formula

export save
export load
end # module
