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
using StatsAPI # for r2

using TimerOutputs # debugging fitting times etc. not strictly needed
#using DSP
using StatsModels
using ProgressMeter
using DocStringExtensions # for Docu
using MLBase # for crossVal

import Term # prettiert output
using OrderedCollections # for Base.Show
import Term: vstack, Panel, tprint
import Term: Tree  # to display Dicts
using PooledArrays
using TypedTables # DataFrames loose the pooled array, so we have to do it differently for now...

using Interpolations # for FIR duration scaling
using ImageTransformations # for FIR duration scaling

#using Tullio
#using BSplineKit # for spline predictors

#using RobustModels # for robust modelling
#using CategoricalArrays
import StatsBase: fit
import StatsBase: coef
import StatsBase: fit!
import StatsBase: coefnames
import StatsBase: modelmatrix
import StatsBase: predict
import StatsModels: width
import StatsModels: terms


import StatsBase.quantile

import Base.show
import Base.+ # overwrite for DesignMatrices
using Distributions: Gamma, pdf # TODO replace this with direct implementation (used in basisfunction.jl)

include("typedefinitions.jl")
include("basisfunctions.jl")
include("timeexpandedterm.jl")
include("designmatrix.jl")
include("fit.jl")
include("utilities.jl")
include("condense.jl")
#include("solver.jl")
include("predict.jl")
include("effects.jl")
include("io.jl")
include("show.jl") # pretty printing

include("solver/main.jl")
include("solver/solvers.jl")
include("solver/prepare.jl")

include("precompile.jl")

#include("plot.jl") # don't include for now
export fit, fit!, designmatrix!
export firbasis, hrfbasis
export AbstractDesignMatrix, AbstractModelFit, UnfoldModel
export UnfoldLinearModel,
    UnfoldLinearMixedModel,
    UnfoldLinearMixedModelContinuousTime,
    UnfoldLinearModelContinuousTime
export FIRBasis, HRFBasis, SplineBasis
export modelmatrix, modelmatrices
export formulas, design, designmatrix, coef
export coeftable, predicttable
export modelfit
export predict, residuals


if !isdefined(Base, :get_extension)
    ## Extension Compatabality with julia  pre 1.9
    include("../ext/UnfoldRobustModelsExt.jl")
    solver_robust = UnfoldRobustModelsExt.solver_robust


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



export effects # effects.jl
import StatsModels.@formula # for exporting
export @formula

export save
export load
end # module
