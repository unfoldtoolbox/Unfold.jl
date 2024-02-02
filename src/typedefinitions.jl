
"""
The main abstract model-type of the toolbox. E.g. `UnfoldLinearModel` is a concrete type of this
"""
abstract type UnfoldModel end


"""
Abstract Type to report modelresults
"""
abstract type ModelFit end



"""
    DesignMatrix
Type that keeps an Array of  `formulas`, designmatrices `Xs` (Array or Array of Arrays in case of MixedModel) and `events`-dataframe 
"""
struct DesignMatrix
    "Array of formulas"
    formulas::Any
    "A concatenated designmatric. In case of Mixed Models an array, where the first one is a FeMat, later ones ReMats. "
    Xs::Any
    "Event table with all events"
    events::Any
end

function DesignMatrix()
    return DesignMatrix([], [], [])
end


"""
Concrete type to implement an Mass-Univariate LinearModel.
`.design` contains the formula + times dict
`.designmatrix` contains a `DesignMatrix`
`modelfit` is a `Any` container for the model results
"""
mutable struct UnfoldLinearModel <: UnfoldModel
    design::Dict
    designmatrix::DesignMatrix
    modelfit::Any
end

UnfoldLinearModel(d::Dict) = UnfoldLinearModel(d, Unfold.DesignMatrix(), [])
UnfoldLinearModel(d::Dict, X::DesignMatrix) = UnfoldLinearModel(d, X, [])

"""
Concrete type to implement an Mass-Univariate LinearMixedModel.
`.design` contains the formula + times dict
`.designmatrix` contains a `DesignMatrix`
`modelfit` is a `Any` container for the model results
"""
mutable struct UnfoldLinearMixedModel <: UnfoldModel
    design::Dict
    designmatrix::DesignMatrix
    modelfit::Any#::Array{UnfoldMixedModelFitCollection} # optional info on the modelfit
end
UnfoldLinearMixedModel(d::Dict) = UnfoldLinearMixedModel(d, Unfold.DesignMatrix(), [])
UnfoldLinearMixedModel(d::Dict, X::DesignMatrix) = UnfoldLinearMixedModel(d, X, [])

"""
Concrete type to implement an deconvolution LinearModel.
`.design` contains the formula + times dict
`.designmatrix` contains a `DesignMatrix`
`modelfit` is a `Any` container for the model results
"""
mutable struct UnfoldLinearModelContinuousTime <: UnfoldModel
    design::Dict
    designmatrix::DesignMatrix
    modelfit::Any
end

UnfoldLinearModelContinuousTime(d::Dict) =
    UnfoldLinearModelContinuousTime(d, Unfold.DesignMatrix(), [])
UnfoldLinearModelContinuousTime(d::Dict, X::DesignMatrix) =
    UnfoldLinearModelContinuousTime(d, X, [])

"""
Concrete type to implement an deconvolution LinearMixedModel.

**Warning** This is to be treated with care, not much testing went into it.

`.design` contains the formula + times dict
`.designmatrix` contains a `DesignMatrix`
`modelfit` is a `Any` container for the model results
"""
mutable struct UnfoldLinearMixedModelContinuousTime <: UnfoldModel
    design::Dict
    designmatrix::DesignMatrix
    modelfit::Any#::UnfoldMixedModelFitCollection
end

UnfoldLinearMixedModelContinuousTime(d::Dict) =
    UnfoldLinearMixedModelContinuousTime(d, Unfold.DesignMatrix(), [])
UnfoldLinearMixedModelContinuousTime(d::Dict, X::DesignMatrix) =
    UnfoldLinearMixedModelContinuousTime(d, X, [])

"""
Contains the results of linearmodels (continuous and not)
"""
struct LinearModelFit <: ModelFit
    estimate::Any
    info::Any
    standarderror::Any
end

LinearModelFit() = LinearModelFit([], [], [])
LinearModelFit(estimate) = LinearModelFit(estimate, [], [])
LinearModelFit(estimate, info) = LinearModelFit(estimate, info, [])

function Base.show(io::IO, obj::UnfoldModel)
    println(io, "Unfold-Type: $(typeof(obj)) \n")
    println(io, "formula: $(obj.design)")
    println(
        io,
        "Useful functions:\n 
    design(uf) \t\t(returns Dict of event => (formula,times/basis))  \n
    designmatrix(uf) \t(returns DesignMatrix with events) \n
    modelfit(uf) \t\t(returns modelfit object) \n
    coeftable(uf) \t\t(returns tidy result dataframe) \n",
    )
end

"""
Abstract non-linear spline term.
Implemented in `UnfoldBSplineKitExt`
"""

abstract type AbstractSplineTerm <: AbstractTerm end
