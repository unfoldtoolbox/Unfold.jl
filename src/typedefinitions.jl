
"""
The main abstract model-type of the toolbox. E.g. `UnfoldLinearModel` is a concrete type of this
"""
abstract type UnfoldModel{T} end


"""
Abstract Type to report modelresults
"""
abstract type AbstractModelFit{T,N} end


abstract type AbstractDesignMatrix{T} end


"""
    DesignMatrix
Type that keeps an Array of  `formulas`, designmatrices `Xs` (Array or Array of Arrays in case of MixedModel) and `events`-dataframe 
"""
struct DesignMatrixLinearModel{T} <: AbstractDesignMatrix{T}
    formulas::Vector{FormulaTerm} # "Array of formulas"
    Xs::Vector{Array{T}} #"A concatenated designmatric. In case of Mixed Models an array, where the first one is a FeMat, later ones ReMats. "
    events::DataFrame #"Event table with all events"
end

struct DesignMatrixLinearModelContinuousTime{T} <: AbstractDesignMatrix{T}
    formulas::Vector{FormulaTerm} # "Array of formulas"
    Xs::SparseMatrixCSC{T} #"A concatenated designmatric. In case of Mixed Models an array, where the first one is a FeMat, later ones ReMats. "
    events::DataFrame #"Event table with all events"
end


function DesignMatrix()
    return DesignMatrix([], [], [])
end


"""
Contains the results of linearmodels (continuous and not)
"""
struct LinearModelFit{T,N} <: AbstractModelFit{T,N}
    estimate::Array{T,N}
    info::Any
    standarderror::Array{T,N}
end

LinearModelFit() = LinearModelFit(Float64[], [], Float64[])
LinearModelFit(estimate) = LinearModelFit(estimate, [])
LinearModelFit(estimate::Array{T,2}, info) where {T} =
    LinearModelFit(estimate, info, similar(Array{T}, 0, 0))
LinearModelFit(estimate::Array{T,3}, info) where {T} =
    LinearModelFit(estimate, info, similar(Array{T}, 0, 0, 0))



"""
Concrete type to implement an Mass-Univariate LinearModel.
`.design` contains the formula + times dict
`.designmatrix` contains a `DesignMatrix`
`modelfit` is a `Any` container for the model results
"""
mutable struct UnfoldLinearModel{T} <: UnfoldModel{T}
    design::Dict
    designmatrix::DesignMatrixLinearModel{T}
    modelfit::LinearModelFit{T,3}
end

UnfoldLinearModel(d::Dict) = UnfoldLinearModel(d, Unfold.DesignMatrixLinearModel(), [])
UnfoldLinearModel(d::Dict, X::AbstractDesignMatrix) = UnfoldLinearModel(d, X, [])


"""
Concrete type to implement an deconvolution LinearModel.
`.design` contains the formula + times dict
`.designmatrix` contains a `DesignMatrix`
`modelfit` is a `Any` container for the model results
"""
mutable struct UnfoldLinearModelContinuousTime{T} <: UnfoldModel{T}
    design::Dict
    designmatrix::DesignMatrixLinearModelContinuousTime{T}
    modelfit::LinearModelFit{T,2}
end

UnfoldLinearModelContinuousTime(d::Dict) =
    UnfoldLinearModelContinuousTime(d, Unfold.DesignMatrixLinearModelContinuousTime(), [])
UnfoldLinearModelContinuousTime(d::Dict, X::AbstractDesignMatrix) =
    UnfoldLinearModelContinuousTime(d, X, [])

#----
# Traits definitions
@traitdef ContinuousTimeTrait{X}

@traitimpl ContinuousTimeTrait{UnfoldLinearModelContinuousTime}
#---


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
