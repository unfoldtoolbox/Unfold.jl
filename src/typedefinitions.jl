
"""
using Base: @deprecate_binding
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
Type that keeps an Array of  `formulas`, designmatrices `modelmatrix` (Array or Array of Arrays in case of MixedModel) and `events`-dataframe 
"""
struct DesignMatrixLinearModel{T} <: AbstractDesignMatrix{T}
    formula::FormulaTerm # "Array of formulas"
    modelmatrix::Array{T,2} #"A concatenated designmatric. In case of Mixed Models an array, where the first one is a FeMat, later ones ReMats. "
    events::DataFrame #"Event table with all events"
end

struct DesignMatrixLinearModelContinuousTime{T} <: AbstractDesignMatrix{T}
    formula::FormulaTerm # "Array of formulas"
    modelmatrix::SparseMatrixCSC{T} #"A concatenated designmatric. In case of Mixed Models an array, where the first one is a FeMat, later ones ReMats. "
    events::DataFrame #"Event table with all events"
end


"""
in case a single FormulaTerm is inputted, this is wrapped into a vector
"""
# DesignMatrixLinearModel(f::FormulaTerm, modelmatrix, events) =
#     DesignMatrixLinearModel([f], modelmatrix, events)
# DesignMatrixLinearModel(f, modelmatrix, events::DataFrame) = DesignMatrixLinearModel(f, modelmatrix, [events])

DesignMatrixLinearModel(args...) = DesignMatrixLinearModel{Float64}()
DesignMatrixLinearModel{T}() where {T} = DesignMatrixLinearModel{T}(
    FormulaTerm(:empty, :empty),
    Array{T,2}(undef, 0, 0),
    DataFrame(),
)
DesignMatrixLinearModelContinuousTime(args...) =
    DesignMatrixLinearModelContinuousTime{Float64}()
DesignMatrixLinearModelContinuousTime{T}() where {T} =
    DesignMatrixLinearModelContinuousTime{T}(
        FormulaTerm(:empty, :empty),
        SparseMatrixCSC(Array{T,2}(undef, 0, 0)),
        DataFrame(),
    )
# DesignMatrixLinearModelContinuousTime(f, modelmatrix, events::DataFrame) =
#     DesignMatrixLinearModelContinuousTime(f, modelmatrix, [events])
# DesignMatrixLinearModelContinuousTime(f, modelmatrix::SparseMatrixCSC, events) =
#     DesignMatrixLinearModelContinuousTime(f, [modelmatrix], events)

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
    design::Vector{<:Pair{<:Any,<:Tuple}}
    designmatrix::Vector{<:DesignMatrixLinearModel{T}}
    modelfit::LinearModelFit{T,3}
end

# empty model
UnfoldLinearModel(args...) = UnfoldLinearModel{Float64}(args...)
UnfoldLinearModel{T}() where {T} = UnfoldLinearModel{T}([:empty => ()])
# backward compatible with dict
UnfoldLinearModel{T}(d::Dict, args...) where {T} =
    UnfoldLinearModel(collect(pairs(d)), args...)
# only design specified
UnfoldLinearModel{T}(d::Vector{<:Pair}) where {T} =
    UnfoldLinearModel{T}(d, [DesignMatrixLinearModel{T}()])
# matrix not a vector of matrix
UnfoldLinearModel{T}(d::Vector{<:Pair}, X::AbstractDesignMatrix) where {T} =
    UnfoldLinearModel{T}(d, [X])
# no modelfit
UnfoldLinearModel{T}(d::Vector{<:Pair}, X::Vector{<:AbstractDesignMatrix{T}}) where {T} =
    UnfoldLinearModel{T}(deepcopy(d), X, LinearModelFit(Array{T,3}(undef, 0, 0, 0)))



"""
Concrete type to implement an deconvolution LinearModel.
`.design` contains the formula + times dict
`.designmatrix` contains a `DesignMatrix`
`modelfit` is a `Any` container for the model results
"""
mutable struct UnfoldLinearModelContinuousTime{T} <: UnfoldModel{T}
    design::Vector{<:Pair{<:Any,<:Tuple}}
    designmatrix::Vector{<:DesignMatrixLinearModelContinuousTime{T}}
    modelfit::LinearModelFit{T,2}
end


# empty model
UnfoldLinearModelContinuousTime(args...) = UnfoldLinearModelContinuousTime{Float64}(args...)
UnfoldLinearModelContinuousTime{T}() where {T} =
    UnfoldLinearModelContinuousTime{T}([:empty => ()])
# backward compatible with dict
UnfoldLinearModelContinuousTime{T}(d::Dict, args...) where {T} =
    UnfoldLinearModelContinuousTime{T}(keys(d) .=> values(d), args...)

# only design specified
UnfoldLinearModelContinuousTime{T}(d::Vector{<:Pair}) where {T} =
    UnfoldLinearModelContinuousTime{T}(d, Unfold.DesignMatrixLinearModelContinuousTime{T}())
# matrix not a vector of matrix
UnfoldLinearModelContinuousTime{T}(d::Vector{<:Pair}, X::AbstractDesignMatrix) where {T} =
    UnfoldLinearModelContinuousTime{T}(d, [X])

# no modelfit
UnfoldLinearModelContinuousTime{T}(
    d::Vector{<:Pair},
    X::Vector{<:AbstractDesignMatrix{T}},
) where {T} = UnfoldLinearModelContinuousTime{T}(
    deepcopy(d),
    X,
    LinearModelFit(Array{T,2}(undef, 0, 0)),
)

#----
# Traits definitions
@traitdef ContinuousTimeTrait{X}

@traitimpl ContinuousTimeTrait{UnfoldLinearModelContinuousTime}
#---



"""
Abstract non-linear spline term.
Implemented in `UnfoldBSplineKitExt`
"""

abstract type AbstractSplineTerm <: AbstractTerm end
