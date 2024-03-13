
struct DesignMatrixLinearMixedModel{T} <: AbstractDesignMatrix{T}
    formulas::Vector{FormulaTerm} # "Array of formulas"
    Xs::Vector{Vector{Union{FeMat{T},ReMat{T}}}} #"A concatenated designmatric. In case of Mixed Models an array, where the first one is a FeMat, later ones ReMats. "
    events::DataFrame #"Event table with all events"
end

struct DesignMatrixLinearMixedModelContinuousTime{T} <: AbstractDesignMatrix{T}
    formulas::Vector{FormulaTerm} # "Array of formulas"
    Xs::Vector{Union{FeMat{T},ReMat{T}}} #"A concatenated designmatric. In case of Mixed Models an array, where the first one is a FeMat, later ones ReMats. "
    events::DataFrame #"Event table with all events"
end



struct intern_LinearMixedModelFitCollection{T} <: MixedModels.MixedModelFitCollection{T}
    fits::Vector
    λ::Vector
    inds::Vector{Vector{Int}}
    lowerbd::Vector{T}
    fcnames::NamedTuple
end



struct UnfoldLinearMixedModelFitCollection{T<:AbstractFloat,N} <: AbstractModelFit{T,N}
    field::intern_LinearMixedModelFitCollection
end





"""
Concrete type to implement an Mass-Univariate LinearMixedModel.
`.design` contains the formula + times dict
`.designmatrix` contains a `DesignMatrix`
`modelfit` is a `Any` container for the model results
"""
mutable struct UnfoldLinearMixedModel{T} <: UnfoldModel{T}
    design::Dict
    designmatrix::DesignMatrixLinearMixedModel{T}
    modelfit::Vector{UnfoldLinearMixedModelFitCollection{T,3}} # optional info on the modelfit
end
UnfoldLinearMixedModel(d::Dict) =
    UnfoldLinearMixedModel(d, Unfold.DesignMatrixLinearMixedModel(), [])
UnfoldLinearMixedModel(d::Dict, X::AbstractDesignMatrix) = UnfoldLinearMixedModel(d, X, [])


"""
Concrete type to implement an deconvolution LinearMixedModel.

**Warning** This is to be treated with care, not much testing went into it.

`.design` contains the formula + times dict
`.designmatrix` contains a `DesignMatrix`
`modelfit` is a `Any` container for the model results
"""
mutable struct UnfoldLinearMixedModelContinuousTime{T} <: UnfoldModel{T}
    design::Dict
    designmatrix::DesignMatrixLinearMixedModelContinuousTime{T}
    modelfit::UnfoldLinearMixedModelFitCollection{T,2}
end

UnfoldLinearMixedModelContinuousTime(d::Dict) = UnfoldLinearMixedModelContinuousTime(
    d,
    Unfold.DesignMatrixLinearMixedModelContinuousTime(),
    [],
)
UnfoldLinearMixedModelContinuousTime(d::Dict, X::AbstractDesignMatrix) =
    UnfoldLinearMixedModelContinuousTime(d, X, [])


@traitimpl ContinuousTimeTrait{UnfoldLinearMixedModelContinuousTime}
