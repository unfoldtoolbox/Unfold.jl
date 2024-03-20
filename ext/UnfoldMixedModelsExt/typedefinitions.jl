
struct DesignMatrixLinearMixedModel{T} <: AbstractDesignMatrix{T}
    formula::FormulaTerm
    modelmatrix::Tuple
    events::DataFrame
end

struct DesignMatrixLinearMixedModelContinuousTime{T} <: AbstractDesignMatrix{T}
    formula::FormulaTerm
    modelmatrix::Tuple
    events::DataFrame
end

DesignMatrixLinearMixedModel{T}() where {T} =
    DesignMatrixLinearMixedModel{T}(FormulaTerm(:empty, :empty), (), DataFrame())
DesignMatrixLinearMixedModelContinuousTime{T}() where {T} =
    DesignMatrixLinearMixedModelContinuousTime{T}(
        FormulaTerm(:empty, :empty),
        (),
        DataFrame(),
    )

struct LinearMixedModelFitCollection{T} <: MixedModels.MixedModelFitCollection{T}
    fits::Vector
    λ::Vector
    inds::Vector{Vector{Int}}
    lowerbd::Vector{T}
    fcnames::NamedTuple
end

LinearMixedModelFitCollection{T}() where {T} =
    LinearMixedModelFitCollection{T}([], [], Vector{Int}[], T[], (;))

struct UnfoldLinearMixedModelFit{T<:AbstractFloat,N} <: AbstractModelFit{T,N}
    collection::LinearMixedModelFitCollection{T}
end


UnfoldLinearMixedModelFit{T,N}() where {T,N} =
    UnfoldLinearMixedModelFit{T,N}(LinearMixedModelFitCollection{T}())

"""
Concrete type to implement an Mass-Univariate LinearMixedModel.
`.design` contains the formula + times dict
`.designmatrix` contains a `DesignMatrix`
`modelfit` is a `Any` container for the model results
"""
mutable struct UnfoldLinearMixedModel{T} <: UnfoldModel{T}
    design::Vector{<:Pair}
    designmatrix::Vector{<:DesignMatrixLinearMixedModel{T}}
    modelfit::UnfoldLinearMixedModelFit{T,3} # optional info on the modelfit
end

UnfoldLinearMixedModel(args...; kwargs...) =
    UnfoldLinearMixedModel{Float64}(args...; kwargs...)
UnfoldLinearMixedModel{T}() where {T} = UnfoldLinearMixedModel{T}(Pair[])
UnfoldLinearMixedModel{T}(d::Vector) where {T} = UnfoldLinearMixedModel{T}(
    d,
    [DesignMatrixLinearMixedModel{T}()],
    UnfoldLinearMixedModelFit{T,3}(),
)
UnfoldLinearMixedModel{T}(d::Vector, X::Vector{<:AbstractDesignMatrix}) where {T} =
    UnfoldLinearMixedModel{T}(d, X, DataFrame())


"""
Concrete type to implement an deconvolution LinearMixedModel.

**Warning** This is to be treated with care, not much testing went into it.

`.design` contains the formula + times dict
`.designmatrix` contains a `DesignMatrix`
`.modelfit` is a `Any` container for the model results
"""
mutable struct UnfoldLinearMixedModelContinuousTime{T} <: UnfoldModel{T}
    design::Vector{<:Pair}
    designmatrix::Vector{<:DesignMatrixLinearMixedModelContinuousTime{T}}
    modelfit::UnfoldLinearMixedModelFit{T,2}
end

UnfoldLinearMixedModelContinuousTime(args...; kwargs...) =
    UnfoldLinearMixedModelContinuousTime{Float64}(args...; kwargs...)
UnfoldLinearMixedModelContinuousTime{T}() where {T} =
    UnfoldLinearMixedModelContinuousTime{T}(Pair[])

UnfoldLinearMixedModelContinuousTime{T}(d::Vector) where {T} =
    UnfoldLinearMixedModelContinuousTime{T}(
        d,
        [DesignMatrixLinearMixedModelContinuousTime{T}()],
        UnfoldLinearMixedModelFit{T,2}(),
    )
UnfoldLinearMixedModelContinuousTime{T}(
    d::Vector{<:Pair},
    X::Vector{<:AbstractDesignMatrix{T}},
) where {T} =
    UnfoldLinearMixedModelContinuousTime{T}(d, X, UnfoldLinearMixedModelFit{T,2}())


@traitimpl Unfold.ContinuousTimeTrait{UnfoldLinearMixedModelContinuousTime}
