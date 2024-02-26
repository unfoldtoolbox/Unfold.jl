
struct UnfoldMixedModelFitCollection{T<:AbstractFloat} <:
       MixedModels.MixedModelFitCollection{T}
    fits::Vector
    λ::Vector
    inds::Vector{Vector{Int}}
    lowerbd::Vector{T}
    fcnames::NamedTuple
end
