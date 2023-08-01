
struct UnfoldMixedModelFitCollection{T<:AbstractFloat} <:
    MixedModels.MixedModelFitCollection{T}
 fits::Vector
 λ::Vector{<:Union{LowerTriangular{T,Matrix{T}},Diagonal{T,Vector{T}}}}
 inds::Vector{Vector{Int}}
 lowerbd::Vector{T}
 fcnames::NamedTuple
end