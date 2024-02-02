
function Unfold.TimeExpandedTerm(
    term::NTuple{N,Union{<:AbstractTerm,<:MixedModels.RandomEffectsTerm}},
    basisfunction,
    eventfields::Array{Symbol,1},
) where {N}
    # for mixed models, apply it to each Term
    Unfold.TimeExpandedTerm.(term, Ref(basisfunction), Ref(eventfields))
end
