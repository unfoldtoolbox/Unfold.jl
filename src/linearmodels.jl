struct UnfoldLinearMixedModel{T<:AbstractMatrix,T2<:AbstractMatrix}
    formula::FormulaTerm
    tbl::DataFrame
    lmm::MixedModel
end

struct UnfoldLinearModel{T<:AbstractMatrix,T2<:AbstractMatrix}
    formula::FormulaTerm
    tbl::DataFrame
    X::T
    beta::T2
    optim::Any
end


function Base.show(io::IO, obj::UnfoldLinearModel)
        println(io, "LinearModelTimeExpanded object")
        println(io, "formula: $(obj.formula)")
end
