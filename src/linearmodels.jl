struct UnfoldModel{M<:Union{AbstractArray{LinearMixedModel},LinearMixedModel}}
    model::M
    formula::FormulaTerm
    tbl::DataFrame
    results::DataFrame
end

function Base.show(io::IO, obj::UnfoldModel)
        println(io, "LinearModelTimeExpanded object")
        println(io, "formula: $(obj.formula)")
end

struct UnfoldLinearModel

end

struct UnfoldLinearMixedModel

end
