struct UnfoldLinearModel
    beta::AbstractArray
    optim
    formula
    X::AbstractArray
end

struct UnfoldLinearMixedModel

end


struct UnfoldModel{M<:Union{AbstractArray{Union{LinearMixedModel,UnfoldLinearModel,UnfoldLinearMixedModel},1}, UnfoldLinearModel,LinearMixedModel,UnfoldLinearMixedModel}}
    model::M
    formula
    tbl::DataFrame
    results::DataFrame
end

function Base.show(io::IO, obj::UnfoldModel)
        println(io, "LinearModelTimeExpanded object")
        println(io, "formula: $(obj.formula)")

end

function Base.show(io::IO, obj::UnfoldModel)
    println(io, "UnfoldModel")
    # TODO Save the original formula without time expansion
    println(io, "Unique Terms: $(unique(obj.results.term))")
    println(io, "Basis Function Columns: $(obj.results.colnames_basis[1]) : $(obj.results.colnames_basis[end]))")
    println(io, "Fields: .model, .formula, .tbl, .results")
end
