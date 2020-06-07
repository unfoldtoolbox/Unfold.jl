abstract type UnfoldModel end

struct UnfoldLinearModel
    beta::AbstractArray # chan x predictor or chan x time x predictor
    modelinfo # optional info on the modelfit
    X::DesignMatrix
end

struct UnfoldLinearMixedModel
    beta::AbstractArray # chan x predictor or chan x time x predictor
    sigma::AbstractArray # chan x ranef or chan x time x ranef
    modelinfo # optional info on the modelfit
    X::DesignMatrix
end

struct UnfoldResult
    model
    results
end


function Base.show(io::IO, obj::UnfoldModel)
        println(io, "LinearModelTimeExpanded object")
        println(io, "formula: $(obj.Xs.formula)")

end

function Base.show(io::IO, obj::UnfoldResult)
    println(io, "UnfoldModel")
    # TODO Save the original formula without time expansion
    println(io, "Unique Terms: $(unique(obj.results.term))")
    println(io, "Basis Function Columns: $(obj.results.colnames_basis[1]) : $(obj.results.colnames_basis[end]))")
    println(io, "Fields: .model ('UnfoldModel' with felds .modelinfo, .beta, .formula, .Xs) \n .results (tidy result table)")
end
