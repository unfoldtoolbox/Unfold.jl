abstract type UnfoldModel end

struct UnfoldLinearModel <:UnfoldModel
    beta::AbstractArray # chan x predictor or chan x time x predictor
    modelinfo # optional info on the modelfit
    X::DesignMatrix
end

struct UnfoldLinearMixedModel <:UnfoldModel
    beta::AbstractArray # chan x predictor or chan x time x predictor
    sigma::AbstractArray # chan x ranef or chan x time x ranef
    modelinfo # optional info on the modelfit
    X::DesignMatrix
end


function Base.show(io::IO, obj::UnfoldModel)
        println(io, "LinearModelTimeExpanded object")
        println(io, "formula: $(obj.X.formulas)")
        println(io, "Fields: .model ('UnfoldModel' with felds .modelinfo, .beta, .formula, .Xs) \n .results (tidy result table)")

end
