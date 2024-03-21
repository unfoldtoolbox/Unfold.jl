function StatsModels.modelmatrix(model::UnfoldLinearMixedModel, bool)
    @assert bool == false "time continuous model matrix is not implemented for a `UnfoldLinearMixedModel`"
    return modelmatrix(model)
end
