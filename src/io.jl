using FileIO
using JLD2
using StatsModels
using Unfold

"""
    StatsModels.modelcols(forms::Vector,events::Vector)
"""
function StatsModels.modelcols(forms::Vector, events::Vector)
    @assert length(forms) == length(events)
    return StatsModels.modelcols.(forms, events)
end


"""
    StatsModels.modelcols(form::FormulaTerm, events)
"""
StatsModels.modelcols(form::FormulaTerm, events) = modelcols(form.rhs, events)


"""
    FileIO.save(file, uf::T; compress=false) where {T<:UnfoldModel}

Save UnfoldModel in a (by default uncompressed) .jld2 file.

For memory efficiency the designmatrix is set to missing.
If needed, it can be reconstructed when loading the model.
"""
function FileIO.save(file, uf::T; compress=false) where {T<:UnfoldModel}
    jldopen(file, "w"; compress=compress) do f
        f["uf"] = T(uf.design, Unfold.DesignMatrix(designmatrix(uf).formulas, missing, designmatrix(uf).events), uf.modelfit)
    end
end


"""
    FileIO.load(file, ::Type{<:UnfoldModel}; generate_Xs=true)
    
Load UnfoldModel from a .jld2 file. 

By default, the designmatrix is reconstructed. If it is not needed set `generate_Xs=false`
which improves time-efficiency.
"""
function FileIO.load(file, ::Type{<:UnfoldModel}; generate_Xs=true)
    f = jldopen(file, "r")
    uf = f["uf"]
    close(f)

    form = designmatrix(uf).formulas
    events = designmatrix(uf).events
    # potentially don't generate Xs, but always generate it for LinearModels as it is small + cheap + we actually need it for many functions
    if generate_Xs || uf isa UnfoldLinearModel
        X = Unfold.modelcols(form, events)
    else
        X = missing
    end


    # reintegrate the designmatrix
    return typeof(uf)(uf.design, Unfold.DesignMatrix(form, X, events), uf.modelfit)
end