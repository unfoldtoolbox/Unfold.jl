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
function FileIO.save(file, uf::T; compress = false) where {T<:UnfoldModel}

    mm = modelmatrix(uf)
    #nd = length(size(mm[1]))
    d = designmatrix(uf)
    jldopen(file, "w"; compress = compress) do f
        f["uf"] = T(
            design(uf),
            [
                typeof(uf.designmatrix[k])(
                    Unfold.formulas(uf)[k],
                    empty_modelmatrix(designmatrix(uf)[k]),
                    Unfold.events(uf)[k],
                ) for k = 1:length(uf.designmatrix)
            ],
            uf.modelfit,
        )
    end
end

function empty_modelmatrix(d::AbstractDesignMatrix)

    return typeof(d)().modelmatrix
end
function FileIO.save(
    file,
    uf::T;
    compress = false,
) where {T<:UnfoldLinearModelContinuousTime}

    mm = modelmatrix(uf)
    nd = length(size(mm))
    jldopen(file, "w"; compress = compress) do f
        f["uf"] = T(
            design(uf),
            [
                typeof(uf.designmatrix[k])(
                    Unfold.formulas(uf)[k],
                    sparse(Array{typeof(mm[k])}(undef, repeat([0], nd)...)),
                    Unfold.events(uf)[k],
                ) for k = 1:length(uf.designmatrix)
            ],
            modelfit(uf),
        )
    end
end

"""
    FileIO.load(file, ::Type{<:UnfoldModel}; generate_Xs=true)
    
Load UnfoldModel from a .jld2 file. 

By default, the designmatrix is reconstructed. If it is not needed set `generate_Xs=false`
which improves time-efficiency.
"""
function FileIO.load(file, ::Type{<:UnfoldModel}; generate_Xs = true)
    f = jldopen(file, "r")
    uf = f["uf"]
    close(f)

    form = formulas(designmatrix(uf))
    events = Unfold.events(designmatrix(uf))

    @debug typeof(form), typeof(events)
    # potentially don't generate Xs, but always generate it for LinearModels as it is small + cheap + we actually need it for many functions
    if generate_Xs || uf isa UnfoldLinearModel
        X = Unfold.modelcols(form, events)
    else
        #nd = length(size(modelmatrix(designmatrix(uf)[1])))
        #@debug eltype(uf), nd, length(form)
        #X = [sparse(Array{eltype(uf)}(undef, repeat([0], nd)...)) for k = 1:length(form)]
        X = [empty_modelmatrix(designmatrix(uf)[k]) for k = 1:length(form)]

    end

    modelfit =
        if isa(uf.modelfit, JLD2.ReconstructedMutable{Symbol("Unfold.LinearModelFit")})
            @warn "old Unfold Model detected, trying to 'upgrade' uf.modelfit"
            mf = uf.modelfit
            T = typeof(mf.estimate)
            if isempty(mf.standarderror)
                LinearModelFit(mf.estimate, mf.info)
            else
                LinearModelFit(mf.estimate, mf.info, T(mf.standarderror))
            end
        else
            uf.modelfit
        end

    @debug typeof(form), typeof(X), typeof(events), size(form), size(X), size(events)
    # reintegrate the designmatrix
    return typeof(uf)(uf.design, typeof(designmatrix(uf)[1]).(form, X, events), modelfit)
end
