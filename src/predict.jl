

using DocStringExtensions: format



# opt-in

# Most important function ;-)
#predict(X::AbstractMatrix, coefs::AbstractMatrix) = coefs * X'
predict(X::AbstractMatrix, coefs::AbstractMatrix) = (X * coefs')'

# 3D-case
function predict(
    X::AbstractMatrix,
    coef::AbstractArray{T,3},
) where {T<:Union{Missing,<:Number}}
    yhat = Array{T}(undef, size(coef, 1), size(X, 1), size(coef, 2))
    for ch = 1:size(coef, 1)
        yhat[ch, :, :] .= X * permutedims(coef[ch, :, :], (2, 1))
    end
    return permutedims(yhat, (1, 3, 2)) # esults in: ch x time x pred
end

function predict(
    X::AbstractMatrix,
    coefs::AbstractMatrix{T},
    onsets::Vector,
    timewindow::NTuple{2,Int},
) where {T}
    # Create an empty array for the residuals
    # Note there is a +1 because start and stop indices are both included in the interval
    yhat = missings(T, size(coefs, 1), timewindow[2] - timewindow[1] + 1, length(onsets))
    for event in eachindex(onsets)

        event_range_temp = onsets[event]+timewindow[1]:onsets[event]+timewindow[2]

        # Clip indices that are outside of the design matrix (i.e. before the start or after the end)
        #indices_inside_data = 0 .< event_range_temp .< size(data, 2)
        indices_inside_data = 0 .< event_range_temp .< size(X, 1)
        event_range = event_range_temp[indices_inside_data]

        # Calculate residuals
        yhat[:, indices_inside_data, event] .= predict(@view(X[event_range, :]), coefs)
    end
    return yhat
end


residuals(uf, data::AbstractArray) = data .- predict(uf)

predict(uf::UnfoldModel, f::FormulaTerm, args...; kwargs...) =
    predict(uf, [f], args...; kwargs...)
predict(uf::UnfoldModel, f::Vector, evts::DataFrame; kwargs...) =
    predict(uf, f, repeat([evts], length(f)); kwargs...)
predict(uf::UnfoldModel; kwargs...) = predict(uf, formulas(uf), events(uf); kwargs...)
predict(uf::UnfoldModel, f::Vector{<:FormulaTerm}; kwargs...) =
    predict(uf, f, events(uf); kwargs...)
predict(uf::UnfoldModel, evts::Vector{<:DataFrame}; overlap = false, kwargs...) =
    predict(uf, Unfold.formulas(uf), evts; overlap, kwargs...)

"""
    function predict(
        uf,
        f::Vector{<:FormulaTerm},
        evts::Vector{<:DataFrame};
        overlap = true,
        kwargs...
    )
Returns a predicted ("y_hat = X*b") `Array`. 

- `uf` is an `<:UnfoldModel`
- `f` is a (vector of) formulas, typically `Unfold.formulas(uf)`, but formulas can be modified e.g. by `effects`.
- `evts` is a (vector of) events, can be `Unfold.events(uf)` to return the (possibly continuous-time) predictions of the model. Can be a custom even


## kwargs:
if `overlap = true` (default), overlap based on the `latency` column of 'evts` will be simulated, or in the case of `!ContinuousTimeTrait` just X*coef is returned. 
if `overlap = false`, predict coefficients without overlap (models with `ContinuousTimeTrait` (=> with basisfunction / deconvolution) only), via `predict_no_overlap`

if `keep_basis` or `exclude_basis` is defined, then `predict_partial_overkap` is called, which allows to selective introduce overlap based on specified (or excluded respective) events/basisfunctions

`epoch_to` and  `epoch_timewindow` currently only defined for partial_overlap, calculated (partial) overlap controlled predictions, but returns them at the specified `epoch_at` event, with the times `epoch_timewindow` in samples

Hint: all vectors can be "single" types, and will be containered in a vector

## Output

- If `overlap=false`, returns a 3D-Array
- If `overlap=true` and `epoch_to=nothing` (default), returns a 2D-array
- If `overlap=true` and `epoch_to!=nothing`, returns a 3D array
"""
function predict(
    uf,
    f::Vector{<:FormulaTerm},
    evts::Vector{<:DataFrame};
    overlap = true,
    keep_basis = [],
    exclude_basis = [],
    epoch_to = nothing,
    epoch_timewindow = nothing,
    #eventcolumn = :event,
)
    @assert !(!isempty(keep_basis) & !isempty(exclude_basis)) "choose either to keep events, or to exclude, but not both"


    coefs = coef(uf)

    if overlap
        if isempty(keep_basis) & isempty(exclude_basis)
            @debug "full-overlap"
            if events(uf) == evts
                @debug "original design predict"
                # off-the-shelf standard continuous predict, as you'd expect
                return predict(modelmatrix(uf), coef(uf))
            else
                @debug "new dataframe predict"
                # predict of new data-table with a latency column. First generate new designmatrix, then off-the-shelf X*b predict
                @assert length(f) == length(evts) "please provide the same amount of formulas (currently $(length(f)) as event-dataframes (currently $(length(evts)))"
                X_new = modelcols.(f, evts) |> Unfold.extend_to_larger
                return predict(X_new, coefs)
            end
        else
            return predict_partial_overlap(
                uf,
                coefs,
                evts;
                keep_basis,
                exclude_basis,
                epoch_to,
                epoch_timewindow,
            )
        end
    else

        @debug "no overlap predict2"
        # no overlap, the "ideal response". We predict each event row independently. user could concatenate if they really want to :)
        return predict_no_overlap(uf, coefs, f, evts)
    end
end

@traitfn function predict_partial_overlap(
    uf::T,
    coefs,
    evts;
    keep_basis = [],
    exclude_basis = [],
    epoch_to = nothing,
    epoch_timewindow = nothing,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}}
    @assert !isempty(keep_basis) & !isempty(exclude_basis) "can't have no overlap & specify keep/exclude at the same time. decide for either case"
    # Partial overlap! we reconstruct with some basisfunctions deactivated
    if !isempty(keep_basis)
        basisnames = keep_basis
    else
        basisnames = basisname(uf)
        basisnames = setdiff(basisnames, exclude_basis)
    end

    X_view = matrix_by_basisname(modelmatrix(uf), uf, basisname)
    coefs_view = matrix_by_basisname(coefs, uf, basisname)


    if isnothing(epoch_to)
        return predict(X_view, coefs_view)

    else
        timewindow =
            isnothing(epoch_timewindow) ? calc_epoch_timewindow(uf, epoch_to) :
            epoch_timewindow
        ix_event = basisname .== epoch_to
        latencies = evts[ix_event].latency
        return predict(X_view, coefs_view, latencies, (timewindow[1], timewindow[end]);)
    end

end

@traitfn function predict_no_overlap(
    uf::T,
    coefs,
    f::Vector,
    evts::Vector,
) where {T <: UnfoldModel; !ContinuousTimeTrait{T}}
    @debug "Not ContinuousTime yhat, Array"
    X = modelcols.(f, evts)

    # figure out which coefficients belong to which event
    Xsizes = size.(X, Ref(2))
    Xsizes_cumsum = vcat(0, cumsum(Xsizes))
    indexes = [(Xsizes_cumsum[ix]+1):Xsizes_cumsum[ix+1] for ix = 1:length(f)]

    # split up the coefs accordingly
    coArray = [@view coefs[:, :, ix] for ix in indexes]

    return predict.(X, coArray)
end

@traitfn function predict_no_overlap(
    uf::T,
    coefs,
    f::Vector,
    evts::Vector,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}}

    yhat = Array{eltype(coefs)}[]
    for (fi, e) in zip(f, evts)

        e.latency .= sum(times(fi) .<= 0)
        X_singles = map(x -> modelcols(fi, DataFrame(x)), eachrow(e))
        #=
        if typeof(fi.rhs.basisfunction) <: FIRBasis
            # this pertains only to FIR-models
            # remove the last time-point because it is attached due to non-integer latency/eventonsets.
            # e.g. x denotes a sample
            # x- - -x- - -x- - -x
            # e- - - -f- - - -g-
            # 
            # e is aligned (integer) with the sample
            # f&g are between two samples, thus the design matrix would interpolate between them. Thus has as a result, that the designmatrix is +1 longer than what would naively be expected
            #
            # because in "predict" we define where samples onset, we can remove the last sample, it s always 0 anyway, but to be sure we test it

            X_singles = map(x -> x[1:end-1, :], X_singles)
        end
        =#
        coefs_view = matrix_by_basisname(coefs, (uf), (basisname([fi])))

        yhat_single = similar(
            coefs_view,
            size(coefs_view, 1),
            size(X_singles[1], 1),
            length(X_singles),
        )

        for ev = 1:length(X_singles)
            p = predict(X_singles[ev], coefs_view)
            yhat_single[:, :, ev] .= p
        end
        #@debug yhat_single
        push!(yhat, yhat_single)

    end
    return yhat
end

"""
Returns a view of the Matrix `y`, according to the indices of the timeexpanded `basisname`
"""
function matrix_by_basisname(y::AbstractMatrix, uf, basisnames::Vector)
    ix = get_basis_indices(uf, basisnames)
    return @view(y[:, ix])
end

"""
returns an integer range with the samples around `epoch_event` as defined in the corresponding basisfunction

"""

function calc_epoch_timewindow(uf, epoch_event)
    basis_ix = findfirst(Unfold.basisname(formulas(uf)) .== epoch_event)
    basisfunction = formulas(uf)[basis_ix].rhs.basisfunction
    return (1:Unfold.height(basisfunction)) .+ Unfold.shift_onset(basisfunction)
end

"""
    get_basis_indices(uf, basisnames::Vector)
returns a boolean vector with length spanning all coefficients, which coefficient is defined by `basisnames` (vector of names)
"""
get_basis_indices(uf, basisnames::Vector) =
    reduce(vcat, Unfold.get_basis_names(uf)) .∈ Ref(basisnames)
get_basis_indices(uf, basisname) = get_basis_indices(uf, [basisname])

"""
    predict_table(model<:UnfoldModel,events=Unfold.events(model),args...;kwargs...)
Shortcut to call efficiently call (pseudocode) `result_to_table(predict(...))`.

Returns a tidy DataFrame with the predicted results. Loops all input to `predict`, but really only makes sense to use if you specify either:

`overlap = false` or `epoch_to = "eventname"`.

"""
function predict_table(
    model::UnfoldModel,
    events = Unfold.events(model),
    args...;
    kwargs...,
)
    p = predict(model, Unfold.formulas(model), events, args...; kwargs...)
    return result_to_table(model, p, events)
end

eventnames(model::UnfoldModel) = first.(design(model))




"""
    times(model<:UnfoldModel)
returns arrays of time-vectors, one for each basisfunction / parallel-fitted-model (MassUnivarite case)
"""
@traitfn times(model::T) where {T <: UnfoldModel; !ContinuousTimeTrait{T}} =
    times(design(model))

@traitfn times(model::T) where {T <: UnfoldModel; ContinuousTimeTrait{T}} =
    times(formulas(model))

times(d::Vector) = times.(d)
times(d::FormulaTerm) = times(d.rhs)
times(d::MatrixTerm) = times(d.term)
times(d::TimeExpandedTerm) = times(d.basisfunction)
function times(
    d::Vector{
        <:Pair{
            <:Union{<:DataType,<:AbstractString,<:Symbol},
            <:Tuple{<:AbstractTerm,<:AbstractVector},
        },
    },
)
    all_times = last.(last.(d)) #[k[2] for k in values(d)] # probably going for steprange would be better
    @debug all_times
    @assert all(all_times .== all_times[1:1]) "all times need to be equal in a mass univariate model"
    return all_times
end

