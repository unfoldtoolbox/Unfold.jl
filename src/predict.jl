

#using DocStringExtensions: format



# opt-in

# Most important function ;-)
#predict(X::AbstractMatrix, coefs::AbstractMatrix) = coefs * X'
predict(X::AbstractMatrix, coefs::AbstractMatrix) = (X * coefs')'

# MassUnivariate Model, we have to extract the coefficients per designmatrix
function predict(X::Vector{<:AbstractMatrix}, coefs::AbstractArray{T,3}) where {T}
    len = size.(X, 2)

    @assert length(unique(len)) == 1 "all designmatrices need to be the same size in this case"
    coefs_views = [@view(coefs[:, :, (i-1)*len[1]+1:i*len[1]]) for i = 1:length(len)]
    predict.(X, coefs_views)
end

# 3D-case
function predict(
    X::AbstractMatrix{Tx},
    coef::AbstractArray{Tc,3},
) where {Tx<:Union{Missing,<:Number},Tc<:Union{Missing,<:Number}}
    yhat = Array{Union{Tc,Tx}}(undef, size(coef, 1), size(X, 1), size(coef, 2))
    for ch = 1:size(coef, 1)
        yhat[ch, :, :] .= X * permutedims(coef[ch, :, :], (2, 1))
    end
    return permutedims(yhat, (1, 3, 2)) # esults in: ch x time x pred
end

function predict(
    X::AbstractMatrix{Tcoef},
    coefs::AbstractMatrix{TX},
    onsets::Vector,
    timewindow::NTuple{2,Int},
) where {Tcoef,TX}
    # Create an empty array for the residuals
    # Note there is a +1 because start and stop indices are both included in the interval
    yhat = Array{Union{TX,Tcoef}}(
        undef,
        size(coefs, 1),
        timewindow[2] - timewindow[1] + 1,
        length(onsets),
    )
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


function _residuals(
    yhat::AbstractArray{Tyhat,3},
    y::AbstractMatrix{Ty},
    onsets::Vector,
    timewindow::NTuple{2,Int},
) where {Ty,Tyhat}

    resid = missings(
        Union{Ty,Tyhat},
        size(y, 1),
        timewindow[2] - timewindow[1] + 1,
        length(onsets),
    )
    for event in eachindex(onsets)

        event_range_temp = onsets[event]+timewindow[1]:onsets[event]+timewindow[2]

        # Clip indices that are outside of the design matrix (i.e. before the start or after the end)

        indices_inside_data = 0 .< event_range_temp .< size(y, 2)
        event_range = event_range_temp[indices_inside_data]


        !isempty(event_range) || @warn "event_range is empty for $event $(onsets[event])"

        # Calculate residuals
        #@debug "_residuals" size(y) size(yhat) size(event_range)
        #@debug event event_range[[1, end]]
        @views resid[:, indices_inside_data, event] .=
            y[:, event_range] .- yhat[:, indices_inside_data, event]

        # if !all(indices_inside_data)
        #    @views resid[:, .!indices_inside_data, event] .=
        #        y[:, event_range_temp[.!indices_inside_data]]
        #end
    end
    return resid
end


@traitfn residuals(
    uf::T,
    data::AbstractArray,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}} =
    _residuals(T, predict(uf), check_data(T, data))
@traitfn function residuals(
    uf::T,
    data::AbstractArray,
) where {T <: UnfoldModel; !ContinuousTimeTrait{T}}
    pred = predict(uf)
    @assert length(pred) == 1 "residuals currently not supported for multi-event MassUnivariateModels. It is not hard to implement, but I ran out of time. Write an issue if you  need this functionality"
    _residuals(T, pred[1], check_data(T, data))

end

_split_data(y::AbstractMatrix, n) = return @view(y[:, 1:n]), @view(y[:, n+1:end])
_split_data(y::AbstractArray, n) = return @view(y[:, 1:n, :]), @view(y[:, n+1:end, :])
#@traitfn 
function _residuals(::Type{T}, yhat, y) where {T<:UnfoldModel}#; ContinuousTimeTrait{T}}

    #_split_data(y::AbstractVector, n) = return @view(y[1:n]), @view(y[n+1:end])

    n_yhat = size(yhat, 2)
    n_y = size(y, 2)
    @debug n_yhat n_y
    if n_yhat >= n_y
        @debug "n_yhat > n_y" size(y) size.(_split_data(yhat, n_y))
        return y .- _split_data(yhat, n_y)[1]
    else
        @debug "n_y < n_yhat"
        yA, yB = _split_data(y, n_yhat)
        @debug size(yA) size(yB)
        res = yA .- n_y
        return cat(res, yB; dims = 2)
    end

end

predict(uf::UnfoldModel, f::FormulaTerm, args...; kwargs...) =
    predict(uf, [f], args...; kwargs...)
predict(uf::UnfoldModel, f::Vector, evts::DataFrame; kwargs...) =
    predict(uf, f, repeat([evts], length(f)); kwargs...)
predict(uf::UnfoldModel; kwargs...) = predict(uf, formulas(uf), events(uf); kwargs...)
predict(uf::UnfoldModel, f::Vector{<:FormulaTerm}; kwargs...) =
    predict(uf, f, events(uf); kwargs...)
predict(uf::UnfoldModel, evts::Vector{<:DataFrame}; overlap = false, kwargs...) =
    predict(uf, Unfold.formulas(uf), evts; overlap, kwargs...)

predict(uf::UnfoldModel, evts::DataFrame; overlap = false, kwargs...) = predict(
    uf,
    Unfold.formulas(uf),
    repeat([evts], length(Unfold.formulas(uf)));
    overlap,
    kwargs...,
)

"""
    function predict(
        uf::UnfoldModel,
        f::Vector{<:FormulaTerm},
        evts::Vector{<:DataFrame};
        overlap::Bool = true,
        kwargs...
    )
Returns a predicted ("y_hat = X*b") `Array`. 

- `uf` is an `<:UnfoldModel`
- `f` is a (vector of) formulas, typically `Unfold.formulas(uf)`, but formulas can be modified e.g. by `effects`.
- `evts` is a (vector of) events, can be `Unfold.events(uf)` to return the (possibly continuous-time) predictions of the model. Can be a custom even


## kwargs:
if `overlap = true` (default), overlap based on the `latency` column of `evts` will be simulated, or in the case of `!ContinuousTimeTrait` just X*coef is returned. 

if `overlap = false`, returns predictions without overlap (models with `ContinuousTimeTrait` (=> with basisfunction / deconvolution) only), via `predict_no_overlap`

if `keep_basis` or `exclude_basis` is defined, then `predict_partial_overlap` is called, which allows to selective introduce overlap based on specified (or excluded respective) events/basisfunctions

`epoch_to` and  `epoch_timewindow`: calculate (partial) overlap controlled predictions, but returns them at the specified `epoch_at` event, with the times `epoch_timewindow` (default is taken from the basisfunction) in samples.

`eventcolumn` can be specified as well if different from the default `event`.

Hint: all `kwargs` can be `Vector`, or if e.g. `string` types are provided, will be put into a `length==1` vector.

## Output

- If `overlap=false`, returns a 3D-Array
- If `overlap=true` and `epoch_to = nothing` (default), returns a 2D-array
- If `overlap=true` and `epoch_to != nothing`, returns a 3D array
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
    kwargs...,
    #eventcolumn = :event,
)
    @assert !(
        (!isempty(keep_basis) || !isempty(exclude_basis)) && (length(formulas(uf)) == 1)
    ) "No way to calculate partial overlap if basisnames is Any; please revise model "

    @assert !(!isempty(keep_basis) & !isempty(exclude_basis)) "choose either to keep events, or to exclude, but not both"


    coefs = coef(uf)

    if overlap
        if isnothing(epoch_to) && (isempty(keep_basis) & isempty(exclude_basis))
            @debug "full-overlap"
            if events(uf) == evts
                @debug "original design predict"
                # off-the-shelf standard continuous predict, as you'd expect
                return predict(modelmatrix(uf), coef(uf))
            else
                @debug "new dataframe predict"
                # predict of new data-table with a latency column. First generate new designmatrix, then off-the-shelf X*b predict
                @assert length(f) == length(evts) "please provide the same amount of formulas (currently $(length(f)) as event-dataframes (currently $(length(evts)))"

                X_new = _modelcols.(f, evts) |> x -> Unfold.extend_to_larger(x)
                @debug typeof(X_new) typeof(modelcols.(f, evts))
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
                kwargs...,
            )
        end
    else

        @debug "no overlap predict2"
        # no overlap, the "ideal response". We predict each event row independently. user could concatenate if they really want to :)
        return predict_no_overlap(uf, coefs, f, evts)
    end
end


"""
    $(SIGNATURES)

Returns predicted time-continuous values, but only for a subset of events. This is achieved by excluding the part of the designmatrix that belongs to the basisfunctions/events you do not want to have in your model.

Typically called via `predict`, for configuration, keyword-arguments and usage see there.

One difference is, that we require the `coefs(uf::UnfoldModel)` already exctracted. 

Due to the time-continuous nature, running it with a model not containing the `ContinuousTimeTrait` it will throw an error.
"""
@traitfn predict_partial_overlap(
    uf::T,
    args;
    kwargs...,
) where {T <: UnfoldModel; !ContinuousTimeTrait{T}} =
    error("can't have partial overlap without Timecontinuous model")

@traitfn function predict_partial_overlap(
    uf::T,
    coefs,
    evts;
    keep_basis = [],
    exclude_basis = [],
    epoch_to = nothing,
    epoch_timewindow = nothing,
    eventcolumn = :event,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}}
    @assert !(!isempty(keep_basis) & !isempty(exclude_basis)) "can't have no overlap & specify keep/exclude at the same time. decide for either case"
    # Partial overlap! we reconstruct with some basisfunctions deactivated
    if !isempty(keep_basis)
        if !isa(keep_basis, Vector) # Check if keep_basis is a vector
            basisnames = [keep_basis]
        else
            basisnames = keep_basis
        end
        @assert !isempty(intersect(basisname(Unfold.formulas(uf)), basisnames)) "Couldn't find (any of) $keep_basis in the models basisnames; you can check which basisnames are available in your model using Unfold.basisname(Unfold.formulas(uf))"
    else
        basisnames = basisname(Unfold.formulas(uf))
        basisnames = setdiff(basisnames, exclude_basis)
    end
    #    @debug basisname
    X_view = matrix_by_basisname(modelmatrix(uf), uf, basisnames)
    coefs_view = matrix_by_basisname(coefs, uf, basisnames)


    if isnothing(epoch_to)
        return predict(X_view, coefs_view)

    else
        #        @debug typeof(evts)
        timewindow =
            isnothing(epoch_timewindow) ? calc_epoch_timewindow(uf, epoch_to) :
            epoch_timewindow

        find_lat =
            x -> x[
                epoch_to == Any ? (1:size(x, 1)) : x[:, eventcolumn] .== epoch_to,
                :latency,
            ]
        latencies = vcat(find_lat.(evts)...)
        #ix_event = vcat(evts...)[:, eventcolumn] .== epoch_to
        #latencies = evts[ix_event].latency
        #@debug latencies[end] size(X_view) timewindow
        return predict(X_view, coefs_view, latencies, (timewindow[1], timewindow[end]);)
    end

end


"""
    $(SIGNATURES)
in ContinuousTime case (typically the deconvolution model), we return idealized predictions without overlap between events.

in the Not-ContinuousTime case (typically the MassUnivariate model), we return predictions for each event independently. In that case, the function is unfortunately a missnomer, as overlap cannot be removed from mass-univariate models.

"""
@traitfn function predict_no_overlap(
    uf::T,
    coefs,
    f::Vector,
    evts::Vector,
) where {T <: UnfoldModel; !ContinuousTimeTrait{T}}
    @debug "Not ContinuousTime yhat, Array"
    X = _modelcols.(f, evts)
    @debug typeof(X)
    # figure out which coefficients belong to which event
    Xsizes = size.(X, Ref(2))
    Xsizes_cumsum = vcat(0, cumsum(Xsizes))
    indexes = [(Xsizes_cumsum[ix]+1):Xsizes_cumsum[ix+1] for ix = 1:length(f)]

    # split up the coefs accordingly
    coArray = [@view coefs[:, :, ix] for ix in indexes]
    @debug typeof(collect(X[1])), typeof(collect(coArray[1])) size(X) size(coArray)
    return predict.(X, coArray)
end

@traitfn function predict_no_overlap(
    uf::T,
    coefs,
    f::Vector,
    evts::Vector,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}}

    has_missings = false
    yhat = Array{eltype(coefs)}[]
    @debug eltype(coefs) typeof(yhat)
    for (fi, e) in zip(f, deepcopy(evts))

        # put latency to 1 in case no basis shift
        e.latency .= -fi.rhs.basisfunction.shift_onset + 1 #sum(times(fi) .<= 0)
        #e.latency .= 1
        X_singles = map(x -> _modelcols(fi, DataFrame(x)), eachrow(e))
        coefs_view = matrix_by_basisname(coefs, (uf), (basisname([fi])))

        yhat_single = missings(
            eltype(coefs_view),
            size(coefs_view, 1),
            maximum(size.(X_singles, 1)),
            length(X_singles),
        )
        @debug coefs_view
        for ev = 1:length(X_singles)
            p = predict(X_singles[ev], coefs_view)
            if has_missings || isa(p, AbstractArray{<:Union{<:Missing,<:Number}})
                # if outside range predictions happen for spline predictors, we have to allow missings
                #                @debug "yeah, missings..."
                has_missings = true
                yhat_single = allowmissing(yhat_single)
            end
            yhat_single[:, 1:length(p), ev] .= p
        end


        yhat = combine_yhat!(yhat, yhat_single)

    end

    @debug typeof(yhat) typeof.(yhat)
    @debug typeof(vcat(yhat))
    return yhat
end



"""
    combine_yhat(list,single)
combines single into list, if either list or single contains missing, automatically casts the respective counter-part to allow missings as well
"""
function combine_yhat!(yhat::Vector{<:Array{T}}, yhat_single::Array{T}) where {T}
    @debug typeof(yhat) typeof.(yhat) typeof(yhat_single)
    push!(yhat, yhat_single)
end
combine_yhat!(yhat::Vector{Array{Union{Missing,<:Number}}}, yhat_single) =
    push!(yhat, allowmissing(yhat_single))
function combine_yhat!(
    yhat::Vector{<:Array{<:Number}},
    yhat_single::Array{T},
) where {T<:Union{Missing,<:Number}}
    @debug "new yhat"
    yhat_new = Array{T}[]
    combine_yhat!.(Ref(yhat_new), yhat)
    combine_yhat!(yhat_new, yhat_single)

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
    return (1:Unfold.height(basisfunction)) .- 1 .+ Unfold.shift_onset(basisfunction) # start at 0
end

"""
    get_basis_indices(uf, basisnames::Vector)
returns a boolean vector with length spanning all coefficients, which coefficient is defined by `basisnames` (vector of names)
"""
get_basis_indices(uf, basisnames::Vector) =
    reduce(vcat, Unfold.get_basis_names(uf)) .∈ Ref(basisnames)
get_basis_indices(uf, basisnames) = get_basis_indices(uf, [basisnames])

"""
    predicttable(model<:UnfoldModel,events=Unfold.events(model),args...;kwargs...)
Shortcut to call efficiently call (pseudocode) `result_to_table(predict(...))`.

Returns a tidy DataFrame with the predicted results. Loops all input to `predict`, but really only makes sense to use if you specify either:

`overlap = false` (the default) or `epoch_to = "eventname"`.

"""
function predicttable(
    model::UnfoldModel,
    events::Vector = Unfold.events(model),
    args...;
    overlap = false,
    kwargs...,
)
    p = predict(
        model,
        Unfold.formulas(model),
        events,
        args...;
        overlap = overlap,
        kwargs...,
    )
    return result_to_table(model, p, events)
end
predicttable(model, events::DataFrame) = predicttable(model, [events])
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
    #    @debug all_times length(all_times) length.(all_times)
    @assert all(all_times .== all_times[1:1]) "all times need to be equal in a mass univariate model"
    return all_times
end

