

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
predict(uf::UnfoldModel, evts::DataFrame; overlap = false, kwargs...) =
    predict(uf, Unfold.formulas(uf), evts; overlap, kwargs...)

# Predict new 
function predict(
    uf,
    f::Vector{<:FormulaTerm},
    evts::Vector{<:DataFrame};
    overlap = true,
    keep_basis = [],
    exclude_basis = [],
    epoch_to = nothing,
    epoch_timewindow = nothing,
    eventcolumn = :event,
)
    @assert !(!isempty(keep_basis) & !isempty(exclude_basis)) "choose either to keep events, or to exclude, but not both"
    #@assert overlap == false & !isempty(keep_basis) & !isempty(exclude_basis) "can't have no overlap & specify keep/exclude at the same time. decide for either case"

    #evts = reduce(vcat, evts) # XXX for now only

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
            return predict_partial_overlap()
        end
    else

        @debug "no overlap predict2"
        # no overlap, the "ideal response". We predict each event row independently. user could concatenate if they really want to :)
        @debug size(f), size(evts)
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
    # Partial overlap! we reconstruct with some basisfunctions deactivated
    if !isempty(keep_basis)
        basisnames = keep_basis
    else
        basisnames = basisname(f)
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
    #    @debug ix
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

#predict_to_table(model, eff::AbstractArray, events::DataFrame) =
#    predict_to_table(eff, repeat([events], length(formulas(model))), times(model))
predict_to_table(model, eff::AbstractArray, events::Vector{<:DataFrame}) =
    predict_to_table(eff, events, times(model), eventnames(model))


eventnames(model::UnfoldModel) = first.(design(model))
function predict_to_table(
    eff::AbstractArray,
    events::Vector,
    times::Vector,
    eventnames::Vector,
)
    @assert length(eventnames) == length(events) == length(times)
    times_vecs = repeat.(poolArray.(times), size.(events, 1))
    # times_vecs = map((x, n) -> repeat(x, outer = n), poolArray.(times), size.(events, 1))
    # init the meta dataframe
    @debug typeof(times_vecs) size(times_vecs)
    data_list = []
    for (single_events, single_eff, single_times, single_eventname) in
        zip(events, eff, times_vecs, eventnames)

        #=
        metadata = FlexTable(
            time = single_times |> poolArray,
            eventname = repeat([single_eventname] |> poolArray, length(single_times)),
        ) # used to be DataFrames, but dataFrames looses type-information of pooledarray upon concatenation with vcat

        for c in names(single_events)
            setproperty!(
                metadata,
                Symbol(c),
                Vector{eltype(single_events[1, c])}(undef, length(metadata.time)),
            )
            #metadata[:, c] .= single_events[1, c] # assign first element in order to have same column type
        end
        =#
        #        @debug metadata.conditionA

        ntimes = size(single_eff, 2)
        evts_tbl = repeat(Table(single_events), inner = (ntimes))
        time_event = Table(
            time = single_times |> poolArray,
            eventname = repeat([single_eventname] |> poolArray, length(single_times)),
        )
        @debug size(time_event) length(single_times) ntimes
        metadata = Table(evts_tbl, time_event)

        #metadata = repeat(metadata, ntimes)
        #for c in names(single_events)
        #    @debug c
        #    for row = 1:size(single_events, 1)
        #        rowIx = (1.0 .+ (row - 1) .* ntimes) .+ range(1.0, length = ntimes) .- 1
        #        getproperty(metadata, Symbol(c))[Int64.(rowIx)] .= single_events[row, c]
        #        #metadata[Int64.(rowIx),c] .= single_events[row, c]
        #    end
        #    @debug metadata.conditionA[end]
        #end

        single_data = Table(
            yhat = single_eff[:],#vec(reshape(single_eff, :, 1)),
            channel = repeat(
                (1:size(single_eff, 1)) |> poolArray,
                outer = size(single_eff, 2) * size(single_eff, 3),
            ),
        )

        # single_data = hcat(single_data, repeat(metadata, size(single_eff, 1)))
        @debug size(metadata) size(single_eff) size(single_data)
        single_data_comb =
            map(merge, single_data, repeat(metadata, inner = size(single_eff, 1)))
        push!(data_list, single_data_comb)
    end
    return reduce(vcat, data_list) |> DataFrame


end

# in case the formula is not an array
#predict(model::UnfoldModel,formulas::AbstractTerm) = predict(model,[formulas])


#=
# special case if one formula is defined but in an array => multiple events
@traitfn function predict(
    model::T,
    formulas::AbstractArray,
    singleevents::DataFrame,
) where {T <: UnfoldModel; !ContinuousTimeTrait{T}}



end
=#

#=
predict(model::UnfoldLinearModel, formulas::FormulaTerm, singleevents::DataFrame) =
    predict(model, formulas.rhs, singleevents)
predict(model::UnfoldLinearModelContinuousTime, formulas::FormulaTerm, singleevents) =
    predict(model, formulas.rhs, singleevents)

=#

#=
@traitfn function predict(
    model::T,
    formulas::AbstractTerm,
    singleevents,
) where {T <: UnfoldModel; !ContinuousTimeTrait{T}}
    @debug "Not ContinuousTime yhat, single Term"
    X = modelcols(formulas, singleevents)
    return predict(model, X)
end
=#
#=
@traitfn function predict(
    m::T,
    f::AbstractArray,
    events,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}}
    @debug f
    predict(m, blockdiag(predict.(Ref(m), f, Ref(events))...))  #  blockdiag([Xsingle1,Xsingle2,Xsingle3]...)
end

@traitfn predict(
    model::T,
    f::MatrixTerm,
    events,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}} = predict(m, f.terms, events) # terms[1] or not?!?

@traitfn predict(
    model::T,
    f::FormulaTerm,
    events,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}} = predict(m, f.rhs, events)

=#

#=
@traitfn function predict(
    model::T,
    f::TimeExpandedTerm,
    events,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}}

    # find out how long each designmatrix is
    n_range = length(times(f.basisfunction))
    # find out how much to shift so that X[1,:] is the the first "sample"
    n_negative = f.basisfunction.shift_onset
    # generate the correct eventfields (default: latencies)
    events = deepcopy(events)

    events[:, f.eventfields[1]] =
        range(-n_negative + 1, step = n_range, length = size(events, 1))
    # get the model
    Xsingle = modelcols(f, events)

    timesSingle = times(f.basisfunction)

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

    if typeof(f.basisfunction) <: FIRBasis
        keep = ones(size(Xsingle, 1))
        keep[range(length(timesSingle), size(Xsingle, 1), step = length(timesSingle))] .= 0
        Xsingle = Xsingle[keep.==1, :]
    end
    return Xsingle
    # don't include, via broadcast: append!(X, [Xsingle])
end

=#


#=
function predict(
    model::UnfoldLinearModelContinuousTime,
    X::AbstractArray{T,2};
    times = nothing,
) where {T<:Union{Missing,<:Number}}

    #@debug size(X), size(coef(model))
    #@debug typeof(X), typeof(coef(model))
    yhat = X * coef(model)'
    #coefs = coef(model)
    #@tullio yhat[i, k] := X[i, j] * coefs[k, j]
    return yhat

end
=#

#=
# kept for backwards compatability
@traitfn predict(
    model::U,
    X::AbstractArray{T,2};
    kwargs...,
) where {T<:Union{Missing,<:Number},U<:UnfoldModel;!ContinuousTimeTrait{U}} =
    predict(coef(model), X; kwargs...)




#function yhat_mult(X::AbstractArray{T,2}, coef) where {T<:Number}
#
#    @tullio yhat[ch, a, b] := X[a, trial] * coef[ch, b, trial]
#    return yhat
# end
function yhat_mult(X::AbstractArray{T,2}, coef) where {T<:Union{Missing,<:Number}}
    yhat = Array{T}(undef, size(coef, 1), size(X, 1), size(coef, 2))
    for ch = 1:size(coef, 1)
        yhat[ch, :, :] .= X * permutedims(coef[ch, :, :], (2, 1))
    end
    return yhat
end
function predict(
    coef::AbstractArray{T1,3},
    X::AbstractArray{T2,2};
    kwargs...,
) where {T1<:Union{Missing,<:Number},T2<:Union{Missing,<:Number}}
    # function that calculates coef*designmat, but in the ch x times x coef vector
    # setup the output matrix, has to be a matrix
    # then transforms it back to 2D matrix times/coef x ch to be compatible with the timecontinuous format
    @debug "type X", typeof(X)

    try

        X = disallowmissing(X)

    catch
    end
    yhat = yhat_mult(X, coef)


    # bring the yhat into a ch x yhat format
    yhat = reshape(permutedims(yhat, (1, 3, 2)), size(yhat, 1), :)
    yhat = permutedims(yhat, (2, 1))
    return yhat
end

=#

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




#=
@traitfn yhat_nranges(
    model::T,
    formulas::Array{<:FormulaTerm},
    events,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}} =
    yhat_nranges.(Ref(model), formulas, Ref(events)) # todo: add hcat? vcat? who knows


@traitfn function yhat_nranges(
    model::T,
    f::FormulaTerm,
    events,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}}
    n_range = length(times(f.rhs.basisfunction))
    if typeof(f.rhs.basisfunction) <: FIRBasis
        n_range = n_range - 1
    end
    return range(1, step = n_range, length = size(events, 1))
end
=#