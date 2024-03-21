import StatsBase.predict

using DocStringExtensions: format



# opt-in

# Most important function ;-)
#predict(X::AbstractMatrix, coefs::AbstractMatrix) = coefs * X'
predict(X::AbstractMatrix, coefs::AbstractMatrix) = (X * coefs')'

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

# predict new data (for effects)
# predict without overlap
# Predict continuos data
residuals(uf, data::AbstractArray) = data .- predict(uf)

# Predict new 
function predict(
    uf,
    evts::DataFrame = reduce(vcat, events(uf));
    overlap = true,
    keep_basis = [],
    exclude_basis = [],
    epoch_to = nothing,
    epoch_timewindow = nothing,
    eventcolumn = :event,
)
    @assert !(!isempty(keep_basis) & !isempty(exclude_basis)) "choose either to keep events, or to exclude, but not both"
    #@assert overlap == false & !isempty(keep_basis) & !isempty(exclude_basis) "can't have no overlap & specify keep/exclude at the same time. decide for either case"


    coefs = coef(uf)
    # generate X_new
    # X_new*b => continuierlich yhats::Matrix, channel x continuostime
    if overlap
        if isempty(keep_basis) & isempty(exclude_basis)
            # return full-overlap
            if events(uf) == evts
                return predict(uf, modelmatrix(uf))
            else
                X_new = modelcols.(formulas(uf), Ref(evts)) |> Unfold.extend_to_larger
                return predict(X_new, coefs)
            end
        else
            if !isempty(keep_basis)
                basisnames = keep_basis
            else
                basisnames = basisname(formulas(uf))
                basisnames = setdiff(basisnames, exclude_basis)
            end

            ix = get_basis_indices(uf, basisnames)



            X_view = @view(modelmatrix(uf)[:, ix])
            coefs_view = @view(coefs[:, ix])
            if isnothing(epoch_to)
                return predict(X_view, coefs_view)

            else
                timewindow =
                    isnothing(epoch_timewindow) ? calc_epoch_timewindow(uf, epoch_to) :
                    epoch_timewindow
                @debug "timewindow" timewindow
                latencies = evts[evts[:, eventcolumn].==epoch_to, :latency]
                return predict(
                    X_view,
                    coefs_view,
                    latencies,
                    (timewindow[1], timewindow[end]);
                )
            end

        end

    else
        # generate [X_single]
        # each X_single * b => "epochiert", ohne overlap [yhats::Matrix] jeweils channel x basistimes

        # XXX improvement possible here, call modelcols only for the appropriate events "subset(evts")"
        X_singles = map(x -> modelcols(formulas(uf), DataFrame(x)), eachrow(evts))
        return predict.(X_singles, Ref(coefs))
    end
end


"""
returns an integer range with the samples around `epoch_event` as defined in the corresponding basisfunction

"""

function calc_epoch_timewindow(uf, epoch_event)
    basis_ix = findfirst(Unfold.basisname(formulas(uf)) .== epoch_event)
    basisfunction = formulas(uf)[basis_ix].rhs.basisfunction

    @debug "calc_epoch_timewindows" basis_ix,
    epoch_event,
    height(basisfunction),
    shift_onset(basisfunction)
    return (1:Unfold.height(basisfunction)) .+ Unfold.shift_onset(basisfunction)
end

get_basis_indices(uf, basisnames) = Unfold.get_basis_name(uf) .∈ Ref(basisnames)




function predicttable(model::UnfoldModel, events)
    # make a copy of it so we don't change it outside the function
    newevents = copy(events)

    formulas = Unfold.formulas(model)
    if typeof(formulas) <: FormulaTerm
        formulas = [formulas]
    end
    if isa(model, UnfoldLinearModel)
        eff = predict(model, formulas[1], newevents)
        timesVec = gen_timeev(times(model), size(newevents, 1))
    else
        @debug size(formulas), typeof(model)


        eff = predict(model, formulas, newevents)
        @debug typeof(model)
        timesVec = yhat_timevec(model, formulas, newevents)
        fromTo = yhat_nranges(model, formulas, newevents) # formerly fromTo == n_ranges
    end


    # init the meta dataframe
    metaData = DataFrame([:time => vcat(timesVec...), :eventname => ""])

    for c in names(newevents)
        metaData[:, c] .= newevents[1, c] # assign first element in order to have same column type
    end
    if isa(model, UnfoldLinearModel)
        # for mass univariate we can make use of the knowledge that all events have the same length :)
        ntimes = size(coef(model), 2)
        for c in names(newevents)
            for row = 1:size(newevents, 1)
                rowIx = (1.0 .+ (row - 1) .* ntimes) .+ range(1.0, length = ntimes) .- 1
                metaData[Int64.(rowIx), c] .= newevents[row, c]
            end
        end

    else


        # shift variable to keep track of multiple basisfunctions
        shift = 0
        # for each basis function
        @debug "fromTo" fromTo
        for (bIx, n_range) in enumerate(fromTo)
            @debug "n_range" typeof(n_range[1:end]), size(events)
            #basistime = range(1, step = n_range, length = size(events, 1))

            # go through all predictors
            for (i, fstart) in enumerate(n_range[1:end])
                fend = n_range[i] + n_range.step .- 1
                #fend = fstart + basistime.step - 1

                # couldn't figure out how to broadcast everything directly (i.e out[fstart:fend,names(newevents)] .= newevents[i,:])
                # copy the correct metadata
                @debug fstart, fend
                for j = fstart:fend
                    metaData[j+shift, names(newevents)] = newevents[i, :]
                end
                # add basisfunction name
                metaData[shift.+(fstart:fend), :eventname] .=
                    formulas[bIx].rhs.basisfunction.name

            end
            # the next meta data has to be at the end
            shift += n_range[end] - 1 + n_range.step

        end# opt-in
    end



    out = DataFrame([:yhat => vec(reshape(eff, :, 1))])
    nchannel = size(eff, 2)

    out.channel = repeat(1:nchannel, inner = size(eff, 1))
    out = hcat(out, repeat(metaData, nchannel))
    return out
end

# in case the formula is not an array
#predict(model::UnfoldModel,formulas::AbstractTerm) = predict(model,[formulas])

# special case if one formula is defined but in an array => multiple events
@traitfn function predict(
    model::T,
    formulas::AbstractArray,
    newevents::DataFrame,
) where {T <: UnfoldModel; !ContinuousTimeTrait{T}}
    @debug "Not ContinuousTime yhat, Array"
    X = modelcols.(formulas, Ref(newevents))

    co = coef(model)
    Xsizes = size.(X, Ref(2))
    Xsizes_cumsum = vcat(0, cumsum(Xsizes))

    indexes = [(Xsizes_cumsum[ix]+1):Xsizes_cumsum[ix+1] for ix = 1:length(formulas)]

    coArray = [co[:, :, ix] for ix in indexes]


    return yhat.(coArray, X)


end


predict(model::UnfoldLinearModel, formulas::FormulaTerm, newevents::DataFrame) =
    predict(model, formulas.rhs, newevents)
predict(model::UnfoldLinearModelContinuousTime, formulas::FormulaTerm, newevents) =
    predict(model, formulas.rhs, newevents)

@traitfn function predict(
    model::T,
    formulas::AbstractTerm,
    newevents,
) where {T <: UnfoldModel; !ContinuousTimeTrait{T}}
    @debug "Not ContinuousTime yhat, single Term"
    X = modelcols(formulas, newevents)
    return predict(model, X)
end


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
    return all_times[1]
end


function gen_timeev(timesVec, nRows)
    timesVec = repeat(timesVec, nRows)
    return timesVec
end




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

@traitfn yhat_timevec(
    model::T,
    formulas::Array{<:FormulaTerm},
    events,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}} =
    yhat_timevec.(Ref(model), formulas, Ref(events))


@traitfn function yhat_timevec(
    model::T,
    f::FormulaTerm,
    events,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}}
    @debug typeof(f)
    # keep track of the times
    timesSingle = times(f.rhs.basisfunction)
    # see yhat blabla why this is necessary
    if typeof(f.rhs.basisfunction) <: FIRBasis
        timesSingle = timesSingle[1:end-1]
    end
    return repeat(timesSingle, size(events, 1))
end
