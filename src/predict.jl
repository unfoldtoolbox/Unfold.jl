import StatsBase.predict

using DocStringExtensions: format


function StatsBase.predict(model::UnfoldModel, events)
    # make a copy of it so we don't change it outside the function
    newevents = copy(events)

    formulas = Unfold.formulas(model)
    if typeof(formulas) <: FormulaTerm
        formulas = [formulas]
    end
    if typeof(model) == UnfoldLinearModel
        eff = yhat(model, formulas[1], newevents)
        timesVec = gen_timeev(times(model), size(newevents, 1))
    else
        @debug size(formulas), typeof(model)
        fromTo, timesVec, eff = yhat(model, formulas, newevents)
    end


    # init the meta dataframe
    metaData = DataFrame([:time => vcat(timesVec...), :eventname => ""])

    for c in names(newevents)
        metaData[:, c] .= newevents[1, c] # assign first element in order to have same column type
    end
    if typeof(model) == UnfoldLinearModel
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

        for (bIx, basisfun) in enumerate(fromTo)

            # go through all predictors
            for (i, fstart) in enumerate(basisfun[1:end])

                fend = basisfun[i] + basisfun.step - 1

                # couldn't figure out how to broadcast everything directly (i.e out[fstart:fend,names(newevents)] .= newevents[i,:])
                # copy the correct metadata

                for j = fstart:fend
                    metaData[j+shift, names(newevents)] = newevents[i, :]
                end
                # add basisfunction name
                metaData[shift.+(fstart:fend), :eventname] .=
                    formulas[bIx].rhs.basisfunction.name

            end
            # the next meta data has to be at the end
            shift += basisfun[end] - 1 + basisfun.step

        end
    end


    out = DataFrame([:yhat => vec(reshape(eff, :, 1))])
    nchannel = size(eff, 2)

    out.channel = repeat(1:nchannel, inner = size(eff, 1))
    out = hcat(out, repeat(metaData, nchannel))
    return out
end

# in case the formula is not an array
#yhat(model::UnfoldModel,formulas::AbstractTerm) = yhat(model,[formulas])

# special case if one formula is defined but in an array => multiple events
@traitfn function yhat(
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


yhat(model::UnfoldLinearModel, formulas::FormulaTerm, newevents::DataFrame) =
    yhat(model, formulas.rhs, newevents)
yhat(model::UnfoldLinearModelContinuousTime, formulas::FormulaTerm, newevents) =
    yhat(model, formulas.rhs, newevents)

@traitfn function yhat(
    model::T,
    formulas::AbstractTerm,
    newevents,
) where {T <: UnfoldModel; !ContinuousTimeTrait{T}}
    @debug "Not ContinuousTime yhat, single Term"
    X = modelcols(formulas, newevents)
    return yhat(model, X)
end


yhat(model::UnfoldLinearModelContinuousTime, formulas::MatrixTerm, events) =
    yhat(model, formulas.terms, events)

@traitfn function yhat(
    model::T,
    formulas,
    events,
) where {T <: UnfoldModel; ContinuousTimeTrait{T}}
    @debug "Yes ContinuousTime yhat"
    X = AbstractArray[]
    fromTo = []
    timesVec = []
    for f in formulas

        # due to many reasons, there can be several different options of how formula arrives here - not easy to catch via multiple dispatch
        if !(isa(f, TimeExpandedTerm))
            if !(isa(f, MatrixTerm))
                f = f.rhs
                #elseif !(isa(f,Effects.TypicalTerm))

            else
                f = f.terms[1]
            end

        end
        # find out how long each designmatrix is
        n_range = length(times(f.basisfunction))
        # find out how much to shift so that X[1,:] is the the first "sample"
        n_negative = f.basisfunction.shiftOnset
        # generate the correct eventfields (default: latencies)
        events[:, f.eventfields[1]] =
            range(-n_negative + 1, step = n_range, length = size(events, 1))



        # get the model
        @debug f, first(events)
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

            keep = fill(true, size(Xsingle, 1))
            keep[range(
                length(timesSingle),
                size(Xsingle, 1),
                step = length(timesSingle),
            )] .= false

            Xsingle = Xsingle[keep, :] # unfortunately @view is not compatible with blockdiag of SparseArrays.
            timesSingle = timesSingle[1:end-1]
            n_range = n_range - 1
        end
        # combine designmats
        append!(X, [Xsingle])

        # keep track of how long each event is
        append!(fromTo, [range(1, step = n_range, length = size(events, 1))])
        # keep track of the times
        append!(timesVec, repeat(timesSingle, size(events, 1)))

    end

    # Concat them, but without introducing any overlap, we want the isolated responses
    Xconcat = blockdiag(X...)

    # calculate yhat
    eff = yhat(model, Xconcat)

    # output is a bit ugly, but we need the other two vectors as well. Should be refactored at some point - but not right now ;-) #XXX
    return fromTo, timesVec, eff

end
function yhat(
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
@traitfn yhat(
    model::U,
    X::AbstractArray{T,2};
    kwargs...,
) where {T<:Union{Missing,<:Number},U<:UnfoldModel;!ContinuousTimeTrait{U}} =
    yhat(coef(model), X; kwargs...)




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
function yhat(
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
