# opt-in

# Most important function ;-)
predict(X, coefs) = X * coefs'

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

        event_range_temp = start_idx[event]:stop_idx[event]

        # Clip indices that are outside of the design matrix (i.e. before the start or after the end)
        #indices_inside_data = 0 .< event_range_temp .< size(data, 2)
        indices_inside_data = 0 .< event_range_temp .< size(X, 1)
        event_range = event_range_temp[indices_inside_data]

        # Calculate residuals
        yhat[:, indices_inside_data, event] .= predict(@view(X[event_range, :]), coefs)
    end
end

# predict new data (for effects)
# predict without overlap
# Predict continuos data
residuals(uf, data::AbstractArray) = data .- predict(uf)

# Predict new 
function predict(
    uf,
    evts::DataFrame = events(uf);
    overlap = true,
    keep_basis = [],
    exclude_basis = [],
    epoch_to = nothing,
    epoch_timewindow = nothing,
    eventcolumn = :event,
)
    @assert !(!isempty(keep_basis) & !isempty(exclude_basis)) "choose either to keep events, or to exclude, but not both"
    @assert overlap == false & !isempty(keep_basis) & !isempty(exclude_basis) "can't have no overlap & specify keep/exclude at the same time. decide for either case"


    coefs = coef(uf)
    # generate X_new
    # X_new*b => continuierlich yhats::Matrix, channel x continuostime
    if overlap
        if isempty(keep_basis) & isempty(exclude_basis)
            # return full-overlap
            if events(uf) == evts
                return predict(uf, modelmatrix(uf))
            else
                X_new = modelcols.(formulas(uf), Ref(evts)) |> Unfold.equalize_lengths
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
            if !isnothing(epoch_to)
                return predict(X_view, coefs_view)

            else
                timewindow =
                    isnothing(epoch_timewindow) ? calc_epoch_timewindow(uf, epoch_event) :
                    epoch_timewindow

                latencies = evts[evts.eventcolumn.==epoch_event, :latencies]
                return predict(X_view, coefs_view, latencies, timewindow;)
            end

        end

    else
        # generate [X_single]
        # each X_single * b => "epochiert", ohne overlap [yhats::Matrix] jeweils channel x basistimes

        # XXX improvement possible here, call modelcols only for the appropriate events "subset(evts")"
        X_singles = modelcols.(Ref(formulas(uf)), eachrow(evts))
        return predict.(X_singles, Ref(coefs))
    end
end


"""
returns an integer range with the samples around `epoch_event` as defined in the corresponding basisfunction

"""

function calc_epoch_timewindow(uf, epoch_event)
    evts = events(uf)
    basis_ix = findfirst(Unfold.basisname(formulas(uf)) .== epoch_event)
    basisfunction = formulas(uf)[basis_ix].rhs.basisfunction

    epoch_timewindow = 1:Unfold.height(basisfunction).-Unfold.shift_onset(basisfunction)
end

get_basis_indices(uf, basisnames) = Unfold.get_basis_name(uf) .∈ Ref(basisnames)

