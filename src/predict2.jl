# opt-in

# Most important function ;-)
predict(X, coefs) = X * coefs'


# predict new data (for effects)
# predict without overlap
# Predict continuos data
residuals(uf, data::AbstractArray) = data .- predict(uf)


# Predict new 
function predict(uf, evts::DataFrame; overlap = false)
    coefs = coef(uf)
    # generate X_new
    # X_new*b => continuierlich yhats::Matrix, channel x continuostime
    if overlap
        X_new = modelcols(formulas(uf), evts)
        return predict(X_new, coefs)
    else

        # generate [X_single]
        # each X_single * b => "epochiert", ohne overlap [yhats::Matrix] jeweils channel x basistimes

        # XXX improvement possible here, call modelcols only for the appropriate events "subset(evts")"
        X_singles = modelcols.(Ref(formulas(uf)), eachrow(evts)) |> equalize_lengths
        return predict.(X_singles, Ref(coefs))
    end
end


# predict partial overlap & full overlap
# predict(uf) = predict(uf, modelmatrix(uf))
function predict(
    uf;
    keep_basis = nothing,
    exclude_basis = nothing,
    epoch_to = nothing,
    epoch_timewindow = nothing,
)
    @assert !(!isnothing(keep_basis) & !isnothing(exclude_basis)) "choose either to keep events, or to exclude, but not both"
    if isnothing(keep_basis) & isnothing(exclude_basis)
        yhat = predict(uf, modelmatrix(uf))
    end

    # some logic to produce the correct basisnames vector out of keep_basis / remove_basis
    # XXX
    # XXX
    #
    basisnames = "xxx"

    yhat = predict(uf, basisnames)

    if !isnothing(epoch_to)
        return yhat
    else
        @assert !is_vector(epoch_event) & isnothing(epoch_timewindow)
        timewindow =
            isnothing(epoch_timewindow) ? length(basisfunction[epoch_event]) :
            epoch_timewindow
        latencies = events(uf)[epoch_event, :latencies]
        return Unfold.epoch(yhat, latencies, timewindow, 1)

    end
end


"""
# predict(uf,basisnames)
basisnames could be [:stimulus,:fixation]
"""
@traitfn function predict(
    uf,
    basisnames::Vector{<:Symbol,<:String},
) where {T <: UnfoldModel; ContinuousTimeTrait{T}}
    coefs = coef(uf)
    X = modelmatrix(uf)
    ix = get_basis_indices(formulas(uf), basisnames)
    return predict(@view(X[:, ix]), @view(coefs[:, ix]))
end


