# misc functions
function epoch(; data, tbl, τ, sfreq, kwargs...)
    epoch(data, tbl, τ, sfreq; kwargs...)
end

"""
    epoch(data::Array{T,1},tbl::DataFrame,τ::Tuple/Vector,sfreq;kwargs...,

Basic function to epoch data; all input also available as kwargs.

Additional kwarg: `eventtime`=:latency, which defines the column in `tbl` that is used to cut the data (in samples). For uneven sample-times we use `round()``
 
"""
function epoch(data::Array{T,1}, tbl, τ, sfreq; kwargs...) where {T<:Union{Missing,Number}}
    data_r = reshape(data, (1, :))
    epoch(data_r, tbl, τ, sfreq; kwargs...)
end

function epoch(data::Matrix, tbl::DataFrame, τ::Vector, sfreq; kwargs...)
    return epoch(data, tbl, (τ[1], τ[2]), sfreq; kwargs...)
end

function epoch(
    data::Matrix,
    tbl::DataFrame,
    τ::Tuple{Number,Number},
    sfreq;
    eventtime::String = "latency",
)
    return epoch(data, tbl, τ, sfreq; eventtime = Symbol(eventtime))
end

function epoch(
    data::Array{T,2},
    tbl::DataFrame,
    τ::Tuple{Number,Number},
    sfreq;
    eventtime::Symbol = :latency,
) where {T<:Union{Missing,Number}}
    # data: channels x times

    # partial taken from EEG.jl

    numEpochs = size(tbl, 1)
    τ = round_times(τ, sfreq)

    times = range(τ[1], stop = τ[2], step = 1 ./ sfreq)
    lenEpochs = length(times)
    numChans = size(data, 1)
    epochs = Array{Union{Missing,Float64}}(
        missing,
        Int(numChans),
        Int(lenEpochs),
        Int(numEpochs),
    )


    # User feedback
    @debug "Creating epochs: $numChans x $lenEpochs x $numEpochs"

    for si = 1:size(tbl, 1)
        #eventonset = tbl[si,eventtime] # in samples
        #d_start = eventonset
        d_start = Int(round(tbl[si, eventtime]) + times[1] .* sfreq)
        d_end = Int(round(tbl[si, eventtime]) + times[end] .* sfreq)

        e_start = 1
        e_end = lenEpochs
        #println("d: $(size(data)),e: $(size(epochs)) | $d_start,$d_end,$e_start,$e_end | $(tbl[si,eventtime])")
        if d_start < 1
            e_start = e_start + (-d_start + 1)
            d_start = 1
        end
        if d_end > size(data, 2)
            e_end = e_end - (d_end - size(data, 2))
            d_end = size(data, 2)
        end
        #println("d: $(size(data)),e: $(size(epochs)) | $d_start,$d_end,$e_start,$e_end | $(tbl[si,eventtime])")
        epochs[:, e_start:e_end, si] = data[:, d_start:d_end]
    end
    return (epochs, times)
end

function round_times(τ, sfreq)
    # function to round τ to sfreq samples. This specifies the epoch length.
    # its a function to be the same for epoch & timeexpanded analyses
    return round.(τ .* sfreq) ./ sfreq

end
function dropMissingEpochs(X, y)
    missingIx = .!any(ismissing.(y), dims = (1, 2))
    goodIx = dropdims(missingIx, dims = (1, 2))
    return X[goodIx, :], Array{Float64}(y[:, :, goodIx])
end


"""
$(SIGNATURES)
Flatten a 1D array from of a 2D/3D array. Also drops the empty dimension
"""
function linearize(x::AbstractArray{T,N}) where {T,N}
    # used in condense_long to generate the long format
    return x[:] #dropdims(reshape(x, :, 1), dims = 2)::AbstractArray{T,1}
end
function linearize(x::String)
    return x
end



"""
$(SIGNATURES)
Equates the length of data and designmatrix by cutting the shorter one

The reason we need this is because when generating the designmatrix, we do not know how long the data actually are. We only assume that event-latencies are synchronized with the data
"""
function equalize_size(
    X::AbstractMatrix,
    data::AbstractArray{T,2},
) where {T<:Union{Missing,<:Number}}
    @debug("2d equalize_size")
    if size(X, 1) > size(data, 2)
        X = @view X[1:size(data, 2), :]
    else
        data = @view data[:, 1:size(X, 1)]
    end
    return X, data
end
function equalize_size(
    X::AbstractMatrix,
    data::AbstractVector{T},
) where {T<:Union{Missing,<:Number}}
    @debug("1d equalize_size")
    if size(X, 1) > length(data)
        X = @view X[1:length(data), :]
    else
        data = @view data[1:size(X, 1)]
    end
    return X, data
end

function equalize_size(
    X::AbstractMatrix,
    data::AbstractArray{T,3},
) where {T<:Union{Missing,<:Number}}
    @debug("3d equalize_size")

    @assert size(X, 1) == size(data, 3) "Your events are not of the same size as your last dimension of data"

    return X, data
end



function clean_data(
    data::AbstractArray{T,2},
    winrej::AbstractArray{<:Number,2},
) where {T<:Union{Float64,Missing}}
    data = Array{Union{Float64,Missing}}(data)
    for row = 1:size(winrej, 1)
        data[:, Int.(winrej[row, 1]:winrej[row, 2])] .= missing
    end
    return data
end


macro maybe_threads(multithreading, code)
    return esc(:(
        if multithreading
            Threads.@threads($code)
        else
            $code
        end
    ))
end


poolArray(x) = PooledArray(x; compress = true)