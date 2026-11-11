"""
    $(SIGNATURES)
convert and permutedim input to follow the following output:

    Ĥ, Y, X = prepare(X, data)

where `Ĥ` is used to save the beta, `Y` is the data in format ch x repeat x time (with size(time) = 1 if data is a Matrix), and `X`.

- if data is a CuArray, everything is transformed to CuArrays as well (via UnfoldCUDAExt.jl, CUDA needs to be loaded)
- same datatype between X and data is enforced
"""
function prepare(X, data::AbstractArray{<:Union{Missing,T},2}) where {T<:Number}
    @info "Missings in data - we remove any timepoint from data and designmatrix"
    ix = .!any.(ismissing, eachslice(data; dims = 2))
    return prepare(disallowmissing(@view(X[ix, :])), disallowmissing(@view(data[:, ix])))##Array{T,2}(@view(data[:, ix])))
end

function prepare(X, data::AbstractArray{<:Union{Missing,T},3}) where {T<:Number}
    @info "Missings in data - we remove any trial from data and designmatrix"
    ix = .!any.(ismissing, eachslice(data; dims = 3))
    return prepare(disallowmissing(@view(X[ix, :])), disallowmissing(@view(data[:, :, ix])))#Array{T,3}(@view(data[:, :, ix])))
end

function prepare(X, data::AbstractArray{<:Number,3})
    Ĥ = zeros(eltype(data), size(data, 1), size(X, 2), size(data, 2))
    return Ĥ, permutedims(data, [3, 1, 2]), (X,)
end

function prepare(X, data::AbstractArray{<:Number,2})
    Ĥ = zeros(eltype(data), size(data, 1), size(X, 2))
    return Ĥ, (data'), (X,)
end

#---
# bychannel: `Xs` is a Vector with one designmatrix per channel (length == size(data, 1))
function prepare(
    Xs::AbstractVector{<:AbstractMatrix},
    data::AbstractArray{<:Union{Missing,T},3},
) where {T<:Number}
    @info "Missings in data - we remove any trial from data and designmatrix"
    ix = .!any.(ismissing, eachslice(data; dims = 3))
    # each channel's designmatrix has its own rows - subset them by the non-missing trials
    return prepare(
        disallowmissing.([@view(X[ix, :]) for X in Xs]),
        disallowmissing(@view(data[:, :, ix])),
    )
end

function prepare(Xs::AbstractVector{<:AbstractMatrix}, data::AbstractArray{<:Number,3})
    @assert length(Xs) == size(data, 1) "bychannel designmatrix has $(length(Xs)) channels, but data has $(size(data,1)) - did you maybe order the events/baseline values differently than the data channels?"
    Ĥ = zeros(eltype(data), size(data, 1), size(Xs[1], 2), size(data, 2))
    # one state-tuple per channel, splatted per channel in `solver_main`
    return Ĥ, permutedims(data, [3, 1, 2]), [(X,) for X in Xs]
end

"""
    $(SIGNATURES)
instead of solving y = Xb, we solve X'Xb = X'y. This function calculates X'X and instantiates X'y to be used in the solver-step, to facilitate X'y calculations later, X' is also calculated.

"""
prepare_XTX(all::Tuple) = prepare_XTX(all...)
prepare_XTX(Ĥ, data::AbstractArray, all::Tuple) = prepare_XTX(Ĥ, data, all...)
prepare_XTX(Ĥ, data::Adjoint, X::Tuple) = prepare_XTX(Ĥ, collect(data), X)
function prepare_XTX(Ĥ, data::T1, X::T2) where {T1,T2}
    Xt = X'
    R_xx = T1(Xt * X)
    R_xy = similar(data, size(X, 2))
    return Ĥ, data, (Xt, R_xx, R_xy)
end

#---
# bychannel: the prepared state is a Vector with one state-tuple per channel, selected per
# channel in `solver_main` via `_channel_prepared`.

function prepare_XTX(Ĥ, data, ps::AbstractVector{<:Tuple})
    # R_xx per channel as a plain Matrix (the shared 3D path instantiates it via the
    # data-type, which fails for 3D data)
    return Ĥ, data, [(X', X' * X, similar(data, size(X, 2))) for X in first.(ps)]
end
prepare_pinv(Ĥ, data::AbstractArray{T,3}, ps::AbstractVector{<:Tuple}) where {T} =
    Ĥ, data, last.(prepare_pinv.(Ref(Ĥ), Ref(data), ps))
prepare_qr(Ĥ, data, ps::AbstractVector{<:Tuple}) =
    Ĥ, data, last.(prepare_qr.(Ref(Ĥ), Ref(data), ps))
prepare_cholesky(Ĥ, data, ps::AbstractVector{<:Tuple}) =
    Ĥ, data, last.(prepare_cholesky.(Ref(Ĥ), Ref(data), ps))
