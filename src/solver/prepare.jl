"""
    $(SIGNATURES)
convert and permutedim input to follow the following output:

    Ĥ, Y, X = prepare(X, data)

where `Ĥ` is used to save the beta, `Y` is the data in format ch x repeat x time (with size(time) = 1 if data is a Matrix), and `X`.

- if data is a CuArray, everything is transformed to CuArrays as well (via UnfoldCUDAExt.jl, CUDA needs to be loaded)
- same datatype between X and data is enforced
"""
function prepare(X, data::AbstractArray{<:Union{Missing,T},2}) where {T<:Number}
    @warn "Missings in data - we remove any timepoint from data and designmatrix"
    ix = .!any.(ismissing, eachslice(data; dims = 2))
    return prepare(disallowmissing(@view(X[ix, :])), disallowmissing(@view(data[:, ix])))##Array{T,2}(@view(data[:, ix])))
end

function prepare(X, data::AbstractArray{<:Union{Missing,T},3}) where {T<:Number}
    @warn "Missings in data - we remove any trial from data and designmatrix"
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
