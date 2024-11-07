module UnfoldCUDAExt
using CUDA
using Unfold
using CUDA.CUSPARSE
using LinearAlgebra

function Unfold.prepare(X, data::CuArray{T,2}) where {T<:Number}
    X_out = CuSparseMatrixCSC{T}(X)
    Y = CuArray{T}(data') # ch last!
    Ĥ = CUDA.zeros(T, size(Y, 2), size(X, 2))
    return Ĥ, Y, (X_out,)
end
prepare_XTX(Ĥ, data::Adjoint{T,CuArray}, X) where {T} = prepare_XTX(Ĥ, CuArray(data), X)

function Unfold.prepare(X, data::CuArray{T,3}) where {T<:Number}
    X_out = CuSparseMatrixCSC{T}(X)
    Y = CuArray{T}(permutedims(data, [3, 1, 2])) # trials - ch - time
    Ĥ = CUDA.zeros(T, size(Y, 2), size(Y, 3), size(X, 2))
    return Ĥ, Y, (X_out,)
end


end