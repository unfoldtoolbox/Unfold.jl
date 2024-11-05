
X, y = Unfold.equalize_size(modelmatrix(m), data)
X = Unfold.SparseMatrixCSC(X) # get rid of view


#@b Unfold.solver_default(X,y)


@time Unfold.lsmr(X, @view(data[1, :]); log = true)


function scale_ls!(A)
    s = ones(size(A, 2))
    for j = 1:size(A, 2)
        i = A.colptr[j]
        k = A.colptr[j+1] - 1
        nj = i <= k ? norm(A.nzval[i:k]) : 0.0
        if nj > 0.0
            A.nzval[i:k] ./= nj
            s[j] = nj
        end
    end
    return s  # s contains the diagonal elements of the (right) diagonal scaling
end

@time begin
    s = scale_ls!(X)
    Unfold.lsmr(X, @view(data[1, :]); log = true)
end