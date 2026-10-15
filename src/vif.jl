"""
    vif(m::UnfoldModel)

Variance Inflation Factor (VIF) of a model's designmatrix, to check whether the model
suffers from (multi-)collinearity.

The VIF of a column is the corresponding diagonal element of the inverse of the
correlation matrix of the designmatrix columns, i.e. `VIF_j = 1 / (1 - R^2_j)`, where
`R^2_j` is the R^2 from regressing column j on all other columns. `VIF = 1` means no
collinearity; large values (rule of thumb: `> 10`) indicates that a combination of other
predictors try to explain the same variability. Colloquial this means the SE are becoming
larger, as the model does cannot uniquely distribute variance.

# Arguments
- `m::UnfoldModel`: Only the designmatrix is used (`modelmatrix(m)`)

# Returns
a `DataFrame` with columns `:coefname` (the coefficient names,
e.g. `"stim : condition : 0.1"`) and `:VIF`, plus an `:eventname` column for
mass-univariate models.

Constant columns (e.g. the intercept) have no defined VIF and are reported as `Inf`.
The same applies to perfectly collinear columns. High VIFs are to be expected for
splines or overlapping basisfunctions and are usually not a problem, as long as
the solver converges.

`vif` is a [`StatsModels`](https://juliastats.org/StatsModels.jl/stable/) function, which Unfold extends for `UnfoldModel`s and re-exports.

# Examples
```julia-repl
julia> using UnfoldSim
julia> data, evts = UnfoldSim.predef_eeg()
julia> f = @formula 0 ~ 1 + condition
julia> m = fit(UnfoldModel, [Any => (f, firbasis(τ = (-1, 1), sfreq = 100))], evts, data)
julia> v = vif(m)
julia> worst = v[v.VIF .> 10, :]
0×2 DataFrame
 Row │ term    VIF
     │ String  Float64
─────┴─────────────────
```
"""
function vif(m::UnfoldModel)
    X = modelmatrix(m)
    if X isa AbstractMatrix
        names = get_coefnames(m)
        @assert length(names) == size(X, 2) "coefnames ($(length(names))) do not match the designmatrix ($(size(X, 2)))"
        return DataFrame(:coefname => String.(names), :VIF => _vif(X))
    end
    if X isa Vector
        keys = first.(design(m))
        eventnames = String[]
        terms = String[]
        vals = Float64[]
        for (k, key) in enumerate(keys)
            names = StatsModels.coefnames(formulas(m)[k].rhs)
            @assert length(names) == size(X[k], 2) "coefnames ($(length(names))) do not match the designmatrix ($(size(X[k], 2)))"
            for nm in names
                push!(eventnames, string(key))
                push!(terms, string(nm))
            end
            vals = vcat(vals, _vif(X[k]))
        end
        return DataFrame(:eventname => eventnames, :coefname => terms, :VIF => vals)
    end
    error("vif is not supported for a modelmatrix of type $(typeof(X))")
end


# VIF of a designmatrix: `1 / (1 - R^2_j)`, i.e. the diagonal of the inverse of the
# correlation matrix of the columns. The column sums and the Gram matrix X'X (at most
# n x n) are computed directly, so sparse designmatrices are never densified.
function _vif(X::AbstractMatrix)
    nr, n = size(X)
    if n > 20_000
        @warn "vif needs an n x n correlation matrix, ~$(round(n^2 * 3 / 1024^3; digits = 1)) GB for your $n-column designmatrix"
    end
    S = Matrix(X' * X) # potentially this could be improved for sparse matrices?
    m = vec(sum(X, dims = 1)) ./ nr
    # centered covariance * nr
    C = S .- nr .* m * m'
    d = sqrt.(max.(diag(C), 0.0))
    # constant columns (e.g. intercept) have no defined VIF
    keep = d .> (sqrt(eps(Float64)) * maximum(d))
    orig = findall(keep)
    vif = fill(Inf, n)
    isempty(orig) && return vif
    R = C[keep, keep] ./ (d[orig] * d[orig]')
    try
        F = cholesky(Hermitian(R))
        vif[orig] = vec(sum((F.L \ I(length(orig))) .^ 2, dims = 1))
    catch
        # R exactly singular: perfectly collinear columns. Rank-revealing QR
        # (R shares the rank of the columns of X) gives the independent subset;
        # the dependent columns keep VIF = Inf
        Fq = qr(R, ColumnNorm())
        Rd = abs.(diag(Fq.R))
        rnk = count(Rd .> (Rd[1] * max(size(Fq.R)...) * eps(Float64)))
        ind = Fq.p[1:rnk]
        Fi = cholesky(Hermitian(R[ind, ind]))
        vif[orig[ind]] = vec(sum((Fi.L \ I(rnk)) .^ 2, dims = 1))
    end
    return vif
end
