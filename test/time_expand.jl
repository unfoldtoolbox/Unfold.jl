tbl = DataFrame([collect(1:4:1000)], [:latency])
X = ones(size(tbl))


basisfunction = firbasis(τ = (-3, 0), sfreq = 1, name = "testing")
term = Unfold.TimeExpandedTerm(FormulaTerm(Term, Term), basisfunction, :latency);

Xdc = Unfold.time_expand(X, term, tbl)

kernel = Unfold.kernel
ncolsBasis = size(kernel(term.basisfunction)(0.), 2)
X = reshape(X, size(X, 1), :)

ncolsX = size(X)[2]
nrowsX = size(X)[1]
ncolsXdc = ncolsBasis * ncolsX
onsets = tbl[!, term.eventfields[1]]

if typeof(term.eventfields) <: Array && length(term.eventfields) == 1
    bases = kernel(term.basisfunction).(tbl[!, term.eventfields[1]])
else
    bases = kernel(term.basisfunction).(eachrow(tbl[!, term.eventfields]))
end

rows = Unfold.timeexpand_rows(onsets,bases,Unfold.shiftOnset(term.basisfunction),ncolsX)
cols = Unfold.timeexpand_cols(term,bases,ncolsBasis,ncolsX)
vals = Unfold.timeexpand_vals(bases,X,size(cols),ncolsX)

@test Unfold.timeexpand_cols_allsamecols(bases,ncolsBasis,ncolsX) == Unfold.timeexpand_cols_generic(bases,ncolsBasis,ncolsX)
