tbl = DataFrame([collect(1:4:1000)], [:latency])
X = ones(size(tbl))

@testset "timeexpand internals" begin
    basisfunction = firbasis(τ = (-3, 0), sfreq = 1, name = "testing")
    term = Unfold.TimeExpandedTerm(FormulaTerm(Term, Term), basisfunction, :latency)

    Xdc = Unfold.time_expand(X, term, tbl)

    kernel = Unfold.kernel
    ncolsBasis = Unfold.width(term.basisfunction)
    ncolsX = size(X, 2)
    onsets = tbl[!, term.eventfields[1]]
    bases = kernel.(Ref(term.basisfunction), onsets)

    rows = Unfold.timeexpand_rows(
        onsets,
        bases,
        Unfold.shift_onset(term.basisfunction),
        ncolsX,
    )
    cols = Unfold.timeexpand_cols(term, bases, ncolsBasis, ncolsX)
    vals = Unfold.timeexpand_vals(bases, X, size(rows), ncolsX)

    @test Unfold.timeexpand_cols_allsamecols(bases, ncolsBasis, ncolsX) ==
          Unfold.timeexpand_cols_generic(bases, ncolsBasis, ncolsX)
end
