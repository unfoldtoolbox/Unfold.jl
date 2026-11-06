"""
    solver_cv(; rng=MersenneTwister(), n_folds=5, shuffle=true, inner_solver=(X,y)->Unfold.solver_default(X,y), fit_test=false, breaks=nothing, fold_assignment=nothing)

A cross-validation solver for Unfold. This function generates a solver that can be passed to `Unfold.fit` to perform k-fold cross-validation on the data.

For **mass-univariate** models (epoched data), the observations are the trials: they are (optionally) shuffled and split into `n_folds` folds, and the model is fitted on the training trials for each fold.

For **deconvolved** (time-expanded / continuous-time) models, responses to single events overlap in time, so the observations are contiguous *segments* of the time course instead of single samples. The segments are either defined by user-provided `breaks` (e.g. the onsets of trials or blocks) or, if `breaks=nothing`, detected automatically as the maximal stretches of the design matrix that contain non-zero rows (auto-detection requires an always-on term, such as an intercept). For each fold, the model is fitted on all samples **outside** the held-out segments, so that responses to held-out events cannot leak into the training data.

# Arguments
- `rng`: A random number generator (e.g., `MersenneTwister(2)`) to use for shuffling the data if `shuffle=true` (default: `MersenneTwister()`! if you do not fix your seed manually, you will get a different shuffle each time).
- `n_folds::Int`: Number of folds for cross-validation (default: 5).
- `shuffle::Bool`: Whether to shuffle the data before splitting into folds (default: true).
- `inner_solver::Function`: A function that takes a design matrix `X` and data `y` and returns a model fit. This is the solver that will be used for each fold (default: `Unfold.solver_default`).
- `fit_test::Bool`: Whether to fit the model on the test set as well and return those estimates (default: false).
- `breaks`: (deconvolved models only) onsets (in samples) that mark the start of a new segment, e.g. the onsets of trials or blocks: a vector of numbers or a `DataFrame` with a `:latency` column (default: `nothing` = auto-detect segments from the design matrix).
- `fold_assignment`: (deconvolved models only) an explicit list of the segment indices to hold out in each fold (one entry per fold, e.g. `[[1, 2], [3], [4, 5]]`). When provided, the segments are not split by `kfolds`; instead, fold `i` is fitted on all segments *except* those listed in `fold_assignment[i]`. The entries must form a partition of the segments `1:n_segments`. This is useful for delta-R² comparisons: fit a full model, then pass `Unfold.fold_assignment(mf_full)` to a reduced model so both are evaluated on identical held-out segments. (default: `nothing` = derive folds from the segments).
# Returns
A function that can be used as a solver in `Unfold.fit`, which performs k-fold cross-validation and returns a `LinearModelFitCV` object containing the estimates, standard errors, and fold indices.
"""
function solver_cv(;
    rng = MersenneTwister(),
    n_folds = 5,
    shuffle = true,
    inner_solver = (X, y) -> Unfold.solver_default(X, y),
    fit_test = false,
    breaks = nothing,
    fold_assignment = nothing,
)

    # constructing the solver-function called by Unfold's fit()
    function cv_kernel(X, y)
        if ndims(y) == 3
            return _cv_mass_univariate(
                X, y; n_folds, shuffle, rng, inner_solver, fit_test,
            )
        elseif ndims(y) == 2
            return _cv_deconvolved(
                X, y; n_folds, shuffle, rng, inner_solver, fit_test, breaks,
                fold_assignment,
            )
        else
            error(
                "solver_cv: unsupported data dimensionality $(ndims(y)) - expected 2 (deconvolved model) or 3 (mass-univariate model)",
            )
        end
    end

    return cv_kernel
end


# ---------------------------------------------------------------------------
# mass-univariate: observations are the trials (rows of X / 3rd dimension of y)
# ---------------------------------------------------------------------------
function _cv_mass_univariate(X, y; n_folds, shuffle, rng, inner_solver, fit_test)
    n_trials = size(y, 3)
    indices = 1:n_trials

    # Generate split indices
    folds_iterator = kfolds(shuffle ? shuffleobs(rng, indices) : indices, k = n_folds)

    train_results = []
    test_results = [] # Only used if fit_test=true
    fold_indices = []

    for (i, (train_idx, test_idx)) in enumerate(folds_iterator)
        # Split the data
        res_train = inner_solver(X[train_idx, :], y[:, :, train_idx])
        push!(train_results, res_train)

        if fit_test
            res_test = inner_solver(X[test_idx, :], y[:, :, test_idx])
            push!(test_results, res_test)
        end
        push!(fold_indices, (train = train_idx, test = test_idx))
    end

    # 3. Concatenate results along a new 4th dimension (the folds)
    # fold_results[i].estimate is (chan, time, coeff)
    all_estimates = cat([r.estimate for r in train_results]..., dims = 4)
    all_stderror = cat([r.standarderror for r in train_results]..., dims = 4)

    # Prepare Info field
    info_dict = (
        inner_infos = [r.info for r in train_results],
        test_estimates = fit_test ?
                         cat([r.estimate for r in test_results]..., dims = 4) :
                         nothing,
    )
    return LinearModelFitCV(
        estimate = all_estimates,
        info = info_dict,
        standarderror = all_stderror,
        folds = fold_indices,
    )
end


# ---------------------------------------------------------------------------
# deconvolved: observations are contiguous segments of the time course
# ---------------------------------------------------------------------------
"""
    _cv_deconvolved(X, y; n_folds, shuffle, rng, inner_solver, fit_test, breaks, fold_assignment)

k-fold cross-validation for deconvolved (continuous-time) models.

The time course is split into segments (at the user-provided `breaks`, or at the
all-zero gaps of the design matrix), the segments are partitioned into `n_folds`
folds and the model is fitted on all samples outside the held-out segments. If
`fold_assignment` is provided, the partition is taken from it instead of from
`kfolds` (each entry lists the segments held out in that fold).
"""
function _cv_deconvolved(
    X, y; n_folds, shuffle, rng, inner_solver, fit_test, breaks, fold_assignment,
)
    n = size(y, 2)
    @assert size(X, 1) == n "solver_cv: design matrix ($(size(X, 1)) rows) and data ($(n) samples) have different lengths"

    segments = isnothing(breaks) ? _cv_segments_autodetect(X) :
                _cv_segments_manual(breaks, n)
    n_segments = length(segments)

    # (train_segments, test_segments) per fold
    fold_split = isnothing(fold_assignment) ?
                 _kfolds_split(n_folds, n_segments; shuffle, rng) :
                 _fold_split_from_assignment(fold_assignment, n_folds, n_segments)

    # which fold (test set) each segment belongs to
    seg_fold = zeros(Int, n_segments)

    train_results = []
    test_results = [] # Only used if fit_test=true
    fold_indices = []

    for (i, (train_seg, test_seg)) in enumerate(fold_split)
        for s in test_seg
            seg_fold[s] = i
        end

        # test rows: all samples inside the held-out segments,
        # fit rows: everything else
        in_test = falses(n)
        for s in test_seg
            in_test[segments[s]] .= true
        end
        test_rows = findall(in_test)
        fit_rows = findall(.!in_test)

        res_train = inner_solver(X[fit_rows, :], y[:, fit_rows])
        push!(train_results, res_train)

        if fit_test
            res_test = inner_solver(X[test_rows, :], y[:, test_rows])
            push!(test_results, res_test)
        end
        push!(fold_indices, (train = train_seg, test = test_seg))
    end

    # Leakage check (only meaningful for user-provided breaks): a non-zero row at
    # the boundary between two segments that belong to different folds means a
    # response is still active at the break and can leak across the
    # train/test boundary. For auto-detected segments the boundaries are
    # all-zero rows of the design matrix, so no response can cross them.
    if !isnothing(breaks)
        warned = Int[]
        for i in 1:(n_segments - 1)
            if seg_fold[i] != seg_fold[i + 1] && sum(abs, X[last(segments[i]), :]) > 0
                push!(warned, last(segments[i]))
            end
        end
        if !isempty(warned)
            @warn "cross-validation boundary at row(s) $warned is non-zero in the design matrix: a response is still active at the break and may leak across the train/test boundary. Consider placing breaks in quiet intervals."
        end
    end

    # Concatenate results along a new 3rd dimension (the folds)
    # train_results[i].estimate is (chan, coeff)
    all_estimates = cat([r.estimate for r in train_results]..., dims = 3)
    all_stderror = cat([r.standarderror for r in train_results]..., dims = 3)

    # Prepare Info field
    info_dict = (
        inner_infos = [r.info for r in train_results],
        segments = segments,
        test_estimates = fit_test ?
                         cat([r.estimate for r in test_results]..., dims = 3) :
                         nothing,
    )
    return LinearModelFitCV(
        estimate = all_estimates,
        info = info_dict,
        standarderror = all_stderror,
        folds = fold_indices,
    )
end


# (train_segments, test_segments) per fold: the standard k-fold partition of the
# segments
function _kfolds_split(n_folds, n_segments; shuffle, rng)
    n_folds <= n_segments || throw(
        ArgumentError(
            "solver_cv: n_folds ($n_folds) exceeds the number of segments ($n_segments) - provide more `breaks` or use fewer folds",
        ),
    )
    return [
        (collect(train_seg), collect(test_seg))
        for (train_seg, test_seg) in
           kfolds(shuffle ? shuffleobs(rng, 1:n_segments) : 1:n_segments, k = n_folds)
    ]
end


# a fixed split: one entry per fold, each a list of the segment indices held out
# in that fold; the entries must partition the segments 1:n_segments. This is
# how the split of a previous fit (e.g. via fold_assignment(modelfit)) can be
# reused in another fit, so both are evaluated on identical held-out segments.
function _fold_split_from_assignment(fold_assignment, n_folds, n_segments)
    length(fold_assignment) == n_folds || throw(
        ArgumentError(
            "solver_cv: fold_assignment has $(length(fold_assignment)) entries but n_folds = $n_folds",
        ),
    )
    sort(vcat(fold_assignment...)) == collect(1:n_segments) || throw(
        ArgumentError(
            "solver_cv: fold_assignment must hold out every segment (1:$n_segments) exactly once - for comparable fits, use the same `breaks`",
        ),
    )
    return [(setdiff(1:n_segments, te), collect(te)) for te in fold_assignment]
end


# segments defined by user-provided breaks: each break onset starts a new segment;
# the prefix before the first break and the suffix after the last break are
# segments as well
function _cv_segments_manual(breaks, n)
    b = breaks isa DataFrame ? breaks[!, :latency] : breaks
    b = sort(unique(Int[round(bi) for bi in b]))
    b = filter(bi -> 1 < bi <= n, b)
    segments = UnitRange{Int}[]
    prev = 1
    for bi in b
        push!(segments, prev:(bi - 1))
        prev = bi
    end
    push!(segments, prev:n)
    return segments
end


# segments detected automatically: maximal stretches of rows with non-zero design
# matrix values. Requires an always-on term (e.g. an intercept) so that event
# windows are non-zero throughout their support.
function _cv_segments_autodetect(X)
    n = size(X, 1)
    active = vec(sum(abs.(X), dims = 2)) .> 0
    any(active) || throw(
        ArgumentError(
            "solver_cv: the design matrix has no non-zero rows - cannot auto-detect segments. Provide `breaks` manually.",
        ),
    )
    segments = UnitRange{Int}[]
    i = 1
    while i <= n
        if active[i]
            j = i
            while j < n && active[j + 1]
                j += 1
            end
            push!(segments, i:j)
            i = j + 1
        else
            i += 1
        end
    end
    return segments
end


# the fold dimension is the last (Nth) dimension of the estimate
coef(m::LinearModelFitCV{T,N}) where {T,N} = dropdims(mean(m.estimate, dims = N), dims = N)

# the observations held out per fold (trials for mass-univariate fits, segment
# indices for deconvolved fits); pass the result of a previous fit as the
# `fold_assignment` keyword to evaluate another model on the very same held-out
# observations (e.g. a reduced model for a delta-R² comparison)
fold_assignment(mf::LinearModelFitCV) = [collect(f.test) for f in mf.folds]

function stderror(m::LinearModelFitCV{T,N}) where {T,N}
    if isempty(m.standarderror)
        stderror = fill(nothing, size(m.estimate)[1:N-1])
    else
        stderror = dropdims(mean(T.(m.standarderror), dims = N), dims = N)
    end
end
