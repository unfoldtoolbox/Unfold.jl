# Dense time-expanded design matrix: multiple-dispatch integration report

Branch `feature-denseDesign` (worktree of `main` @ v0.8.13).
Goal: a new `<:UnfoldModel` that builds the time-expanded (deconvolution) design matrix as a
dense `Matrix` instead of the `SparseMatrixCSC` hardcoded in
`DesignMatrixLinearModelContinuousTime`, flowing through the existing
fit / designmatrix / solver / predict machinery by multiple dispatch, with minimal surgical
relaxations of the spots that block new model types.

## What was added

| Where | What |
|---|---|
| `src/denseDesign.jl` | *All* DenseDesign-specific code in one file (~75 lines): the `DesignMatrixLinearModelContinuousTimeDense{T}` struct (like `DesignMatrixLinearModelContinuousTime` but `modelmatrix::Matrix{T}`; the sparse time-expanded matrix produced by the existing `time_expand` machinery is converted on construction — a plain struct-field conversion, no code path change), the `UnfoldLinearModelContinuousTimeDenseDesign{T} <: UnfoldModel{T}` model with the same 5-constructor chain as `UnfoldLinearModelContinuousTime`, the `@traitimpl ContinuousTimeTrait{...}` (one line makes *all* `@traitfn`-based machinery work), and one `modelmatrix` routing method. Everything else (design, fit, solver, predict, effects, io) works via the trait or generic `AbstractMatrix` dispatch. |
| `src/Unfold.jl` | `include("denseDesign.jl")` + export of the new model type. |
| `test/densedesign.jl` | new tests: construction & type checks, sparse↔dense estimate equivalence (single- and multi-event), user-facing interface (coeftable, predict, residuals/r², effects, modelmatrix-without-basisfunction, times, manual `designmatrix!`/`fit!` pathway), end-to-end on `UnfoldSim.predef_eeg` with the known r²≈0.74 value measured with the sparse model, save/load roundtrip. |

## Inventory: hardcoded spots that had to be touched

1. **`src/designmatrix.jl` — `extend_to_larger` keys on the matrix *element type*, not the model type.**
   `extend_to_larger(::Vector{<:SparseMatrixCSC})` splats and pairwise-concatenates; the catch-all
   `extend_to_larger(::Vector{<:AbstractMatrix}) = modelmatrix` keeps the vector — which is the
   correct behavior for the *mass-univariate* modelmatrix (a vector of per-event matrices, used
   element-wise by `fit!`, `coeftable`, `get_basis_names`, `predict`). This is the single hardest
   blocker: a dense time-expanded modelmatrix is also a `Vector{<:Matrix}`, so dispatching
   "concatenate vs. keep vector" on the matrix element type *conflates* "mass-univariate" with
   "dense". **Resolution:** the pairwise sparse methods were generalized to
   `extend_to_larger(::AbstractMatrix, ::AbstractMatrix[, args...])` (per the code review), with a
   small `padzeros` helper dispatching the padding (sparse: reuse the stored `colptr`/`rowval`
   structure, no copy; dense: `vcat` zero rows) — one varargs-reduction and one base method serve
   both storages. The routing, however, cannot key on the matrix element type: it is moved up one
   level to dispatch on the *designmatrix* type (in `src/denseDesign.jl`):
   `StatsModels.modelmatrix(d::Vector{<:DesignMatrixLinearModelContinuousTimeDense}) =
   extend_to_larger(modelmatrices.(d)...)`. The generic `Vector{<:AbstractMatrix}` method is
   untouched, so mass-univariate behavior is provably unchanged.

2. **`src/designmatrix.jl` — `modelmatrix(uf, basisfunction)` was bound to the concrete sparse type.**
   Called by `effects` (`modelmatrix(model, false)`) to get the design matrix without the
   time-expansion. Per the code review, the `Union{...}` of CT model types was replaced with a
   `@traitfn ... where {T<:UnfoldModel;ContinuousTimeTrait{T}}` — future CT models are covered by
   the trait instead of being added to a `Union`. The default `= true` was dropped: it created a
   1-arg method competing (definition-order-dependently) with the generic
   `modelmatrix(uf::UnfoldModel)`; one-argument calls now always go through the generic (identical
   behavior — its body is the same `modelmatrix(designmatrix(uf))`), two-argument calls through
   the traitfn.

3. **`src/fit.jl` — two `isa(uf, UnfoldLinearModel)` checks in `fit!`.**
   Replaced with `!SimpleTraits.istrait(ContinuousTimeTrait{typeof(uf)})` — the same trait check
   already used a few lines above for input validation. For every type reachable today the branch
   decisions are identical (`UnfoldLinearModel` is not CT → branch taken, as with `isa`;
   `UnfoldLinearModelContinuousTime` is CT → branch skipped, as with `isa`). The mass-univariate
   per-event fitting loop and the times-vs-data assertion are now trait-driven, so any future
   non-CT model flows through them and any future CT model flows around them.

4. **`src/fit.jl` — the `fit!` split point itself.** Not touched, but worth recording: `fit!`
   has exactly two behavior paths (mass-univariate per-event loop vs. one `solver(X, data)` call
   after `equalize_size`). The CT path is fully storage-agnostic: `solver_default` →
   `solver_main` → `prepare`/`prepare_XTX`/`solver_lsmr!` accept any `AbstractMatrix`, so the
   dense model needed **zero** solver changes. The `SparseMatrixCSC` specializations in
   `src/solver.jl` are dead code (the file is commented out of the build in `Unfold.jl`).

## Inventory: spots that already worked untouched

- **Design-matrix construction:** `designmatrix(type, f, tbl, basisfunction)` builds the sparse
  `Xdc` via `StatsModels.modelcols(::TimeExpandedTerm, tbl)` → `time_expand`. The dense model
  receives its dense storage purely through struct-field conversion — no branch, no dispatch
  change. The FIR-basis machinery (`firbasis`, `firkernel`, fast `time_expand_firdiag`,
  interpolate/scale_duration paths) is reused as-is.
- **All `@traitfn` machinery** (`ContinuousTimeTrait`): `check_data`, `coeftable`, `times`,
  `get_basis_names`, `residuals`, `predict_no_overlap`, `predict_partial_overlap`,
  `effects`/`_typify` — covered by the single `@traitimpl`.
- **`predict.jl`:** `predict(X::AbstractMatrix, coefs::AbstractMatrix) = (X*coefs')'`,
  `matrix_by_basisname`, and `non_zero_rows` (already has an `AbstractArray` fallback) — generic.
- **`solver/` (main, solvers, prepare, solver_cv):** no matrix-type hardcoding in the live code.
- **`io.jl` save/load:** generic via `typeof(uf.designmatrix[k])(...)` /
  `typeof(designmatrix(uf)[1]).(form, X, events)`; the loaded dense model reconstructs to dense
  by the same field-conversion. `empty_modelmatrix` works via the new empty constructors.
- **`vif.jl`, `show.jl`, `condense.jl`:** generic (size-based or trait-based).

## Inventory: spots that resisted conversion (left hardcoded, on purpose)

1. **`design_to_modeltype(f, t::BasisFunction) = UnfoldLinearModelContinuousTime` (`src/fit.jl`).**
   Autodetection always returns the sparse CT model. Converting would require a *policy* (which
   storage should autodetection choose?), i.e. a user-visible decision — a keyword such as
   `fit(...; dense_design = true)` or a `design_to_modeltype` hook per basisfunction. For now the
   dense model is requested explicitly, exactly like `UnfoldLinearModel` vs.
   `UnfoldLinearModelContinuousTime` differ today. *Relaxation candidate for a follow-up.*

2. **`ContinuousTimeTrait` is a boolean trait.** Two models with identical *dynamics* but
   different *storage* need one `@traitimpl` each, and every `@traitfn` is duplicated per trait
   side. A richer trait (e.g. `ContinuousTimeTrait` + a `StorageStyle{Sparse/Dense}` axis) would
   remove the duplication — but touching every `@traitfn` for that is beyond "surgical" and was
   not done. *Relaxation candidate if a third storage appears.*

3. **Constructor boilerplate per model type (`src/typedefinitions.jl`).** Every new
   `<:UnfoldModel` repeats ~5 constructors + an empty `DesignMatrix` struct. A parametric
   `UnfoldLinearModelContinuousTime{T,S<:Storage}` would unify this but changes existing public
   types (rejected for this branch; see verdict below).

4. **`io.jl` line 73: `generate_Xs || uf isa UnfoldLinearModel`.** An io-specific optimization
   (mass-univariate design matrices are small, so always regenerate). Left as-is: for the dense
   model `generate_Xs=false` leaves an empty designmatrix, exactly like the sparse CT model. Not
   a blocker; belongs to the same "key on trait, not concrete type" cleanup as fit!.

5. **`modelmatrix(uf, basisfunction)` `Union` widening** — *resolved* in the code review: the
   method is now a `@traitfn` on `ContinuousTimeTrait`, so future CT model types are covered
   without touching this signature.

## Verdict: how modular is `UnfoldModel` already?

**Quite modular along the fit axis, weak along the storage axis.**

- The abstract type + one boolean trait + `@traitfn` covers *dynamics* (epoched vs.
  continuous-time) cleanly: a new model needs one traitimpl and everything downstream
  (coeftable, effects, predict, residuals, r², times) works.
- The solver pipeline is already storage-agnostic end to end; nothing there blocks any
  `AbstractMatrix` design matrix.
- The weak spot is code keyed on the matrix *element type* instead of the model/design type —
  `extend_to_larger` being the clearest case, where "sparse ⇒ concatenate" and "dense ⇒ keep
  vector" were conflated. The fix applied here (dispatch on the designmatrix type at the
  `modelmatrix(d::Vector{<:AbstractDesignMatrix})` boundary) is the pattern to reuse.
- The remaining friction is boilerplate, not blockers: each new model repeats constructors and
  needs one-line additions to 1–2 Unions/signatures.

**What should be relaxed/unified for future model types (recommended, not implemented here):**
1. Route *all* storage-vs-vector decisions through the designmatrix type. The pairwise
   `extend_to_larger` methods are now storage-shared (`AbstractMatrix` + `padzeros`); the
   remaining element-type-dispatched behavior after this branch is only the
   `Vector{<:SparseMatrixCSC}` splat-routing, which the sparse CT model could move to the same
   designmatrix-type pattern for symmetry.
2. Replace the two remaining `isa(uf, UnfoldLinearModel)`-style checks (io.jl) with the trait.
3. Consider a parametric storage axis (`UnfoldLinearModelContinuousTime{T,Storage}`) or a
   `StorageStyle` trait if a third representation (e.g. banded/woodbury, or per-channel) appears;
   until then the boolean trait + explicit types is fine.
4. Give `design_to_modeltype` a user-facing hook (keyword) so autodetection can choose storage
   without hardcoding.

## Design notes: how full TRF with time-continuous predictors could plug in later

*Analysis only — not implemented on this branch.*

"Full TRF" here means a linear time-response model whose predictors are *continuous input
streams* (e.g. the speech envelope, pixel luminance) that are lagged/convolved against the
continuous EEG, rather than event-locked regressors from an events table.

- **Candidate A — build on `feature-bychannel` (`bychannel(sym)` / `ChannelVaryingTerm`).**
  That branch already introduces per-channel predictors: an event-table column holding one value
  per channel per row, producing per-channel design matrices and a per-channel fit. It is
  currently mass-univariate-only and `predict` is unimplemented. For time-continuous TRF with
  per-channel predictors, the natural extension is a time-expanded `ChannelVaryingTerm` (or
  `bychannel` applied *inside* the basis-function kernel): `modelcols` would then yield a
  per-channel `Xdc` (nch × time × lags), one dense `Matrix` per channel — precisely the dense
  storage this branch introduces. The blockers are known and now partly removed by this branch:
  `fit!` is trait-driven (a per-channel TRF needs its own `fit!` dispatch, like the
  mass-univariate bychannel fit), `extend_to_larger`/`modelmatrix` dispatch on the designmatrix
  type (a per-channel designmatrix type would get its own method), and `time_expand` needs a
  variant that accepts channel-indexed regressors. Sparse storage would make per-channel
  slicing awkward; dense (this branch) keeps per-channel views trivial — so the dense model is
  the natural base for A.

- **Candidate B — stream-based interface.** Let the basis function be evaluated at *every*
  sample of a predictor stream instead of per event: a new `BasisFunction` subtype (e.g.
  `LagBasis`) whose "events" are samples, an `apply_basisfunction` variant that wraps a
  stream term into `TimeExpandedTerm`-like structure, and a `modelcols` method that builds the
  lag matrix directly (dense, ideally banded). The events table would shrink to metadata
  (onsets/segments) or disappear into `eventfields`. Insertion points: `designmatrix(...)` for
  accepting a predictor `Matrix` alongside `evts`, and a new `predict` variant that consumes a
  stream instead of an events table (`predict(uf, f, evts)` is events-centric today). The
  solver, `fit!` (post-relaxation), `coeftable` and `r2` would work unchanged if the term
  yields the same dense `Xdc` layout the dense model consumes — which is the main argument for
  keeping the concrete `Xdc` contract (`modelmatrix(uf)` = one concatenated matrix, columns =
  basis × predictor-lags) stable.

- **Recommendation:** whichever candidate wins, define the TRF term so `modelcols` emits the
  same layout the dense model already consumes; then fit/solve/predict/effects come for free.
  Candidate A is the smaller step (reuses the events-table interface and the existing
  per-channel work); candidate B is the cleaner long-term interface but needs new
  `BasisFunction` + `predict` machinery. An alternative interface worth considering is
  `UnfoldSim`-style predictor matrices with explicit `sfreq` metadata, making the stream
  first-class and the events table optional.

## Test results

Full suite via `Pkg.test()`, run after all changes: **all testsets pass** — BasisFunctions 24,
Fitting 39, Designmatrix 18, **DenseDesign 28 (new)**, Splines 13, Predict 40, Effects 44,
VIF 28, Utilities 44, IO 16, solver_cv 23. That is the pre-existing 289 tests (identical to the
baseline on pristine `main`, verified before the change) plus the 28 new `DenseDesign` tests,
`exit: 0`. Sparse↔dense estimates agree within `rtol=1e-5, atol=1e-8` on identical
data with the same default solver; the dense model recovers the known simulated effects
(`UnfoldSim.predef_eeg`, r² ≈ 0.74, face > car condition effect).

## Deferred

- Cross-validation (`solver_cv`) and regularized solvers for the dense model (explicitly out of
  scope; the default `solver_default` path is shared and works).
- Autodetection of the dense model from the design (see inventory item 1).
- Full TRF implementation (analysis above only).
