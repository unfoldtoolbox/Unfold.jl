#  # [Cross-validated Unfold models](@id solver_cv)
# This tutorial shows how to run an Unfold model with k-fold cross-validation.
# # Setup

using Unfold
using UnfoldSim
using Random
using Statistics
using UnfoldMakie, CairoMakie

# Define formula and basis function for a mass-univariate model.
f = @formula 0 ~ 1 + condition

# # Mass-univariate models
# For epoched (mass-univariate) data, the cross-validation units are the trials: they are (optionally) shuffled and split into folds.
eeg, evts = UnfoldSim.predef_eeg(; return_epoched = true, noiselevel = 15)

# `solver_cv` wraps the default solver, but runs it for each cross-validation fold
cv_solver = solver_cv(n_folds = 5, shuffle = true) # shuffle is true by default

# Now we can fit the model with the CV solver.
m_cv = fit(UnfoldModel, f, evts, eeg, 1:size(eeg, 1); solver = cv_solver)

# The 4th dimension contains the CV-fold (channel, time, coefficient, fold).
size(modelfit(m_cv).estimate)

# You also get train/test indices for each fold
# (we have to index once into the first "event", e.g. you could run multiple events / formulas in one model)
length(m_cv.modelfit.folds[1])
# we can also access e.g. the third fold and check the train/test indices. Let's display only the first 6 indices
first(m_cv.modelfit.folds[1][3].train, 6)

# # `coef` and `coeftable``
# For `LinearModelFitCV`, `coef(m)` returns the mean over folds.
# This means `coeftable(m_cv)` reports the fold-averaged estimates.
first(coeftable(m_cv), 6)

# You can access fold-specific estimates directly from `modelfit.estimate`.
fold_1_estimate = modelfit(m_cv).estimate[:, :, :, 1]
size(fold_1_estimate)

# # Finally let's plot our estimates

f, ax, h = series(modelfit(m_cv).estimate[1, :, 1, :]')
lines!(coef(m_cv)[1, :, 1], color = :black, linestyle = :dash)
ax.xlabel = "Time (samples)"
f

# The colored lines are the fold-specific estimates, the dashed black line is the mean across folds.

# # Deconvolved models
# For deconvolved (time-expanded) models, the response to a single event lasts longer than the inter-stimulus interval, so individual samples cannot be cross-validated without leaking into each other.
# Instead, the cross-validation units are contiguous *segments* of the time course, and for each fold the model is fitted on all samples **outside** the held-out segments.
#
# The segments are either defined manually with `breaks` (the onsets of each new segment, e.g. the onsets of trials or blocks: a vector of numbers or a `DataFrame` with a `:latency` column) or, by default, detected automatically as the maximal stretches of the design matrix that contain non-zero rows (auto-detection requires an always-on term, such as an intercept).
# If a response is still active at a fold boundary, `solver_cv` warns about it.
eeg_cont, evts = UnfoldSim.predef_eeg(; return_epoched = false, noiselevel = 15)
bf = firbasis(τ = (-0.1, 0.15), sfreq = 100)
cv_solver = solver_cv(n_folds = 5) # segments are auto-detected from the design matrix
m_cv = fit(UnfoldModel, [Any => (f, bf)], evts, eeg_cont; solver = cv_solver)

# The 3rd dimension contains the CV-fold (channel, coefficient, fold).
size(modelfit(m_cv).estimate)

# `coef` and `coeftable` again report the fold-averaged estimates.
first(coeftable(m_cv), 6)

# # Reusing a fold split
# `fold_assignment` returns, for each fold, the observations (trials, or segments for deconvolved models) held out in that fold.
# Passing the result of a previous fit back as the `fold_assignment` keyword evaluates a different model on the *identical* held-out segments — useful for delta-R² comparisons between a full model and a reduced model (with the reduced model, auto-detection would not produce the same segments, so give both fits the same manual `breaks`).
fold_assignment(modelfit(m_cv))

# # Reusing a split: `fold_assignment`
# For model comparisons — e.g. a delta-R² analysis in which an effect of interest is removed from the model — the full and the reduced model need to be evaluated on the *same* held-out segments.
# Instead of constructing the split manually, extract it from the finished fit with `fold_assignment` (one list of held-out segments per fold) and pass it to the next fit:
f_red = @formula 0 ~ 1 # reduced model: the effect of interest is removed
fa = fold_assignment(modelfit(m_cv))

# The entries of `fa` must partition the segments of the new fit (which are detected the same way as for the full model, so the split applies directly).
m_red = fit(UnfoldModel, [Any => (f_red, bf)], evts, eeg_cont; solver = solver_cv(n_folds = 5, fold_assignment = fa))

# Both fits now hold out identical segments per fold, so the fold columns of the two fits can be compared directly (e.g. to compute a per-fold delta-R²).
fold_assignment(modelfit(m_red)) == fa
