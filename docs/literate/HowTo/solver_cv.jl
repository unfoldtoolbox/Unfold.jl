#  # [Cross-validated Unfold models](@id solver_cv)
# This tutorial shows how to run an Unfold model with k-fold cross-validation.
# !!! important
#      cross validation is not yet implemented for deconvolution models - not because it is hard, but because I didnt do it yet
# # Setup

using Unfold
using UnfoldSim
using Random
using Statistics
using UnfoldMakie, CairoMakie
eeg, evts = UnfoldSim.predef_eeg(; return_epoched = true, noiselevel = 15)

# Define formula and basis function for a mass-univariate model.
f = @formula 0 ~ 1 + condition

# # Cross-validation solver
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
