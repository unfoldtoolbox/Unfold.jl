#  # [FIR-Scaled duration predictors](@id hassall-duration)
using Unfold
using Interpolations

using UnfoldSim
using UnfoldMakie, CairoMakie
using DataFrames
using DisplayAs # hide

data, evts = UnfoldSim.predef_eeg(sfreq = 10, n_repeats = 1)

evts.duration = 5:24


# putting `scale_duration = Interpolation.Linear()` will introduce a Cameron-Hassall 2022 PNAS- Style basisfunction, that scales with the `:duration` column
basisfunction = firbasis(τ = (-1, 2), sfreq = 5, scale_duration = Interpolations.Linear())

# Two examples with `duration = 10`
Unfold.kernel(basisfunction, [0, 10])
# and `duration = 20`
Unfold.kernel(basisfunction, [0, 20])

# let's fit a model
f = @formula 0 ~ 1 + condition
bf_vec = [Any => (f, basisfunction)]
m = fit(UnfoldModel, bf_vec, evts, data; eventfields = [:latency, :duration]);


## currently bugged for small matrices
## plot_designmatrix(designmatrix(m))
## thus using
heatmap(Matrix(modelmatrix(m))')
# As one can see, the designmatrix is nicely scaled

# We can predict overlap-corrected results
p = predict(m; overlap = false)[1]
heatmap(p[1, :, :])
# note the `missings` which are displayed as white pixels.

# ## Block-design predictors
# In contrast, it is also possible to put `scale_duration = true` - which wil not scale the matrix as before, but introduce a step-function.

# putting `scale_duration = Interpolation.Linear()` will introduce a Cameron-Hassall 2022 PNAS- Style basisfunction, that scales with the `:duration` column
basisfunction = firbasis(τ = (-1, 2), sfreq = 5, scale_duration = true)
# Two examples with `duration = 10`
Unfold.kernel(basisfunction, [0, 10])
# and `duration = 20`
Unfold.kernel(basisfunction, [0, 20])

# let's fit a model
f = @formula 0 ~ 1 + condition
bf_vec = [Any => (f, basisfunction)]
m = fit(UnfoldModel, bf_vec, evts, data; eventfields = [:latency, :duration]);


heatmap(Matrix(modelmatrix(m))')
# as one can see, now the designmatrix is not stretched - but rather "block"-ed
p = predict(m; overlap = false)[1]
heatmap(p[1, :, :])
