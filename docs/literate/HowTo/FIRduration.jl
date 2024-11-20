using Unfold

using UnfoldSim
using UnfoldMakie, CairoMakie
using DataFrames
using DisplayAs # hide

data, evts = UnfoldSim.predef_eeg(sfreq = 5, n_repeats = 1)

evts.duration = 5:24 #rand(2:8, size(evts, 1))

# putting `scale_duration = true` will introduce a Cameron-Hassel Style basisfunction, that scales with the `:duration` column
basisfunction = firbasis(τ = (-1, 2), sfreq = 5, scale_duration = true)

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