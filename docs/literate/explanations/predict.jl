# # The predict-family
## Setup
using Unfold
using UnfoldSim
using CairoMakie

dat, evts = UnfoldSim.predef_eeg(noiselevel = 5)
design = [
    "car" => (@formula(0 ~ 1 + continuous), firbasis(τ = (-0.5, 1), sfreq = 100)),
    "face" => (@formula(0 ~ 1 + continuous), firbasis(τ = (-0.3, 0.5), sfreq = 100)),
]

m = fit(UnfoldModel, design, evts, dat; eventcolumn = :condition);

# # Overview
# In a linear model $EEG = Xβ + e$, predictions boil down to finding $\hat{EEG} = Xβ$, thus EEG data without any error term. 
# Different types of predictions can be generated by modifying the $X$ accordingly.
#
# !!! note
#     We simulated only a single channel, all results generalize to the multi channel case
#
# # Different types of predictions
# ## Time-Continuous case
# Let's start with the cases, where the EEG was not epoched before using Unfold, i.e. the EEG was analysed with e.g. FIR-deconvolution

# ### Continuous EEG
# This returns $EEG = Xβ$ - the continuous modelled EEG
p = predict(m) # same as predict(m, overlap = true)
lines(p[1, 1:1000])


# ### No-overlap
# This results in one prediction Array per event
p = predict(m, overlap = false)
size(p)
# Each has the size:
size(p[1])

# Visualizing the 1000 events
series(range(-0.5, 1, step = 1 / 100), p[1][1, :, :]', solid_color = :orange)
series!(range(-0.3, 0.5, step = 1 / 100), p[2][1, :, :]', solid_color = :teal)
current_figure()

# At ~0.3s we can see a split between the predicted EEG single trials into 10 "strands" - this is the granularity of our continuous predictor. You could use `effects` to improve upon this granularity / customize it.
# ### with-overlap, epoched
# Sometimes helpful is to add in the overlap we removed via the deconvolution.
p = predict(m, epoch_to = ["car"], eventcolumn = :condition)
series(range(-0.5, 1, step = 1 / 100), p[1, :, :]', solid_color = :orange)

# ### partial-overlap
# We can also include/exclude certain events with "partial-overlap"
p_car = predict(m, keep_basis = ["car"], eventcolumn = :condition)
p_face = predict(m, exclude_basis = ["car"], eventcolumn = :condition) # same as keep_basis=["face"]
f = lines(p_car[1, 1:1000])
lines!(p_face[1, 1:1000])
f

# In the plot, we see the two partial predictions for car and face. They are respectively "0" outsie the basisfunction windows

# !!! note
#      The above options can be combined as well, e.g. to get an `epoch_to`, `exclude_basis` version. `epoch_timewindow` can be specified as well.