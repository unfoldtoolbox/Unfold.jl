# # Removing BAD data
#
# Sometimes one has bad data (artefacts etc.) in the data that one wants to remove prior to fitting.
# ### Setup
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```
## Load required packages
using UnfoldSim
using Unfold
using CairoMakie, UnfoldMakie
using Random
# ```@raw html
# </details >
# ```
# ## Simulation
# Let's start with a simulation of some data, and add some bad-data.
rng = MersenneTwister(1)
data, evts = UnfoldSim.predef_eeg(n_repeats = 1)
ix = 100:500
data[ix] .+= 50 .* (rand(rng, length(ix)) .- 0.5)

f = lines(data[1:1000])
f

# Clearly this data is bad and we should rather remove it!
# We can use Julias `missing` data-type to indicate those portions.
using Missings
data_missing = allowmissing(data)
data_missing[ix] .= missing

lines(f.figure[2, 1], data_missing[1:1000])
f
# Nice, the data is sucessfully marked as missing, without changing the overall size / timing of the data
#
# We can now fit a model with, and without the noise data. Unfold automatically detects missing data, and does not use it in the modelfit.
m = fit(UnfoldModel, @formula(0 ~ 1), evts, data, firbasis((-0.3, 0.5), 100))
m_missing =
    fit(UnfoldModel, @formula(0 ~ 1), evts, data_missing, firbasis((-0.3, 0.5), 100))
c = coeftable(m)
c.group .= "with bad data"
c_missing = coeftable(m_missing)
c_missing.group .= "bad data removed"

plot_erp(vcat(c, c_missing), mapping = (; color = :group))


# Currently no helper function exists to translate e.g. MNE BAD segments/annotations automatically, but pull requests are very welcome!
