# # [Save and load Unfold models](@id unfold_io)

# Unfold.jl allows storing Unfold models in a memory-efficient way using (compressed) .jld2 files.

# ## Simulate EEG data and fit an Unfold model
# ```@raw html
# <details>
# <summary>Click to expand</summary>
# ```

# ### Simulate some example data using UnfoldSim.jl
using UnfoldSim
data, events = UnfoldSim.predef_eeg(; n_repeats=10)
first(events,5)

# ### Fit an Unfold model
using Unfold
basisfunction = firbasis(τ = (-0.5,1.0), sfreq = 100, name = "stimulus")
f = @formula 0 ~ 1 + condition + continuous 
bfDict = Dict(Any => (f,basisfunction))
m = fit(UnfoldModel, bfDict, events, data);

# ```@raw html
# </details >
# ```

# ## Save and load the fitted Unfold model

# The following code saves the model in a compressed .jld2 file. The default option of the `save` function is `compress=false`.
# For memory efficiency the designmatrix is set to missing. If needed, it can be reconstructed when loading the model.
save_path = mktempdir(; cleanup=false) # create a temporary directory for the example
save(joinpath(save_path, "m_compressed.jld2"), m; compress=true);

# The `load` function allows to retrieve the model again.
# By default, the designmatrix is reconstructed. If it is not needed set `generate_Xs=false`` which improves time-efficiency.
m_loaded = load(joinpath(save_path, "m_compressed.jld2"), UnfoldModel, generate_Xs=true);