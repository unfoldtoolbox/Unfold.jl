## Effects
# Effects are super useful to understand the actual modelfits. If you are an EEG-Researcher, you can think of effects as the "modelled ERPs", and the coefficients as the "difference-waves".
# In some way, we are fitting a model with coefficients and then try to get back the "original" ERPs - of course typically with some effect adjusted, overlap removed or similar - else why bother ;)

# Define some packages

using Unfold
using DataFrames
using Random
using CSV


# # Setup things
# Let's load some data and fit a model of a 2-level categorical and a continuous predictor with interaction.
include(joinpath(dirname(pathof(Unfold)), "../test/test_utilities.jl") ) # to load data

data, evts = loadtestdata("test_case_3a")
basisfunction = firbasis(τ = (-0.5, 1.5), sfreq = 20, name = "basisA")

f = @formula 0 ~ 1+conditionA*continuousA # 1

m = fit(UnfoldModel, Dict(Any=>(f,basisfunction)), evts, data,eventcolumn="type")

# Plot the results
plot_results(coeftable(m))

# As expected, we get four lines - the interaction is flat, the slope of the continuous is around 4, the categorical effect is at 3 and the intercept at 0 (everything is dummy coded by default)
#
# A convenience function is [effects](@ref).

eff = Unfold.effects(Dict(:conditionA => [0,1],:continuousA=>[0]),m)