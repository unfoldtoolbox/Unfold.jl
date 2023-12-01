#  # [Effects](@id effects)
# Effects are super useful to understand the actual modelfits. If you are an EEG-Researcher, you can think of effects as the "modelled ERPs", and the coefficients as the "difference-waves".
# In some way, we are fitting a model with coefficients and then try to get back the "original" ERPs - of course typically with some effect adjusted, overlap removed or similar - else why bother ;)

# Setup some packages

using Unfold
using DataFrames
using Random
using CSV
using UnfoldMakie
using UnfoldSim
using UnfoldMakie

# # Setup things
# Let's generate some data and fit a model of a 2-level categorical and a continuous predictor with interaction.
data,evts = UnfoldSim.predef_eeg(;noiselevel=8)

basisfunction = firbasis(τ = (-0.1, 0.5), sfreq = 100, name = "basisA")


f = @formula 0 ~ 1+condition+continuous # 1

m = fit(UnfoldModel, Dict(Any=>(f,basisfunction)), evts, data,eventcolumn="type")

# Plot the results
plot_erp(coeftable(m))

# As expected, we get four lines - the interaction is flat, the slope of the continuous is around 4, the categorical effect is at 3 and the intercept at 0 (everything is dummy coded by default)
# ### Effects
# A convenience function is `effects`. It allows to specify effects on specific levels, while setting non-specified ones to a typical value (usually the mean)

eff = effects(Dict(:condition => ["car","face"]),m)
plot_erp(eff;mapping=(;color=:condition,))

# We can also generate continuous predictions
eff = effects(Dict(:continuous => -5:0.5:5),m)
plot_erp(eff;mapping=(;color=:continuous,group=:continuous=>nonnumeric),categorical_color=false,categorical_group=false)

# or split it up by condition
eff = effects(Dict(:condition=>["car","face"],:continuous => -5:2:5),m)
plot_erp(eff;mapping=(;color=:condition,col=:continuous=>nonnumeric))

# ## What is typical anyway?
# The user can specify the typical function applied to the covariates/factors that are marginalized over. This offers even greater flexibility.
# Note that this is rarely necessary, in most applications the mean will be assumed.
eff_max = effects(Dict(:condition=>["car","face"]),m;typical=maximum) # typical continuous value fixed to 1
eff_max.typical .= :maximum
eff = effects(Dict(:condition=>["car","face"]),m)
eff.typical .= :mean # mean is the default

plot_erp(vcat(eff,eff_max);mapping=(;color=:condition,col=:typical))
