#  # [Effects](@id effects)
# Effects are super useful to understand the actual modelfits. If you are an EEG researcher, you can think of the coefficients as the "difference waves" and the (marginal) effects as the "modelled ERP evaluated at a certain predictor value combination".
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
# Let's generate some data and fit a model of a 2-level categorical and a continuous predictor without interaction.
data, evts = UnfoldSim.predef_eeg(; noiselevel = 8)

basisfunction = firbasis(τ = (-0.1, 0.5), sfreq = 100, name = "basisA")


f = @formula 0 ~ 1 + condition + continuous # 1

m = fit(UnfoldModel, [Any => (f, basisfunction)], evts, data, eventcolumn = "type")

# Plot the results
plot_erp(coeftable(m))

# As expected, we get three lines representing the coefficients - the slope of the continuous peaks at around 1µV / 1-unit-change, the categorical effect around 3µV and the intercept shows the reference-category with the typical p1/n1/p3 complex.
#
# ### Effects
# In order to better understand the actual predicted ERP curves, often researchers had to do manual contrasts. Remember that a linear model is y = X*b, which allows (after `b` was estimated) to input a so-called `contrast` vector for X. You might know this in the form of `[1,0,-1,1]` or similar form. Quite errorprone for larger models!
# Here the convenience function `effects` comes into play. It allows to specify the contrast-vectors by providing actual levels of the design. If multiple variables are provided, it calculates all possible combinations. If a variable is skipped, 
# it automatically sets it to it's `typical value` (usually the `mean`, but for categorical variables could also be others - the emmeans package has quite some discussion on this).

eff = effects(Dict(:condition => ["car", "face"]), m)
plot_erp(eff; mapping = (; color = :condition,))

# We can also generate continuous predictions
eff = effects(Dict(:continuous => -5:0.5:5), m)
plot_erp(
    eff;
    mapping = (; color = :continuous, group = :continuous => nonnumeric),
    categorical_color = false,
    categorical_group = false,
)

# or split it up by condition and calculate all combinations automagically.

eff = effects(Dict(:condition => ["car", "face"], :continuous => -5:2:5), m)
plot_erp(eff; mapping = (; color = :condition, col = :continuous => nonnumeric))

# ## What is typical anyway?
# The user can specify the `typical function` applied to the covariates/factors that are marginalized over. This offers even greater flexibility in defining what is "typical", rather than only take the mean over a predictor as `typical`
#
# Note that this is rarely necessary, in most applications, readers would assume the mean. But for e.g. skewed distributions, it could be interesting to look at e.g. the `mode`, or with outliers, the `median`, or `mean(winsor)`
#
# As an illustration we will employ the maximum over the `continuous` predictor

eff_max = effects(Dict(:condition => ["car", "face"]), m; typical = maximum)
eff_max.typical .= :maximum
eff = effects(Dict(:condition => ["car", "face"]), m)
eff.typical .= :mean # mean is the default

plot_erp(vcat(eff, eff_max); mapping = (; color = :condition, col = :typical))
