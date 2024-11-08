#  # [Marginal effects](@id effects)
# [Marginal effect plots](https://library.virginia.edu/data/articles/a-beginners-guide-to-marginal-effects) are useful for understanding model fits. 

# If you are an EEG researcher, you can think of the coefficients as the 'difference waves' and the (marginal) effects as the 'modelled ERP evaluated at a certain predictor value combination'.
# In some way, we are fitting a model with coefficients, receiving intercepts and slopes, and then try to recover the 'classical' ERPs in their "data-domain", typically with some effect adjustment, overlap removal, or similar.

# # Setup things
# Setup some packages

using Unfold
using DataFrames
using Random
using CSV
using UnfoldMakie
using UnfoldSim
using UnfoldMakie
using DisplayAs # hide

# Generate data and fit a model with a 2-level categorical predictor and a continuous predictor without interaction.
data, evts = UnfoldSim.predef_eeg(; noiselevel = 8)

basisfunction = firbasis(τ = (-0.1, 0.5), sfreq = 100; interpolate = false)

f = @formula 0 ~ 1 + condition + continuous # 1

m = fit(UnfoldModel, [Any => (f, basisfunction)], evts, data, eventcolumn = "type")
m |> DisplayAs.withcontext(:is_pluto => true) # hide

# Plot the results
plot_erp(coeftable(m))

#= 
The coefficients are represented by three lines on a figure:
- the intercept showing the reference category for a typical p1/n1/p3 ERP components;
- the slope of continuous variables with 1µV range;
- the effect of categorical variabe with 3µV range. 
=#

# ### Effects function
# In order to better understand the actual predicted ERP curves, often researchers had to do manual contrasts. Remember that a linear model is `y = X * b`, which allows (after `b` was estimated) to input a so-called `contrast` vector for X. You might know this in the form of `[1, 0, -1, 1]` or similar form. However, for larger models, this method can be prone to errors.

# The `effects` function is a convenient way to specify contrast vectors by providing the actual levels of the experimental design. It can be used to calculate all possible combinations of multiple variables. 

# If a predictor-variable is not specified here, the function will automatically set it to its typical value. This value is usually the `mean`, but for categorical variables, it could be something else. The R package `emmeans` has a lot of discussion on this topic.

eff = effects(Dict(:condition => ["car", "face"]), m)
plot_erp(eff; mapping = (; color = :condition,))

# We can also generate continuous predictions:
eff = effects(Dict(:continuous => -5:0.5:5), m)
plot_erp(
    eff;
    mapping = (; color = :continuous, group = :continuous => nonnumeric),
    categorical_color = false,
    categorical_group = false,
)

# Or we can split our marginal effects by condition and calculate all combinations "automagically".

eff = effects(Dict(:condition => ["car", "face"], :continuous => -5:2:5), m)
plot_erp(eff; mapping = (; color = :condition, col = :continuous))

# ## What is typical anyway?
# The `effects` function includes an argument called `typical`, which specifies the function applied to the marginalized covariates/factors. The default value is `mean`, which is usually sufficient for analysis.
#
# However, for skewed distributions, it may be more appropriate to use the `mode`, while for outliers, the `median` or `winsor` mean may be more appropriate.
#
# To illustrate, we will use the `maximum` function on the `continuous` predictor.

eff_max = effects(Dict(:condition => ["car", "face"]), m; typical = maximum)
eff_max.typical .= :maximum
eff = effects(Dict(:condition => ["car", "face"]), m)
eff.typical .= :mean # mean is the default

plot_erp(vcat(eff, eff_max); mapping = (; color = :condition, col = :typical))
