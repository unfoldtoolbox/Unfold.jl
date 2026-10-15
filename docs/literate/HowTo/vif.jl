# # [Variance Inflation Factors (VIF)](@id vif)
# @warn "This page was drafted with the help of an LLM. Please [report any errors or inconsistencies](https://github.com/unfoldtoolbox/Unfold.jl/issues)."
#
# The **Variance Inflation Factor (VIF)** is a diagnostic for *collinearity* in a linear model's designmatrix. For each coefficient, it measures how much the variance of its estimate is inflated compared to the case where that predictor would be orthogonal to all other predictors:
#
# $$\text{VIF}_j = \frac{1}{1 - R^2_j}$$
#
# where $R^2_j$ is the coefficient of determination from regressing predictor $j$ on all other predictors. `VIF = 1` means no collinearity; the larger the VIF, the more inflated the standard error of the corresponding coefficient. A common rule of thumb is to flag predictors with `VIF > 10`.
#
# Unfold.jl extends [`vif`](@ref) for `UnfoldModel`s, and it works for both *time-expanded* and *mass-univariate* models. Since it only uses the designmatrix (`modelmatrix(m)`), it also works for models that were only built, never fitted (`fit(...; fit = false)`).
#
# ```@raw html
# <details>
# <summary>Set things up</summary>
# ```

using Unfold
using UnfoldSim
using DataFrames
using Random

# ```@raw html
# </details >
# ```

# We simulate some example data using `UnfoldSim.jl` and extract epochs, i.e. one data point per event and time point:
data, evts = UnfoldSim.predef_eeg()
data_e, times = Unfold.epoch(data = data, tbl = evts, τ = (-1, 1.9), sfreq = 100)
f = @formula 0 ~ 1 + condition + continuous

# # Mass-univariate models
# `vif(m)` returns a `DataFrame` with one row per column of the designmatrix (i.e. per term). Note that we don't even need to fit the model: the VIF only requires the designmatrix.
m_e = fit(UnfoldModel, [Any => (f, times)], evts, data_e; fit = false)
vif(m_e)

# The `(Intercept)` column is constant, so it has no defined VIF and is reported as `Inf`. The other two predictors are uncorrelated, so their VIFs are (nearly) 1.

# # Detecting collinearity
# Let's create a collinearity problem: we duplicate the `continuous` predictor and add a bit of independent noise to both copies:
rng = MersenneTwister(1)
evts2 = deepcopy(evts)
true_c = evts2.continuous
evts2.continuous .= true_c .+ 0.01 .* randn(rng, nrow(evts2))
evts2.continuous2 = true_c .+ 0.01 .* randn(rng, nrow(evts2))

# The VIFs of the two near-duplicate terms now explode (~50,000), far beyond the rule of thumb of 10, while the unrelated `condition` stays at 1 — the diagnostic both flags the problem and tells us which coefficients are unstable. In time-expanded models, the same problem would show up on every tap of the affected predictors.
f2 = @formula 0 ~ 1 + condition + continuous + continuous2
m2 = fit(UnfoldModel, [Any => (f2, times)], evts2, data_e; fit = false)
vif(m2)

# # Multiple events
# For models with multiple events, the output gains an `:eventname` column:
evts3 = deepcopy(evts)
evts3.type .= repeat(["A", "B"], nrow(evts3) ÷ 2)
m_ab = fit(
    UnfoldModel,
    ["A" => (f, times), "B" => (f, times)],
    evts3,
    data_e;
    eventcolumn = "type",
    fit = false,
)
vif(m_ab)

# # Time-expanded models
# The same function works for time-expanded (FIR) models:
basisfunction = firbasis(τ = (-1, 1), sfreq = 100)
m = fit(UnfoldModel, [Any => (f, basisfunction)], evts, data)
v = vif(m)
first(v, 6)

# For FIR basisfunctions, neighboring taps are strongly correlated, so moderate VIFs (up to ~3.4 here) are to be expected and not a problem. For higher-order splines or strongly overlapping basisfunctions, VIFs can be much larger, but are usually not an issue as long as the solver converges.
#
# Using the rule of thumb, we can filter for the "worst offenders":
v[v.VIF .> 10, :]

# Here the result is an empty `0×2 DataFrame`, i.e. there is no collinearity problem.
