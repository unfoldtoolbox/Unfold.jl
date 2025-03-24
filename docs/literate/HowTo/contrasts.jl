#  # [Contrast Coding](@id contrasts)
using CairoMakie
using Unfold
using UnfoldMakie
using UnfoldSim


## Contrast coding
# Unfold.jl uses the `StatsModels` package for the formula interface. This allows for a wide range of contrast coding schemes. For a full tutorial, please see [the StatsModels docs](https://juliastats.org/StatsModels.jl/stable/contrasts/).
# Please read their tutorial, as a motivation of why one would change the contrast coding scheme is outside of the realms of this package and more a basic linear regression question.


# !!! hint
#      Given we have a nice `effects` implementation (mimicking `emmeans` and similar packages), coding schema is typically less important.

# Here we will show a simple example of how to change the contrast coding scheme. We will use the `condition` variable, which has two levels, `A` and `B`. We will change the contrast coding from `Dummy` aka `Reference` aka 0/1 coding to `Sum` coding, which is the default in R.
eeg, evts = UnfoldSim.predef_eeg(noiselevel = 0)
f = @formula 0 ~ 1 + condition
basis = firbasis((-0.1, 0.6), 100)
m_dummy = fit(UnfoldModel, f, evts, eeg, basis)
m_effec =
    fit(UnfoldModel, f, evts, eeg, basis; contrasts = Dict(:condition => EffectsCoding()))


# we could directly inspect the designmatrix
modelmatrix(m_dummy, false)[1][1:5, :]

# and the effects coding
modelmatrix(m_effec, false)[1][1:5, :]

# To confirm the difference in the actual fit, let's visualize them
c_d = coeftable(m_dummy)
c_e = coeftable(m_effec)
c_d.group .= "Dummy Coding"
c_e.group .= "Effects Coding"
c = vcat(c_d, c_e)

plot_erp(c; mapping = (; color = :coefname, col = :group))

# As expected, the effects-coding slope of `condition: face` is half the size of the dummy-coding one (because -1/1 coding was used).
