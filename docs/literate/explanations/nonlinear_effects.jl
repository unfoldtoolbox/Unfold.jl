# # [Non-linear effects]](@id nonlinear)



using BSplineKit,Unfold
using CairoMakie
using DataFrames
using Random
using Colors



# ## Generating a non-linear signal
# We start with generating 
rng = MersenneTwister(2) # make repeatable
n = 20 # datapoints
evts = DataFrame(:x=>rand(rng,n))
signal = -(3*(evts.x .-0.5)).^2 .+ 0.5 .* rand(rng,n)

plot(evts.x,signal)
#
# Looks perfectly non-linear - great! 
#
# # Compare linear & non-linear fit
# First we have to reshape it to fit to Unfold format, a 3d array,  1 channel x 1 timepoint x 20 datapoints
signal = reshape(signal,length(signal),1,1)
signal = permutedims(signal,[3,2,1])
size(signal)

# Next we define three different models. **linear** **4 splines** and **10 splines**
# note the different formulas, one `x` the other `spl(x,4)`
design_linear = Dict(Any=>	(@formula(0~1+x),[0]))
design_spl3   = Dict(Any=>	(@formula(0~1+spl(x,4)),[0]))
design_spl10  = Dict(Any=>	(@formula(0~1+spl(x,10)),[0])) #hide

# and we fit the parameters
uf_linear = fit(UnfoldModel,design_linear,evts,signal);
uf_spl3 = fit(UnfoldModel,design_spl3,evts,signal);
uf_spl10 = fit(UnfoldModel,design_spl10,evts,signal); #hide

# extract the fitted values using Unfold.effects

p_linear = Unfold.effects(Dict(:x => range(0,stop=1,length=100)),uf_linear);
p_spl3 = Unfold.effects(Dict(:x => range(0,stop=1,length=100)),uf_spl3);
p_spl10 = Unfold.effects(Dict(:x => range(0,stop=1,length=100)),uf_spl10);
first(p_linear,5)

# And plot them
pl = plot(evts.x,signal[1,1,:])
	lines!(p_linear.x,p_linear.yhat)
	lines!(p_spl3.x,p_spl3.yhat)
	lines!(p_spl10.x,p_spl10.yhat)
pl

# What we can clearly see here, that the linear effect (blue) underfits the data, the `spl(x,10)` overfits it, but the `spl(x,4)` fits it perfectly.


# ## Looking under the hood
# Let's have a brief look how the splines manage what they are managing. 
#
# The most important bit to understand is, that we are replacing `x` by a set of coefficients `spl(x)`.
# These new coefficients each tile the range of `x` (in our case from [0-1]) in overlapping areas, each which will be fit by one coefficient.
# Because the ranges are overlapping, we get a smooth function.
#
# Maybe this becomes clear after looking at such a basisfunction:

term_spl = Unfold.formula(uf_spl10).rhs.terms[2]

# This is the spline term.
typeof(term_spl)

# a special type generated in Unfold.jl
#
# it has a field .fun which is the spline function. We can evaluate it at a point
const splFunction =  Base.get_extension(Unfold,:UnfoldBSplineKitExt).splFunction
splFunction([0.2],term_spl)

# each column of this 1 row matrix is a coefficient for our regression model.
lines(disallowmissing(splFunction([0.2],term_spl))[1,:])

# Note: We have to use disallowmissing, because our splines return a `missing` whenever we ask it to return a value outside it's defined range, e.g.:
splFunction([-0.2],term_spl)

# because it never has seen any data outside and can't extrapolate!

# Okay back to our main issue: Let's plot the whole basis set
basisSet = splFunction(0.:0.01:1,term_spl)
basisSet = disallowmissing(basisSet[.!any(ismissing.(basisSet),dims=2)[:,1],:]) # remove missings
ax = Axis(Figure()[1,1])
[lines!(ax,basisSet[:,k]) for k in 1:size(basisSet,2)]
current_figure()

# notice how we flipped the plot around, i.e. now on the x-axis we do not plot the coefficients, but the `x`-values
# Now each line is one basis-function of the spline. 
#
# Unfold returns us one coefficient per basis-function
β=coef(uf_spl10)[1,1,:]
β = Float64.(disallowmissing(β))

# but because we used an intercept, we have to do some remodelling in the basisSet
X = hcat(ones(size(basisSet,1)), basisSet[:,1:5], basisSet[:,7:end])

# now we can weight the spline by the basisfunction
weighted = (β .*X')

# Plotting them makes for a nice looking plot!
ax = Axis(Figure()[1,1])
[lines!(weighted[k,:]) for k = 1:10]
current_figure()


# and sum them up 
lines(sum(weighted,dims=1)[1,:])
plot!(X*β,color="gray") #(same as matrixproduct X*β directly!)
current_figure()

# And this is how you can think about splines
