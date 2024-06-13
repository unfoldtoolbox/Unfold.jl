# # Using Unfold.jl from Python
# it is straight forward to call Unfold from Python using `JuliaCall`.

# ## Quick start

# Create a Python environment and install JuliaCall.
# ```bash
# pip install juliacall
# ```

# Create a Julia environment and install Unfold
# ```python
# # Import the Julia package manager
# from juliacall import Pkg as jlPkg

# # Activate the environment in the current folder
# jlPkg.activate(".")

# # Install Unfold (in the activated environment)
# jlPkg.add("Unfold")
# ```

# Import Julia's main module and Unfold
# ```python
# # Import Julia's Main module
# from juliacall import Main as jl

# # Import Unfold
# # The function seval() can be used to evaluate a piece of Julia code given as a string
# jl.seval("using Unfold")
# Unfold = jl.Unfold # simplify name
# ```

# Now you can use all Unfold functions as for example
# ```python
# dummy_model = Unfold.UnfoldLinearModel(jl.Dict())
# ```

# ## Example: Unfold model fitting from Python
# In this [notebook](https://github.com/unfoldtoolbox/Unfold.jl/blob/main/docs/src/HowTo/juliacall_unfold.ipynb), you can find a more detailed example of how to use Unfold from Python to load data, fit an Unfold model and visualise the results in Python.

# ## Important limitations
# Python doesnt not offer the full expressions that are available in Julia. So there are some things you need to give special attention:
#
# `@formula` - we havent found a way to call macros yet, even though we think it should be possible. For now please use `f = jl.seval("@formula(0~1+my+cool+design)")`. Later versions might support something like f = `@formula("0~1+my+cool+design)"` directly
#
# specifying the design: Since Unfold 0.7 we officially switched to the `["eventtypeA"=>(formula,basisfunction),"eventtypeB"=>(otherformula,otherbasisfunction)]` syntax from a Dict-based syntax. Unfortunately, `=>` (a pair) is not supported in Python and one needs to do some rewriting:
`jl.convert(jl.Pair,(formula,basisfunction))` which makes the code less readable. We are thinking of ways to remedy this - but right now there is now way around. For now, it is also possible to use the old syntax {"eventtypeA"=>(formula,basisfunction),"eventtypeB"=>(otherformula,otherbasisfunction)} which is clearly easier to read :)
#
# `UnfoldSim.design` special case in UnfoldSim design, we need a `Dict` with a `Symbol` , one has to do something like `condition_dict_jl = {convert(jl.Symbol,"condA"):["car", "face"]}` to do so. We will try to allow strings here as well, removing this constraint.
