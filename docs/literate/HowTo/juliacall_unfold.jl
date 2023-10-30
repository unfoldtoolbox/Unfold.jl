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