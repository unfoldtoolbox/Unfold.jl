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
# **`@formula`**: we havent found a way to call macros yet, even though we think it should be possible. For now please use `f = jl.seval("@formula(0~1+my+cool+design)")`. Later versions might support something like `f = @formula("0~1+my+cool+design)"` directly
#
# **Specifying the design**: Since Unfold 0.7 we officially switched to the 
# ```julia
# ["eventtypeA"=>(formula,basisfunction),
# "eventtypeB"=>(otherformula,otherbasisfunction)]
# ```
# Array-based syntax, from a Dict-based syntax. Unfortunately, `=>` (a pair) is not supported in Python and one needs to do some rewriting:
# ```python
# jl.convert(jl.Pair,(formula,basisfunction))
# ```
# which makes the code less readable. We are thinking of ways to remedy this - but right now there is now way around. For now, it is also possible to use the old syntax e.g. in python 
# ```python
# {"eventtypeA"=>(formula,basisfunction),"eventtypeB"=>(otherformula,otherbasisfunction)}
# ```
# which is clearly easier to read :)
#
# **`UnfoldSim.design`**: we need a `Dict` with a `Symbol` , one has to do something like `condition_dict_jl = {convert(jl.Symbol,"condA"):["car", "face"]}` to do so. We will [try to allow strings}(https://github.com/unfoldtoolbox/UnfoldSim.jl/issues/96) here as well, removing this constraint.
#
# When preprocessing your raw data through MNE Python, take the following into consideration:
# The [Raw object](https://mne.tools/stable/generated/mne.io.Raw.html) contains the [first_samp](https://mne.tools/stable/documentation/glossary.html#term-first_samp) attribute which is an integer representing the number of time samples that passed between the onset of the hardware acquisition system and the time when data recording started.
# The Raw data doesn't include these time samples, meaning that the first sample is the beginning of the data aquisition.
# From the Raw object you can obtain an events array from the annotations through [mne.events_from_annotations()](https://mne.tools/stable/generated/mne.events_from_annotations.html).
# The events array, however, does include first_samp, meaning that the annotated events in events array don't match the Raw object anymore.
# Alternatively, it might be easier to convert the annotations to a pandas dataframe directly (`to_data_frame()`), or even better, load the "*_events.tsv" from a BIDS dataset. In the latter case, all columns will be preserved, which MNE's read_annotation drops.
