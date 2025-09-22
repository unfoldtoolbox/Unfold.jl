# Calling Unfold.jl from R

Julia code can be called from within R using `JuliaCall`. Here is a **barebone**  tutorial how to do so. If you ever run into issues, please write us an github issue - we have not extensively used this bridge in the past.

## Install JuliaCall

```{R}
install.packages("JuliaCall")
install_julia() # installs Julia
```

Once JuliaCall and Julia are installed, we can start installing Unfold and calling the library.

Two ways to do so: The "proper" way, with a reproducible environment:

```{R}
julia <- julia_setup() # setup julia
path_to_env = '/tmp/my_julia_env' # could be any path to a Julia Project.toml
julia_eval(paste('import Pkg;','Pkg.activate("',path_to_env,'");Pkg.instantiate()'))

# if Unfold is not yet installed
julia_eval('Pkg.add("Unfold");')

```

The fast way:

```{R}
julia_install_package_if_needed("Unfold")
```

We are now ready to create some fake data (which you'd replace with your own of course), and call Unfold.jl. We will only provide an example for an epoched, mass-univariate analysis.

```{R}

R_data <- matrix(rnorm(200), nrow = 10, ncol = 20) # 10 timepoints, 20 trials
R_df <- data.frame(trial=1:20)
julia_assign("jl_data",R_data) # move data to julia
julia_assign("jl_df",R_df) # move DF to julia

julia_library("Unfold") # load Unfold - takes 2-3s

julia_eval("m = fit(UnfoldModel,[Any=>(@formula(0~1+trial),1:10)],jl_df,jl_data)") # first time slow, next time fast!

eff = julia_eval("effects(Dict(:trial=>[1,5,10]),m)") # calculate marginal means and return DataFrame

```

Wonderful - you now have a tidy DataFrame in R, with the marginal effect of your unfold model evaluated at `trial = 1,5,10`.

Note that the [Julia's formula syntax](https://juliastats.org/StatsModels.jl/stable/api/#Formulae-and-terms) is similar, but not identical to R's formulas. ``
