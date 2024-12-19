# # Solver implementation
# This document describes how the `solver_main` is implemented and how to add custom solvers.

# some setup
using Unfold, UnfoldSim, CairoMakie
using LinearAlgebra: cholesky
# ## Solver main
# This function gis a eneral purpose solver-wrapper function. It calls  `prepare_fun` and iterates over the first dimension of `data`, repeatedly calling the `solver_fun`.
#
# Without any bells and whistles (progress, history etc.) the function roughly looks like this:

function _solver_min(X, data; prepare_fun, solver_fun!, stderror = false)
    Ĥ, dataP, prepared = prepare_fun(X, data)
    for ch = 1:size(dataP, 2)
        for t = 1:size(dataP, 3)
            ch == 1 || copyto!(view(Ĥ, ch, :, t), view(Ĥ, ch - 1, :, t))
            solver_fun!(view(Ĥ, ch, :, t), view(dataP, :, ch, t), prepared...)
        end
    end
    modelfit = stderror ? calculate_stderror(X, data, Ĥ) : nothing

    return modelfit

end

# Before diving into the `prepare_fun` and `solver_fun!` functions, let's discuss first the inner loop `t=1:size(dataP,3)`.
# This loop really only comes alife (that is `size(dataP,3)!=1`) if a mass-univariate model is fitted, that is, when ndims(data)==3`. We still have it around for 2D, un-epoched data, to have exactly the same code in both cases.

# ## `prepare_fun``
# This function is the setup / prepare function. It is typically a chain of functions with similar input / output characteristica.
# The first fuction of the chain/pipeline should be a function taking `(X,data)`and returning `(Ĥ::AbstractArray, dataP::AbstractArray, prepared::Tuple)`.
# - `Ĥ` is used to save the beta/parameters inplace
# - `dataP` is the data in format ch x repeat x time (with size(time) = 1 if data initially was a Matrix/2D-array)
# - `prepared` is a tuple of all the other variables needed in the solver-step, e.g. the `pinv(X)` or `X'X` or simply `X`
# The `prepare` function which is typiclly the first, just permutes the data & converts everything to GPU in case `data::CuArray`.
#
# The next function in a pipeline then would take this `(Ĥ::AbstractArray, dataP::AbstractArray, prepared::Tuple)` inputs and process it further.`
#
# ## `solver_fun!`
# This function actually performs the fitting. It takes the inputs `(Ĥ::view(Matrix),data::view(Array),prepared::Tuple)`
# - `Ĥ` is the current beta/parameters view, a vector/slice for one channel and one timepoint
# - `data` is similarly the current data view, a vector/slice for one channel and one timepoint
# - `prepared` is the tuple-output of the `prepare` function.
#
# The `solver_fun!` can output some history of the solver, e.g. a log for iterative solvers.
# ## Example (simple)
# let's setup our own solver:

_my_solver!(Ĥ, data, X) = Ĥ .= Matrix(X) \ data

# let's simulate some data and see this in action
data, evts = UnfoldSim.predef_eeg()
m = fit(
    UnfoldModel,
    @formula(0 ~ 1 + condition),
    evts,
    data,
    firbasis((-0.1, 0.5), 100);
    solver = (x, y) ->
        Unfold.solver_main(x, y; solver_fun! = _my_solver!, show_time = true),
)

# Remember from this table the time for one solve (~700ms on my test-computer) this is the time per channel.
series(coef(m))


# ## Cholesky Example
# !!! note
#     the following function is already implemented in Unfold.jl as well. See `?Unfold.solver_predefined`
# Given that the `prepare` function returns all necessary ingredients, this is a bit simple. So let's make it more complex
#
# for nicety, we need some unpacking wrappers
_prepare_cholesky(all::Tuple) = _prepare_cholesky(all...)
_prepare_cholesky(Ĥ, data, all::Tuple) = _prepare_cholesky(Ĥ, data, all...)
# this function effectively only pre-calculates the cholesky decomposition
_prepare_cholesky(Ĥ, data, Xt, R_xx, R_xy) = (Ĥ, data, (Xt, cholesky(R_xx), R_xy))


# now we have everything to put together our solver-pipeline
_my_prepare =
    (x, y) -> Unfold.prepare(collect(x), y) |> Unfold.prepare_XTX |> _prepare_cholesky

# let's test
# (note we have to reshape the data)
@time _my_prepare(modelmatrix(m), reshape(data, 1, :))

# finally, we need a solver
# this is how we solve the single-channel equation
function _my_cholesky!(beta, data, Xt, XtX_cholesky, R_xy)
    @time Unfold.calc_Rxy!(R_xy, Xt, data)
    @time beta .= XtX_cholesky \ R_xy
end

m = fit(
    UnfoldModel,
    @formula(0 ~ 1 + condition),
    evts,
    data,
    firbasis((-0.1, 0.5), 100);
    solver = (x, y) -> Unfold.solver_main(
        x,
        y;
        prepare_fun = _my_prepare,
        solver_fun! = _my_cholesky!,
        show_time = true,
    ),
)

# This (on my test-computer) took only 97ms per channel, so it is ~7x faster per channel.

series(coef(m))
