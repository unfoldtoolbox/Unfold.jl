"""
    $(SIGNATURES)
default solvers.
- If data is continuous (2D), we solve Xb = y via lsmr
- If data is epoched (3D) we solve Xb = y via pinv

We highly recommend to check out `solver_predefined` for faster options by rather solving X'Xb = X'y via QR, cholesky, pinv or `\`-solver. A benchmark is available in the online documentation.

Please see `?solver_main` for keyword arguments of the solver (like `stderror`, `multithreading`, `show_time`, `show_progress`)

"""
solver_default(X, y::AbstractMatrix; kwargs...) = solver_main(X, y; kwargs...)
function solver_default(X, y::AbstractArray{T,3}; kwargs...) where {T}
    return solver_main(
        X,
        y;
        prepare_fun = (x, y) -> prepare(x, y) |> prepare_pinv,
        solver_fun! = solver_pinv!,
        kwargs...,
    )
end
#---

"""
    $(SIGNATURES)

helper function that returns solver with appropriate prepare-pipelines and fitting solver-functions. X is a (typically sparse) designmatrix, y is a 2D or 3D array.

`solver` : one of `:cg`, `:pinv`, `:intern`, `:qr`, `:cholesky`, `:lsmr` (default)

Only `lsmr` solves Xb = y via an iterative solver and should be more accurate in principle.

The other predefined-solvers solve X'Xb = X'y which is often computationally much cheaper, and because X'X can be precalculated, it should be cheaper to apply.

Testing this empirically is somewhat complicated, as depending on your sparsity structure (≈ your design) and the size of your data (sfreq & minutes) the best solver and the reached accuracy can change quite a bit.

##  GPU
All solvers except :lsmr support GPU calculations. For lsmr on the GPU try `solver_krylov` instead

"""
function solver_predefined(X, y_in::AbstractMatrix; solver = :lsmr, kwargs...)
    prepare_fun = if solver ∈ [:cg, :intern]
        (x, y) -> prepare(x, y) |> prepare_XTX
    elseif solver == :pinv
        (x, y) -> prepare(x, y) |> prepare_XTX |> prepare_pinv
    elseif solver == :qr
        (x, y) -> prepare(x, y) |> prepare_XTX |> prepare_qr
    elseif solver == :cholesky
        (x, y) -> prepare(x, y) |> prepare_XTX |> prepare_cholesky
    elseif solver == :lsmr
        (x, y) -> prepare(x, y)
    else
        error("unknown predefined solver")
    end

    solver_fun! = if solver == :cg
        solver_cg!
    elseif solver == :pinv
        solver_pinv!
    elseif solver == :intern
        solver_intern!
    elseif solver == :qr
        solver_qr!
    elseif solver == :cholesky
        solver_cholesky!
    elseif solver == :lsmr
        solver_lsmr!
    else
        error("unknown predefined solver")
    end

    solver_main(X, y_in; prepare_fun, solver_fun!, kwargs...)
end
"""
    $(SIGNATURES)
general purpose solver function. Calls  `prepare_fun` and iterates over the first dimension of `data`, repeatedly calling the `solver_fun`.

- Output of the `solver_fun!` is saved as solver_history
- Solutions of channels are initialized from previous channels (starting from ch>=2)

## Keyword Arguments
- `prepare_fun` : `(X,data)->prepare(X,data)` by default. Can be easily used to chain multiple preparation steps and thereby overload the solver.
Example: `prepare_fun = (X,data)-> prepare(X,data) |> prepare_XTX |> prepare_qr`. 
Each `prepare` function needs to `return Ĥ, data, (Xt, R_xx)` where `Ĥ` is used to save the solver results, data is a permuted version of the data (following X x ch x time, with X being continuous time, and length(time) = 1 for 2D arrays, and X being trials for 3D arrays and length(time) = length of epochs)

- `solver_fun!` : a function taking the inputs `(Ĥ::vectorview(Matrix),data::vectorview(Array),prepared::Tuple)`, the output of the `prepare_fun`, but along a single channel & `time` (in case of 3D data)
- `multithreading` : `true` by default. Multithread over channels? Not recommended for GPU calculations
- `show_time` : `false` by default. Timer-Outputs for prepare and channel stage. Helpful to select best performing solver
- `show_progress` : `true` by default. Show progressbar over channels / timepoints (for 3D)
- `stderror` : `false` by default. Calculate the stderror. Warning: for 2D arrays, the stderror is very likely too small as autocorrelation is not taking into account! 

"""

function solver_main(
    X,
    data::AbstractArray{T,N};
    prepare_fun = prepare,
    solver_fun! = solver_lsmr!,
    show_time = false,
    show_progress = true,
    multithreading = false,
    stderror = false,
) where {T,N}
    to = TimerOutput()
    solver_history = Array{Any}(undef, size(data, 1))



    @timeit to "prepare" Ĥ, dataP, prepared = prepare_fun(X, data)
    # dataP are permuted data, with repetition x channel [x time]


    p = Progress(size(dataP, 2) * size(dataP, 3); enabled = show_progress)

    @debug "solver types" typeof(Ĥ) typeof(dataP) typeof(prepared[1]) typeof(data) typeof(X)
    # we have to put the solve-timer outside, because it clashes with multithreading
    @timeit to "solve" Unfold.@maybe_threads multithreading for ch = 1:size(dataP, 2)
        # init solver with previous fit
        ch == 1 || multithreading || copyto!(view(Ĥ, ch, :, :), view(Ĥ, ch - 1, :, :))

        for t = 1:size(dataP, 3) # for 2D, this will aways be one; for 3D this will traverse timepoints
            view_Ĥ = view(Ĥ, ch, :, t)
            view_dataP = view(dataP, :, ch, t)
            # run the solver function
            h = solver_fun!(view_Ĥ, view_dataP, prepared...)

            # save the history
            solver_history[ch] = h

            next!(p)
        end
    end
    finish!(p)

    fit_info = [prepare_fun, solver_fun!, solver_history]
    @timeit to "stderror&modelfit" if stderror
        _Ĥ = _permuteback(Ĥ)
        _stderror = calculate_stderror(X, data, _Ĥ)
        modelfit = Unfold.LinearModelFit{T,N}(Array{T}(_Ĥ), fit_info, Array{T}(_stderror))
    else
        modelfit = Unfold.LinearModelFit{T,N}(Array{T}(_permuteback(Ĥ)), fit_info)
    end

    show_time ? begin
        println(to)
        #        println("")
    end : nothing
    return modelfit

end
_permuteback(Ĥ::AbstractArray) = permutedims(Ĥ, (1, 3, 2))
_permuteback(Ĥ::AbstractMatrix) = Ĥ


#----
using StatsBase: var
function calculate_stderror(
    Xdc,
    data::AbstractMatrix,
    beta::AbstractArray{T},
) where {T<:Union{Missing,<:Number}}

    # remove missings
    ix = any(.!ismissing.(data), dims = 1)[1, :]
    if length(ix) != size(data, 2)
        @warn(
            "Limitation: Missing data are calculated over all channels for standard error"
        )
    end

    data = data[:, ix]
    Xdc = Xdc[ix, :]

    # Hat matrix only once
    hat_prime = inv(disallowmissing(Matrix(Xdc' * Xdc)))
    # Calculate residual variance
    @warn(
        "Autocorrelation was NOT taken into account. Therefore SE are UNRELIABLE. Use at your own discretion"
    )

    se = Array{T}(undef, size(data, 1), size(Xdc, 2))
    for ch = 1:size(data, 1)
        residualVar = var(data[ch, :] .- Xdc * beta[ch, :])
        @assert(!isnan(residualVar), "residual Variance was NaN")
        hat = hat_prime .* residualVar
        #se = sqrt(diag(cfg.contrast(:,:)*hat*cfg.contrast(:,:)'));
        se[ch, :] = sqrt.(diag(hat))
    end
    return se
end
function calculate_stderror(
    X,
    data::AbstractArray{T1,3},
    beta::AbstractArray{T2},
) where {T1<:Union{Missing,<:Number},T2<:Union{Missing,<:Number}}
    #function calculate_stderror(Xdc,data::AbstractArray{T,2},beta) where {T<:Union{Missing, <:Number}}  
    X = disallowmissing(X)
    # Hat matrix
    hat_prime = inv(Matrix(X' * X))

    se = Array{T2}(undef, size(data, 1), size(data, 2), size(X, 2))
    for ch = 1:size(data, 1)
        for t = 1:size(data, 2)
            ix = .!ismissing.(data[ch, t, :])
            # Calculate residual variance
            #@show size(data) size(X) size(beta)
            residualVar = var(data[ch, t, ix] .- X[ix, :] * beta[ch, t, :])
            @assert(!isnan(residualVar), "residual Variance was NaN")
            hat = hat_prime .* residualVar
            #se = sqrt(diag(cfg.contrast(:,:)*hat*cfg.contrast(:,:)'));
            se[ch, t, :] = sqrt.(diag(hat))
        end
    end
    return se
end


