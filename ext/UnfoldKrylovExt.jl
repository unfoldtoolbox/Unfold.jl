module UnfoldKrylovExt

import Krylov
using CUDA, CUDA.CUSPARSE, SparseArrays
using Unfold
using Missings
using ProgressMeter

"""
Alternative implementation of LSMR using Krylov.jl

The fastest solver (x30 in some cases) with GPU = true

To activate you have to do:
> using Krylov,CUDA 

Difference to solver_default: No suport for per-channel missings. If one sample is missing in any channel, whole channel is removed due a lack of support for Missings in Krylov.

"""
function solver_krylov(
    X,
    data::AbstractArray{T,2};
    GPU = false,
    history = true,
    multithreading = GPU ? false : true,
    show_progress = true,
    stderror = false,
) where {T<:Union{Missing,<:Number}}

    @assert !(multithreading && GPU) "currently no support for both GPU and multi-threading"
    minfo = Array{Any,1}(undef, size(data, 1))

    ix = any(@. !ismissing(data); dims = 1)[1, :]
    X_loop = disallowmissing(X[ix, :])
    data = disallowmissing(view(data, :, ix))


    if GPU
        X_loop = CuSparseMatrixCSC(X_loop)
        lsmr_solver = Krylov.LsmrSolver(size(X_loop)..., CuVector{Float64})
        data = CuArray(data)
    else
        lsmr_solver = Krylov.LsmrSolver(size(X_loop)..., Vector{Float64})
    end

    p = Progress(size(data, 1); enabled = show_progress)

    beta = zeros(T, size(data, 1), size(X, 2)) # had issues with undef
    Unfold.@maybe_threads multithreading for ch = 1:size(data, 1)
        @debug ch
        data_loop = view(data, ch, :)
        # use the previous channel as a starting point
        #ch == 1 || Krylov.warm_start!(lsmr_solver,beta[ch-1,:])#copyto!(view(beta, ch, :), view(beta, ch-1, :))

        Krylov.solve!(lsmr_solver, X_loop, data_loop; history = history)
        beta[ch, :] = Vector(lsmr_solver.x)

        #beta[ch,:],h = Krylov.lsmr(X_loop,data_loop,history=history)


        minfo[ch] = deepcopy(lsmr_solver.stats)
        next!(p)
    end
    finish!(p)
    if stderror
        stderror = Unfold.calculate_stderror(X, data, beta)
        modelfit = Unfold.LinearModelFit{T,2}(beta, ["krylov_lsmr", minfo], stderror)
    else
        modelfit = Unfold.LinearModelFit{T,2}(beta, ["krylov_lsmr", minfo])
    end
    return modelfit
end


end
