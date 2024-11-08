
function _lsmr!(beta, X::SparseMatrixCSC, data::AbstractArray{<:Number}, ch)
    _, h = lsmr!(@view(beta[ch, :]), X, @view(data[ch, :]); log = true)
    return h
end
function _lsmr!(
    beta,
    X::SparseMatrixCSC,
    data::AbstractArray{<:Union{<:Number,Missing}},
    ch,
)
    dd = view(data, ch, :)
    ix = @. !ismissing(dd)
    _, h = lsmr!(@view(beta[ch, :]), (X[ix, :]), @view(data[ch, ix]); log = true)
    return h
end

function solver_default(
    X,
    data::AbstractArray{T,2};
    stderror = false,
    multithreading = true,
    show_progress = true,
) where {T<:Union{Missing,<:Number}}
    minfo = Array{IterativeSolvers.ConvergenceHistory,1}(undef, size(data, 1))

    beta = zeros(T, size(data, 1), size(X, 2)) # had issues with undef

    p = Progress(size(data, 1); enabled = show_progress)
    X = SparseMatrixCSC(X) # X s often a SubArray, lsmr really doesnt like indexing into subarrays, one copy needed.
    @maybe_threads multithreading for ch = 1:size(data, 1)

        # use the previous channel as a starting point
        ch == 1 || copyto!(view(beta, ch, :), view(beta, ch - 1, :))

        h = _lsmr!(beta, X, data, ch)
        minfo[ch] = h
        next!(p)
    end
    finish!(p)

    if stderror
        stderror = calculate_stderror(X, data, beta)
        modelfit = Unfold.LinearModelFit(beta, ["lsmr", minfo], stderror)
    else
        modelfit = Unfold.LinearModelFit(beta, ["lsmr", minfo])
    end
    return modelfit
end

function solver_default(
    X,
    data::AbstractArray{T,3};
    stderror = true,
    multithreading = true,
    show_progress = true,
) where {T<:Union{Missing,<:Number}}
    #beta = Array{Union{Missing,Number}}(undef, size(data, 1), size(data, 2), size(X, 2))
    beta = zeros(T, size(data, 1), size(data, 2), size(X, 2))
    p = Progress(size(data, 1); enabled = show_progress)
    @maybe_threads multithreading for ch = 1:size(data, 1)
        for t = 1:size(data, 2)
            #            @debug("$(ndims(data,)),$t,$ch")

            dd = view(data, ch, t, :)
            ix = @. !ismissing(dd)

            beta[ch, t, :] = @view(X[ix, :]) \ @view(data[ch, t, ix])
            # qr(X) was slower on Februar 2022
        end
        next!(p)
    end
    finish!(p)
    if stderror
        stderror = calculate_stderror(X, data, beta)
        modelfit = LinearModelFit(beta, ["solver_default"], stderror)
    else
        modelfit = LinearModelFit(beta, ["solver_default"])
    end
    return modelfit
end



