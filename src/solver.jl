using StatsBase: var
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
    @maybe_threads multithreading for ch = 1:size(data, 1)
        dd = view(data, ch, :)
        ix = @. !ismissing(dd)
        # use the previous channel as a starting point
        ch == 1 || copyto!(view(beta, ch, :), view(beta, ch - 1, :))

        beta[ch, :], h =
            lsmr!(@view(beta[ch, :]), (X[ix, :]), @view(data[ch, ix]), log = true)

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
    stderror = false,
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


function calculate_stderror(
    Xdc,
    data::Matrix,
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
    # Calculate residual variance
    @warn(
        "Autocorrelation was NOT taken into account. Therefore SE are UNRELIABLE. Use at your own discretion"
    )

    se = Array{T2}(undef, size(data, 1), size(data, 2), size(X, 2))
    for ch = 1:size(data, 1)
        for t = 1:size(data, 2)
            ix = .!ismissing.(data[ch, t, :])
            residualVar = var(data[ch, t, ix] .- X[ix, :] * beta[ch, t, :])
            @assert(!isnan(residualVar), "residual Variance was NaN")
            hat = hat_prime .* residualVar
            #se = sqrt(diag(cfg.contrast(:,:)*hat*cfg.contrast(:,:)'));
            se[ch, t, :] = sqrt.(diag(hat))
        end
    end
    return se
end

solver_b2b(X, data, cross_val_reps) = solver_b2b(X, data, cross_val_reps = cross_val_reps)
function solver_b2b(
    X,
    data::AbstractArray{T,3};
    cross_val_reps = 10,
    multithreading = true,
    showprogress = true,
) where {T<:Union{Missing,<:Number}}

    X, data = dropMissingEpochs(X, data)


    E = zeros(size(data, 2), size(X, 2), size(X, 2))
    W = Array{Float64}(undef, size(data, 2), size(X, 2), size(data, 1))

    prog = Progress(size(data, 2) * cross_val_reps, 0.1; enabled = showprogress)
    @maybe_threads multithreading for m = 1:cross_val_reps
        k_ix = collect(Kfold(size(data, 3), 2))
        X1 = @view X[k_ix[1], :]
        X2 = @view X[k_ix[2], :]

        for t = 1:size(data, 2)

            Y1 = @view data[:, t, k_ix[1]]
            Y2 = @view data[:, t, k_ix[2]]


            G = (Y1' \ X1)
            H = X2 \ (Y2' * G)

            E[t, :, :] += Diagonal(H[diagind(H)])
            ProgressMeter.next!(prog; showvalues = [(:time, t), (:cross_val_rep, m)])
        end
        E[t, :, :] = E[t, :, :] ./ cross_val_reps
        W[t, :, :] = (X * E[t, :, :])' / data[:, t, :]

    end

    # extract diagonal
    beta = mapslices(diag, E, dims = [2, 3])
    # reshape to conform to ch x time x pred
    beta = permutedims(beta, [3 1 2])
    modelinfo = Dict("W" => W, "E" => E, "cross_val_reps" => cross_val_reps) # no history implemented (yet?)
    return LinearModelFit(beta, modelinfo)
end
