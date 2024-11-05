using UnfoldSim
using Krylov, CUDA
using TimerOutputs, ProgressMeter
using Unfold
using Random
using DataFrames
using StableRNGs

using BSplineKit
using Chairmarks
using LinearAlgebra
using CUDA.CUSPARSE
using SparseArrays

#--- Generate Data


function benchmark_data(;
    sfreq = 100,
    n_repeats = 100,
    n_splines = 10,
    n_channels = 50,
    overlap = (0.2, 0.2),
)
    data, evts = UnfoldSim.predef_eeg(StableRNG(1); n_repeats, sfreq)
    evts = evts[1:end-2, :]
    data = reshape(data, 1, :)
    evts.type = rand(StableRNG(1), [0, 1], nrow(evts))

    ba1 = firbasis(τ = (-0.1, 1), sfreq = sfreq)
    ba2 = firbasis(τ = (-0.3, 1), sfreq = sfreq)
    if n_splines == 0
        f1 = @formula 0 ~ 1 + condition
        f2 = @formula 0 ~ 1 + condition
    elseif isa(n_splines, Int)
        f1 = @eval @formula 0 ~ 1 + spl(continuous, $n_splines) + condition
        f2 = @eval @formula 0 ~ 1 + spl(continuous, $n_splines) + condition
    elseif isa(n_splines, Tuple)
        @assert length(n_splines) >= 2
        f1 = @eval @formula 0 ~ 1 + condition + spl(continuous, $(n_splines[1]))
        f2 = @eval @formula 0 ~ 1 + condition + spl(continuous, $(n_splines[1]))
        s_basic = f1.rhs[3]
        k = 1
        for s in n_splines[2:end]
            k = k + 1
            spl_name = Symbol("cont_$k")
            evts[:, spl_name] .= rand(MersenneTwister(k), size(evts, 1))
            s_basic.args[1] = Term(spl_name)
            s_basic.args[2] = ConstantTerm(s)
            #            s_basic.exorig = 
            s_basic = FunctionTerm(s_basic.f, s_basic.args, :(spl($spl_name, $s)))
            f1 = FormulaTerm(f1.lhs, f1.rhs + s_basic)
            f2 = FormulaTerm(f2.lhs, f2.rhs + s_basic)
            @info f1
        end

    end

    dict_lin = [0 => (f1, ba1), 1 => (f2, ba2)]

    X = modelmatrix(
        designmatrix(UnfoldLinearModelContinuousTime, dict_lin, evts; eventcolumn = :type),
    )

    data_one = data[1:1, 1:size(X, 1)] # cute the data to have same length
    data20 = repeat(data_one, n_channels)
    data20 .= data20 .+ rand(StableRNG(1), size(data20)...) .* 20

    return X, data20
end

function unfold_benchmarks(;
    n_channels,
    sfreq,
    n_repeats,
    n_splines,
    max_seconds = 1000,
    overlap = (0.5, 0.2),
)
    X, y = benchmark_data(; n_channels, sfreq, n_repeats, n_splines, overlap)

    #y = b_data[1:1, :]
    #max_seconds = 10

    #default_single =
    #    @be Unfold.solver_default(X, y; multithreading = false) seconds = max_seconds
    default_multi =
        @be Unfold.solver_default(X, y; multithreading = true) seconds = max_seconds
    #krylov_single =
    #    @be Unfold.solver_krylov(X, y; multithreading = false) seconds = max_seconds
    #krylov_multi = @be Unfold.solver_krylov(X, y; multithreading = true) seconds = max_seconds
    krylov_gpu = @be Unfold.solver_krylov(X, y; GPU = true) seconds = max_seconds # multithreading = false automatically
    df_res = DataFrame(
        [
            #("default_single", default_single, false, minimum(default_single)),
            ("default_multi", default_multi, false, minimum(default_multi)),
            #       ("krylov_single", krylov_single, false, minimum(krylov_single)),
            #("krylov_multi", minimum(krylov_multi)),
            ("krylov_gpu", krylov_gpu, true, minimum(krylov_gpu)),
        ],
        [:method, :res, :gpu, :min],
    )

    #df_res = DataFrame()
    res = []
    for gpu in [false, true]
        for s in [:cg, :pinv, :intern, :qr]#, :krylov_cg]
            @info "solver $s - gpu:$gpu"
            #if gpu == true && s == :pinv
            # scalar indexing issue
            #    continue
            #end
            try
                y_solver = gpu ? cu(y) : y
                res = @be _ Unfold.solver_predefined(
                    X,
                    y_solver;
                    solver = s,
                    multithreading = false,
                ) begin
                    GC.gc()
                    CUDA.reclaim()
                end seconds = max_seconds
            catch err
                @info "ERROR!!!" typeof(err)
                continue
            end
            df_res = vcat(
                df_res,
                DataFrame(:method => s, :gpu => gpu, :res => res, :min => minimum(res)),
            )

        end
    end
    df_res.min_time = map(x -> x.time, df_res.min)
    df_res.min_allocs = map(x -> x.allocs, df_res.min)
    df_res.min_bytes = map(x -> x.bytes, df_res.min)
    df_res.sizeDesign .= Ref(size(X))
    df_res.sizeEEG .= Ref(size(y))
    df_res.n_channels .= n_channels
    df_res.sfreq .= sfreq
    df_res.n_repeats .= n_repeats
    df_res.n_splines .= Ref(n_splines)
    df_res.overlap .= Ref(overlap)
    df_res.nnz .= nnz(X)
    sort!(df_res, :min_time)
    return df_res
end

function XX_solver(X, data::AbstractArray{T}; solver = :cg, gpu = false) where {T}

    if gpu
        X = CuSparseMatrixCSC{T}(X)

        Xt = CuSparseMatrixCSC(X')
        R_xx = CuArray{T}(X'X) # no longer sparse


        Y = CuArray{T}(data)
        R_xy = CuVector{T}(undef, size(X, 2))
        Ĥ = CUDA.zeros(T, size(Y, 1), size(X, 2))
        cg_solver = Krylov.CgSolver(size(R_xx)..., CuVector{T})
    else
        R_xx = Matrix(X'X)
        Y = data
        Xt = Unfold.SparseArrays.SparseMatrixCSC(X') # this didnt give the expected boost, but in older Julia versions, X'y was not very performant. Might be interesting for really large matrices?
        R_xy = Vector{T}(undef, size(X, 2))
        Ĥ = zeros(T, size(Y, 1), size(X, 2))
        cg_solver = Krylov.CgSolver(size(R_xx)..., Vector{T})

    end
    if solver == :pinv
        R_xx_inv = pinv(R_xx)
    elseif solver == :qr
        R_xx_qr = qr(R_xx)
    end
    for ch = 1:size(Y, 1)

        # using a view here I get a sever punishment in speed -really crazy
        #y = @view(Y[ch,:])

        # I dont think the speed improveoment of replacing this with the inplace mul! version is thaaaat strong, but I have it working now, so it certainly won't hurt ;)
        #R_xy .= Xt * Y[ch, :]
        mul!(R_xy, Xt, Y[ch, :])

        # set a more reasonable starting point
        ch == 1 || copyto!(view(Ĥ, ch, :), view(Ĥ, ch - 1, :))

        if solver == :cg
            _, h = Unfold.IterativeSolvers.cg!(@view(Ĥ[ch, :]), R_xx, R_xy, log = true)
        elseif solver == :krylov_cg
            Krylov.solve!(cg_solver, R_xx, R_xy; history = true)
            Ĥ[ch, :] .= cg_solver.x

        elseif solver == :intern
            Ĥ[ch, :] .= R_xx \ R_xy
        elseif solver == :pinv
            Ĥ[ch, :] .= R_xx_inv * R_xy
        elseif solver == :qr
            Ĥ[ch, :] .= R_xx_qr \ R_xy
        else
            error("not implemented")
        end
    end
    return Ĥ

end

#---

case_small = unfold_benchmarks(
    n_channels = 1,
    sfreq = 10,
    n_splines = 4,
    n_repeats = 10;
    max_seconds = 1,
)
case_default = unfold_benchmarks(
    n_channels = 1,
    sfreq = 100,
    n_splines = 4,
    n_repeats = 200;
    max_seconds = 5,
)
case_sfreq = unfold_benchmarks(
    n_channels = 1,
    sfreq = 1000,
    n_splines = 4,
    n_repeats = 200;
    max_seconds = 10,
)
case_tall = unfold_benchmarks(
    n_channels = 1,
    sfreq = 100,
    n_splines = 4,
    n_repeats = 2000;
    max_seconds = 1,
)
#case_wide = unfold_benchmarks(
#    n_channels = 1,
#    sfreq = 100,
#    n_splines = (10, 10, 10, 10, 10),
#    n_repeats = 400;
#    max_seconds = 1,
#)
case_multichannel = unfold_benchmarks(
    n_channels = 128,
    sfreq = 100,
    n_splines = 4,
    n_repeats = 200,
    max_seconds = 10,
)

case_large = unfold_benchmarks(
    n_channels = 128,
    sfreq = 500,
    n_splines = (4, 4),
    n_repeats = 500,
    max_seconds = 10,
)


case_wide = unfold_benchmarks(
    n_channels = 128,
    sfreq = 100,
    n_splines = 4,
    n_repeats = 200,
    overlap = (0.5, 0.2),
    max_seconds = 2,
)


#---
X, y = benchmark_data(
    n_channels = 10,
    sfreq = 100,
    n_splines = 4,
    n_repeats = 100,
    overlap = (0.2, 0.2),
)


Unfold.solver_predefined(X, cu(y); solver = :qr, show_time = true, multithreading = false)

#function solver_cg_krylov!(beta, data, Xt, R_xx, R_xy)
#    calc_Rxy!(R_xy, Xt, data)
#    _, h = Unfold.IterativeSolvers.cg!(beta, R_xx, R_xy, log = true)
#end
