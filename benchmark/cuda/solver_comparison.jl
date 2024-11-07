using DataFrames: Markdown, MarkdownHighlighter, MarkdownDecoration
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
include("../generate_data.jl")

function unfold_benchmarks(;
    n_channels,
    sfreq,
    n_repeats,
    n_splines,
    max_seconds = 1000,
    overlap = (0.2, 0.2),
)
    X, y = benchmark_data(; n_channels, sfreq, n_repeats, n_splines, overlap)

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
    df_res.comment .= ""
    #df_res = DataFrame()
    res = []
    for gpu in [false, true]
        for s in [:cg, :pinv, :intern, :qr, :cholesky]#, :krylov_cg]
            @info "solver $s - gpu:$gpu"
            comment = ""
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
                @info "error!!" "gpu:$gpu solver:$s"
                @info "ERROR!!!" typeof(err)
                comment = string(err)
                #continue
                res = @be sum([1, 2])
            end
            df_res = vcat(
                df_res,
                DataFrame(
                    :method => s,
                    :gpu => gpu,
                    :res => res,
                    :min => minimum(res),
                    :comment => comment,
                ),
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
#=case_sfreq = unfold_benchmarks(
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
    =#
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



#---- print markdown

function df_to_md(data)
    data_subset = data[
        :,
        [
            :method,
            :gpu,
            :min_time,
            :min_bytes,
            :sizeDesign,
            :n_channels,
            :sfreq,
            :n_repeats,
            :n_splines,
            :overlap,
            :nnz,
            :comment,
        ],
    ]
    data_subset.min_bytes = data_subset.min_bytes / 1024 / 1024 / 1024
    data_subset.percent_X_filled = data_subset.nnz ./ prod.(data_subset.sizeDesign)
    rename!(data_subset, :min_bytes => :GB, :min_time => :time)

    data_subset = data_subset[:, Not([:nnz, :sfreq, :n_repeats, :n_splines])]
    hl_time = MarkdownHighlighter(
        (d, i, j) ->
            isa(d[i, j], Float64) &&
                (d[i, j] ≈ minimum(data_subset.time[data_subset[:, :comment].==""])),
        MarkdownDecoration(bold = true),
    )
    hl_alloc = MarkdownHighlighter(
        (d, i, j) ->
            isa(d[i, j], Float64) &&
                (d[i, j] ≈ minimum(data_subset.GB[data_subset[:, :comment].==""])),
        MarkdownDecoration(bold = true),
    )
    pretty_table(
        data_subset,
        backend = Val(:markdown),
        header = names(data_subset),
        highlighters = (hl_alloc, hl_time),
    )
end


df_to_md.([case_small, case_default, case_multichannel, case_large])