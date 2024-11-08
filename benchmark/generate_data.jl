using Random
using UnfoldSim
using StatsModels
using StableRNGs

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
