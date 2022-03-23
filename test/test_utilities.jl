using CSV
using DelimitedFiles
using DSP
using Random
using LinearAlgebra
using DataFrames

function loadtestdata(
    testCase::String;
    dataPath::String = (@__DIR__) * "/data_new_testcases",
)
    #println(pwd()) # to debug github action
    data = readdlm(joinpath(dataPath, "$(testCase)_data.csv"), ',', Float64, '\n')
    data = dropdims(data, dims = 1) # convert to vector
    evts = CSV.read(joinpath(dataPath, "$(testCase)_events.csv"), DataFrame)
    return data, evts
end


function gridexpand(conditionA = [0, 1.0], continuousA = [-1.0, 0, 1.0])

    tmp = reshape(
        [[x, y] for x in conditionA, y in continuousA],
        length(conditionA) * length(continuousA),
    )
    evts_grid = DataFrame(hcat(tmp...)')
    rename!(evts_grid, ["conditionA", "continuousA"])
    return evts_grid
end




simulate_lmm(args...; kwargs...) = simulate_lmm(Random.GLOBAL_RNG, args...; kwargs...)

function simulate_lmm(
    rng::AbstractRNG,
    τ = 1.5,
    fs = 12;
    β = [0.0, -1.0],
    σs = [1, 1, 1],
    σ = 1,
    n_sub = 20,
    n_item = 30,
    noise_type = "AR-exponential",
)
    rng_copy = deepcopy(rng)


    subj_btwn = item_btwn = both_win = nothing
    #subj_btwn = Dict("age" => ["O", "Y"])

    # there are no between-item factors in this design so you can omit it or set it to nothing
    item_btwn = Dict("stimType" => ["I", "II"])

    # put within-subject/item factors in a Dict
    #both_win = Dict("condition" => ["A", "B"])

    # simulate data
    evt = DataFrame(
        simdat_crossed(
            n_sub,
            n_item,
            subj_btwn = subj_btwn,
            item_btwn = item_btwn,
            both_win = both_win,
        ),
    )

    #    f1 = @formula dv ~ 1 + age * condition  + (1+condition|item) + (1+condition|subj);
    f1 = @formula dv ~ 1 + stimType + (1 + stimType | subj) + (1 | item)
    m = MixedModels.fit(MixedModel, f1, evt)

    # set the random effects

    #gen_han(τ,fs,1)
    basis = gen_han(τ, fs, 2)

    epoch_dat = zeros(Int(τ * fs), size(evt, 1))

    #MixedModels doesnt really support σ==0 because all Ranef are scaled residual variance
    σ_lmm = 0.0001
    σs = σs ./ σ_lmm

    for t = 1:size(epoch_dat, 1)
        b = basis[t]

        MixedModelsSim.update!(m, create_re(b .* σs[1], b .* σs[2]), create_re(b .* σs[3]))
        simulate!(deepcopy(rng_copy), m, β = [b .* β[1], b .* β[2]], σ = σ_lmm)


        epoch_dat[t, :] = m.y
    end

    epoch_dat = reshape(epoch_dat, (1, size(epoch_dat)...))

    # add some noise
    if noise_type == "normal"
        epoch_dat = epoch_dat .+ randn(rng, size(epoch_dat)) .* σ
    elseif noise_type == "AR-exponential"
        epoch_dat = epoch_dat .+ gen_noise_exp(rng, size(epoch_dat)) .* σ

    end

    return evt, epoch_dat
end


function gen_han(τ, fs, peak)
    hanLen = Int(τ * fs / 3)
    han = hanning(hanLen, zerophase = false)
    sig = zeros(Int(τ * fs))
    sig[1+hanLen*(peak-1):hanLen*peak] .= han
    return sig
end





function circulant(x)
    # Author: Jaromil Frossard
    # returns a symmetric matrix where X was circ-shifted.
    lx = length(x)
    ids = [1:1:(lx-1);]
    a = Array{Float64,2}(undef, lx, lx)
    for i = 1:length(x)
        if i == 1
            a[i, :] = x
        else
            a[i, :] = vcat(x[i], a[i-1, ids])
        end
    end
    return Symmetric(a)
end


function exponentialCorrelation(x; nu = 1, length_ratio = 1)
    # Author: Jaromil Frossard
    # generate exponential function
    R = length(x) * length_ratio
    return exp.(-3 * (x / R) .^ nu)
end

function noise_exp(rng, n_t)
    Σ = circulant(exponentialCorrelation([0:1:(n_t-1);], nu = 1.5))
    Ut = LinearAlgebra.cholesky(Σ).U'
    return (randn(rng, n_t)'*Ut')[1, :]
end

function gen_noise_exp(rng, si)
    n = hcat(map(x -> noise_exp(rng, si[2]), 1:si[3])...)
    return reshape(n, si)
end
