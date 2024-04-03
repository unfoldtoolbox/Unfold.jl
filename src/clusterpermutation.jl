using Logging
using Random
using ProgressMeter

function cluster_permutation_test(mres, dat, times, coeffOfInterest::Array; kwargs...)
    df_p = DataFrame()
    for coef in coeffOfInterest
        df_p_single = cluster_permutation_test(mres, dat, times, coef; kwargs...)
        append!(df_p, df_p_single)
    end
    return df_p
end


function cluster_permutation_test(
    mres::UnfoldLinearMixedModel,
    dat::Array,
    times::StepRangeLen,
    coeffOfInterest::Int;
    clusterFormingThreshold = 2,
    adjacency = nothing,
    nPerm = 500,
    test_times = (minimum(times), maximum(times)),
    stat = :β,
)



    # from where to where to actually test (can be subrange)
    @assert(test_times[1] < test_times[2])
    tRange = findfirst(times .>= test_times[1]):findfirst(times .>= test_times[2])

    # extract observed z-values
    stat_obs = [
        getproperty(m, stat) for m in coefpvalues(mres.modelfit) if
        String(m.coefname) == coefnames(formula(designmatrix(mres)))[2][coeffOfInterest]
    ]


    # cluster observed
    obs_cluster = pymne_cluster(
        stat_obs,
        clusterFormingThreshold;
        tRange = tRange,
        adjacency = adjacency,
    )

    # if there is no cluster in the real data, return empty p_df
    if isempty(obs_cluster[2])
        p_df = DataFrame(:from => [], :to => [], :pval => [], :coefname => [])
        return p_df
    end

    # this is the meat of the function.
    permDat = cluster_permutation(mres, dat, tRange, coeffOfInterest, nPerm; stat = stat)


    # cluster permuted
    perm_cluster = [
        pymne_cluster(permDat[:, p], clusterFormingThreshold; adjacency = adjacency) for
        p = 1:nPerm
    ]

    # get maximal cluster value
    perm_H0 = maximum.([abs.(p[2]) for p in perm_cluster if ~isempty(p[2])])
    perm_H0 = vcat(perm_H0, zeros(nPerm - length(perm_H0)))

    # get pvalues
    p_vals = PyMNE.stats.cluster_level._pval_from_histogram(obs_cluster[2], perm_H0, 0)

    fs = 1 ./ Float64(times.step)

    # get p-values, we have to shift by the test_times[1] though
    if clusterFormingThreshold == "tfce"
        # not sure how to handle this nicer...
        fromList = [o.start ./ fs + times[1] for o in obs_cluster[1]]
        toList = [o.stop ./ fs + times[1] for o in obs_cluster[1]]
    else
        fromList = [o[1].start ./ fs + times[1] for o in obs_cluster[1]]
        toList = [o[1].stop ./ fs + times[1] for o in obs_cluster[1]]
    end
    p_df = DataFrame(:from => fromList, :to => toList, :pval => p_vals)

    p_df[!, :coefname] .= coefnames(formula(mres))[2][coeffOfInterest] # yeah this might break :S
    return p_df
    #println(p_vals) 
    #h = hist(perm_H0,bins=100)
    #vlines!(h.axis,abs.(obs_cluster[2]),color="red")

end



cluster_permutation(args...; kwargs...) =
    cluster_permutation(MersenneTwister(1), args...; kwargs...)

function cluster_permutation(
    rng::AbstractRNG,
    mres,
    dat,
    tRange,
    coeffOfInterest,
    nPerm;
    stat = :β,
)
    permDat = Matrix{Float64}(undef, length(tRange), nPerm)
    mm_outer = Unfold.LinearMixedModel_wrapper(formula(mres), dat[1, 1, :], modelmatrix(mres))
    mm_outer.optsum.maxtime = 0.01 # 

    chIx = 1 # for now we only support 1 channel anyway
    #
    #p = Progress(length(tRange))
    #Threads.@threads for tIx =1:length(tRange)
    @showprogress "Processing Timepoints" for tIx = 1:length(tRange)
        #println(string(tIx))
        # splice in the correct data for residual calculation
        mm = deepcopy(mm_outer)
        mm.y .= dat[chIx, tRange[tIx], :]

        # set the previous calculated model-fit
        updateL!(setθ!(mm, Vector(mres.modelfit.θ[tRange[tIx]])))

        # get the coefficient 
        H0 = coef(mm)
        # set the one of interest to 0
        H0[coeffOfInterest] = 0
        # run the permutation
        
        perm = undef
        Logging.with_logger(NullLogger()) do   # remove NLopt warnings  
         perm = permutation(
            deepcopy(rng), # important here is to set the same seed to keep flip all time-points the same
            nPerm,
            mm;
            β = H0,
            blup_method = olsranef,
            use_threads = false,
            hide_progress = true,
            inflation_method= (x,y,z) ->I(length(H0)),#MixedModelsPermutations.inflation_factor(mm, olsranef(mm), residuals(mm))

        ) # constant rng to keep autocorr & olsranef for singular models
        end

        # extract the test-statistic

        #perm_z = [m.z for m in perm.coefpvalues if String(m.coefname)==coefnames(mm)[coeffOfInterest]]
        permDat[tIx, :] = [
            getproperty(m, stat) for m in perm.coefpvalues if
            String(m.coefname) == coefnames(mm)[coeffOfInterest]
        ]

        #next!(p)
    end
    return permDat
end


# function to call pymne -> _find_cluster function
function pymne_cluster(
    data,
    clusterFormingThreshold;
    tRange = 1:size(data, 1),
    adjacency = nothing,
)


    #return PyMNE.stats.cluster_level._find_clusters(data,clusterFormingThreshold,adjacency=adjacency,include=tRange)
    # this raised an error due to the named arguments, but was working before. not sure whats going on

    # need to convert to boolean
    tRange_bool = zeros(size(data, 1))
    tRange_bool[tRange] .= 1
    return PyMNE.stats.cluster_level._find_clusters(
        data,
        clusterFormingThreshold,
        0,
        adjacency,
        1,
        tRange_bool,
    )
end

pymne_cluster(data, clusterFormingThreshold::String; kwargs...) =
    pymne_cluster(data, Dict(:start => 0, :step => 0.2); kwargs...)
