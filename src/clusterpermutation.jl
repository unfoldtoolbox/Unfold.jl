
using Random
using ProgressMeter
using PyMNE # not in Unfold.jl currently
using MixedModelsPermutations
#using BlockDiagonals # for olsranef
using StatsModels  # for olsranef

function cluster_permutation_test(mres::UnfoldLinearMixedModel,
                    dat::Array,
                    times::StepRangeLen,
                    coeffOfInterest = 2;
                    clusterFormingThreshold = 2,
                    adjacency = nothing,
                    nPerm = 500,
                    test_times = (minimum(times),maximum(times))

                    )



# from where to where to actually test (can be subrange)
@assert(test_times[1]<test_times[2])
tRange = findfirst(times.>=test_times[1]):findfirst(times.>=test_times[2])


# extract observed z-values
z_obs = [m.z for m in coefpvalues(mres.modelinfo) if String(m.coefname)==coefnames(mres.X.formulas)[2][coeffOfInterest]]

# cluster observed
obs_cluster = pymne_cluster(z_obs,clusterFormingThreshold;tRange = tRange,adjacency=adjacency)

# if there is no cluster in the real data, return empty p_df
if isempty(obs_cluster[2])
    p_df = DataFrame(:from=>[],:to=>[],:pval=>[],:coefname=>[])
    return p_df
end

# this is the meat of the function.
permDat = cluster_permutation(mres,dat,tRange,coeffOfInterest,nPerm)


# cluster permuted
perm_cluster = [pymne_cluster(permDat[:,p],clusterFormingThreshold;adjacency=adjacency) for p = 1:nPerm]

# get maximal cluster value
perm_H0 = maximum.([abs.(p[2]) for p in perm_cluster if ~isempty(p[2])])
perm_H0 = vcat(perm_H0,zeros(nPerm-length(perm_H0)))

# get pvalues
p_vals = PyMNE.stats.cluster_level._pval_from_histogram(obs_cluster[2], perm_H0, 0)

fs = 1 ./Float64(times.step)
# get p-values, we have to shift by the test_times[1] though
if clusterFormingThreshold == "tfce"
    # not sure how to handle this nicer...
    fromList = [o.start ./ fs - test_times[1] for o in obs_cluster[1]]
    toList = [o.stop ./fs - test_times[1] for o in obs_cluster[1]]
else
    fromList = [o[1].start ./ fs - test_times[1] for o in obs_cluster[1]]
toList = [o[1].stop ./fs - test_times[1] for o in obs_cluster[1]]
end
p_df = DataFrame(:from=>fromList,:to=>toList,:pval=>p_vals)

p_df[!,:coefname] .= coefnames(mres.X.formulas)[2][coeffOfInterest] # yeah this might break :S
return p_df
#println(p_vals) 
#h = hist(perm_H0,bins=100)
#vlines!(h.axis,abs.(obs_cluster[2]),color="red")

end



cluster_permutation(args...;kwargs...) = cluster_permutation(MersenneTwister(1),args...;kwargs...)

function cluster_permutation(rng::AbstractRNG,mres,dat,tRange,coeffOfInterest,nPerm)
    permDat = Matrix{Float64}(undef,length(tRange),nPerm)
    mm_outer = Unfold.LinearMixedModel_wrapper(mres.X.formulas,dat[1,1,:],mres.X.Xs)
    
    chIx = 1 # for now we only support 1 channel anyway
    #
    #p = Progress(length(tRange))
    #Threads.@threads for tIx =1:length(tRange)
    #@showprogress "Processing Timepoints" 
    for tIx =1:length(tRange)
        println(string(tIx))
        # splice in the correct data for residual calculation
        mm = deepcopy(mm_outer)
        mm.y .= dat[chIx,tRange[tIx],:]
        
        # set the previous calculated model-fit
        updateL!(setθ!(mm,Vector(mres.modelinfo.θ[tRange[tIx]])))
    
        # get the coefficient 
        H0 = coef(mm)
        # set the one of interest to 0
        H0[coeffOfInterest] = 0
        # run the permutation
        # important here is to set the same seed to keep flip all time-points the same
        perm = permutation(deepcopy(rng),nPerm,mm;β=H0,blup_method=olsranef,use_threads=false); # constant rng to keep autocorr & olsranef for singular models
        
        # extract the z-value
        perm_z = [m.z for m in perm.coefpvalues if String(m.coefname)==coefnames(mm)[coeffOfInterest]]
    
        # save it
        permDat[tIx,:] = perm_z
        #next!(p)
    end
    return permDat
    end


    # function to call pymne -> _find_cluster function
    function pymne_cluster(data,clusterFormingThreshold;tRange = 1:size(data,1), adjacency=nothing)

        return PyMNE.stats.cluster_level._find_clusters(data,clusterFormingThreshold,adjacency=adjacency,include=tRange)
    end
    
    pymne_cluster(data,clusterFormingThreshold::String;kwargs...) = pymne_cluster(data,Dict(:start=>0,:step=>0.2);kwargs...)



function olsranefjf(model::LinearMixedModel{T}) where {T}
    
    l = size(model.reterms)[1]

    mat = Array{Any}(undef, l);
    code = Array{Any}(undef, l);
    ### I get the contrasts 
    for i in 1:l
        trm = model.reterms[i];
        dim = size(trm.z)[1];
        cd = StatsModels.ContrastsMatrix(EffectsCoding(), trm.levels).matrix;
        cd = kron(cd,I(dim));
        code[i] = cd;
        mat[i] = trm*cd;
    end
    mat
    X = hcat(mat...)
    X1 = hcat(ones(size(X)[1]), X)  
    fixef_res = response(model) - model.X*model.β;
    flatblups = X1'X1 \ X1'fixef_res;
    flatblups = deleteat!(flatblups,1)
    
    code_all = BlockDiagonal([x for x in code])
    
    flatblups = code_all*flatblups

    blups = Vector{Matrix{T}}()

    offset = 1
    for trm in model.reterms
        chunksize = size(trm, 2)
        ngrps = length(trm.levels)
        npreds = length(trm.cnames)
        re = Matrix{T}(reshape(view(flatblups, offset:(offset+chunksize-1)),
                               npreds, ngrps))
        offset += chunksize
        push!(blups, re)
    end
    return blups,dummy_scalings(model.reterms)
end
