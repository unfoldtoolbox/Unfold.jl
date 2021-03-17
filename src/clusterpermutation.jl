
using Random
using ProgressMeter
using PyMNE # not in Unfold.jl currently
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


# this is the meat of the function.
permDat = cluster_permutation(mres,dat,tRange,coeffOfInterest,nPerm)

# extract observed z-values
z_obs = [m.z for m in coefpvalues(mres.modelinfo) if String(m.coefname)==coefnames(mres.X.formulas)[2][coeffOfInterest]]

# cluster observed
obs_cluster = pymne_cluster(z_obs,clusterFormingThreshold;adjacency=adjacency)

# cluster permuted
perm_cluster = [pymne_cluster(permDat[:,p],clusterFormingThreshold;adjacency=adjacency) for p = 1:nPerm]

# get maximal cluster value
perm_H0 = maximum.([abs.(p[2]) for p in perm_cluster if ~isempty(p[2])])
perm_H0 = vcat(perm_H0,zeros(nPerm-length(perm_H0)))

# get pvalues
p_vals = PyMNE.stats.cluster_level._pval_from_histogram(obs_cluster[2], perm_H0, 0)

fs = 1 ./Float64(times.step)
# get p-values, we have to shift by the test_times[1] though
p_df = DataFrame(:from=>[o[1].start ./ fs - test_times[1] for o in obs_cluster[1]],:to=>[o[1].stop ./fs - test_times[1] for o in obs_cluster[1]],:pval=>p_vals)

p_df[!,:coefname] .= coefnames(mm)[coeffOfInterest]
return p_df
#println(p_vals)
#h = hist(perm_H0,bins=100)
#vlines!(h.axis,abs.(obs_cluster[2]),color="red")

end




function cluster_permutation(mres,dat,tRange,coeffOfInterest,nPerm)
    permDat = Matrix{Float64}(undef,length(tRange),nPerm)
    mm = unfold.LinearMixedModel_wrapper(mres.X.formulas,dat[1,1,:],mres.X.Xs)
    
    chIx = 1 # for now we only support 1 channel anyway
    @showprogress "Permuting over Time" for tIx =1:length(tRange)
        println(tIx)
        # splice in the correct data for residual calculation
        mm.y .= dat[chIx,tRange[tIx],:]
        
        # set the previous calculated model-fit
        updateL!(setθ!(mm,Vector(mres.modelinfo.θ[tRange[tIx]])))
    
        # get the coefficient 
        H0 = coef(mm)
        # set the one of interest to 0
        H0[coeffOfInterest] = 0
        # run the permutation
        # important here is to set the same seed to keep flip all time-points the same
        perm = permutation(MersenneTwister(1),nPerm,mm;β=H0,blup_method=olsranef); # constant rng to keep autocorr & olsranef for singular models
        
        # extract the z-value
        perm_z = [m.z for m in perm.coefpvalues if String(m.coefname)==coefnames(mm)[coeffOfInterest]]
    
        # save it
        permDat[tIx,:] = perm_z
    end
    return permDat
    end


    # function to call pymne -> _find_cluster function
    function pymne_cluster(data,clusterFormingThreshold;adjacency=nothing)
        return PyMNE.stats.cluster_level._find_clusters(data,clusterFormingThreshold,adjacency=adjacency)
    end
    