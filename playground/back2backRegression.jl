using unfold, DataFrames
using StatsModels
using ProgressMeter
using CSV
using Tables
using Statistics
using DataFramesMeta
using MLBase
using LinearAlgebra
include("../test/debug_readEEGlab.jl")
include("dataset-CORE_helper.jl")
sub = 1

# load data
filename = "C:/Users/behinger/Downloads/p3_clean/$(sub)_P3_shifted_ds_reref_ucbip_hpfilt_ica_weighted_clean.set"
data,srate,evts_df,chanlocs_df,EEG = import_eeglab(filename)
data = data[1:29,:] # remove EOGs
# reref to average ref
data = data .- mean(data,dims=1)
# clean data
data = unfold.clean_data(data,EEG["reject"]["rejmanual"]) # we marked bad segments in matlab/eeglab
evts= parse_trigger_p3(evts_df)
evts = @where(evts, .!:invalidresponse)
# Mass Univariate
f_stim = @formula 0~0+trialtype
f_butt = @formula 0~0+trialtype#+answer
contrast = Dict(:trialtype => DummyCoding(base="distractor"),:answer=>DummyCoding(base="distractor"))
evts_stim = filter(x->x.eventtype=="stimulus",evts)
data_e,times = unfold.epoch(data=data,tbl=evts_stim,τ=(-0.8,1.3),sfreq=srate)
Xstim = designmatrix(UnfoldLinearModel,f_stim,filter(x->(x.eventtype=="stimulus"),evts))



##

##

using MLJ
function solver_b2bcv(X,data::AbstractArray{T,3},cross_val_reps = 10;solver=(a,b,c)->ridge(a,b,c)) where {T<:Union{Missing, <:Number}}
    
    #
    X,data = unfold.dropMissingEpochs(X,data)

    # Open MLJ Tuner
    @load RidgeRegressor pkg=MLJLinearModels
    ##

    tm = TunedModel(model=RidgeRegressor(fit_intercept=false),
                resampling = CV(nfolds=5),
                tuning = Grid(resolution=10),
                range = range(RidgeRegressor(), :lambda, lower=1e-2, upper=1000, scale=:log10),
                measure = rms)
    

    E = zeros(size(data,2),size(X,2),size(X,2))
    W = Array{Float64}(undef,size(data,2),size(X,2),size(data,1))
    println("n = samples = $(size(X,1)) = $(size(data,3))")
    @showprogress 0.1 for t in 1:size(data,2)        
        for m in 1:cross_val_reps
            k_ix = collect(Kfold(size(data,3),2))
            Y1 = data[:,t,k_ix[1]]
            Y2 = data[:,t,k_ix[2]]
            X1 = X[k_ix[1],:]
            X2 = X[k_ix[2],:]
       
            G = solver(tm,Y1',X1)
            H = solver(tm,X2, (Y2'*G))

            E[t,:,:] = E[t,:,:]+Diagonal(H[diagind(H)])

        end
        E[t,:,:] = E[t,:,:] ./ cross_val_reps
        W[t,:,:] = (X*E[t,:,:])' / data[:,t,:]

    end

    # extract diagonal
    beta = mapslices(diag,E,dims=[2,3])
    # reshape to conform to ch x time x pred
    beta = permutedims(beta,[3 1 2])
    modelinfo = Dict("W"=>W,"E"=>E,"cross_val_reps"=>cross_val_reps) # no history implemented (yet?)
    return beta, modelinfo
end



function ridge(tm,data,X)
    G = Array{Float64}(undef,size(data,2),size(X,2))
    for pred in 1:size(X,2)
        #println(elscitype(data))
        mtm = machine(tm,table(data),X[:,pred])
        fit!(mtm,verbosity=0)
        G[:,pred] = Tables.matrix(fitted_params(mtm).best_fitted_params.coefs)[:,2]
    end
    return G
end

using GLMNet
function ridge_glmnet(tm,data,X)
    G = Array{Float64}(undef,size(data,2),size(X,2))
    for pred in 1:size(X,2)
        #println(elscitype(data))
        cv = glmnetcv(data,X[:,pred],intercept=false)
        G[:,pred] =cv.path.betas[:,argmin(cv.meanloss)]
    end
    return G
end

import ScikitLearn
using ScikitLearn.GridSearch: GridSearchCV

@ScikitLearn.sk_import linear_model: Ridge
function ridge_sklearn(tm,data,X)
    G = Array{Float64}(undef,size(data,2),size(X,2))
    D = Dict(:C => 10 .^range(log10(tm.range.lower),stop=log10(tm.range.upper),length=10))

    cv = GridSearchCV(Ridge(),Dict(:alpha => (10 .^range(log10(tm.range.lower),stop=log10(tm.range.upper),length=10))))
    ScikitLearn.fit!(cv,data,X)

    G = cv.best_estimator_.coef_'
    
    return G
    

end

##

##
@time b,modelinfo = unfold.solver_b2b(Xstim.Xs,data_e); #1.2s
@time b_cv,modelinfo_cv = solver_b2bcv(Xstim.Xs,data_e); #~9min
#@time b_cvglmnet,modelinfo_cvglmnet = solver_b2bcv(Xstim.Xs,data_e;solver=(a,b,c)->ridge_glmnet(a,b,c)) # ~21 min
@time b_cvsk,modelinfo_cvsk = solver_b2bcv(Xstim.Xs,data_e;solver=(a,b,c)->ridge_sklearn(a,b,c)); # ~6min




##
um,res_tmp = fit(UnfoldLinearModel,f_stim, evts_stim,data_e,times,solver=(a,b)->unfold.solver_b2b(a,b))

###

X,data = unfold.dropMissingEpochs(Xstim.Xs,data_e)

D = Dict(:C => 10 .^range(log10(tm.range.lower),stop=log10(tm.range.upper),length=10))

cv = GridSearchCV(Ridge(fit_intercept=false,),Dict(:alpha => (10 .^range(log10(tm.range.lower),stop=log10(tm.range.upper),length=10))))

ScikitLearn.fit!(cv, data[:,1,:]', X)

##
X,data_m = unfold.dropMissingEpochs(Xstim.Xs,data_e)

using MLBase
k_ix = collect(Kfold(size(data_m,3),2))
Y1 = data_m[:,1,k_ix[1]]
Y2 = data_m[:,1,k_ix[2]]
X1 = X[k_ix[1],:]
X2 = X[k_ix[2],:]

XL = collect(Y1')
y = X1[:,1]
