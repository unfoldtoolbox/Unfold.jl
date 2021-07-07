using StatsBase: var
function solver_default(X,data::AbstractArray{T,2};stderror=false) where {T<:Union{Missing, <:Number}}
    # likely much larger matrix, using lsqr
    minfo = []
    beta = Array{Float64}(undef,size(data,1),size(X,2))
   
    @showprogress .1 for ch in 1:size(data,1)
        ix = .!ismissing.(data[ch,:])
        beta[ch,:],h = lsmr(X[ix,:],data[ch,ix],log=true)
        push!(minfo,h)   
    end

    if stderror
        stderror = calculate_stderror(X,data,beta)
        minfo = modelinfo(modelinfo,stderror)
    end
    return beta, minfo
end
struct modelinfo
    history
    stderror
end

function calculate_stderror(Xdc,data::Matrix{T},beta) where {T<:Union{Missing, <:Number}}  

    # remove missings
    ix = any(.!ismissing.(data),dims=1)[1,:]
    if length(ix)!=size(data,2)
        @warn("Limitation: Missing data are calculated over all channels for standard error")
    end
    data = data[:,ix]
    Xdc = Xdc[ix,:]
    
    # Hat matrix only once
    hat_prime = inv(Matrix( Xdc'*Xdc))
    # Calculate residual variance
    @warn("Autocorrelation was NOT taken into account. Therefore SE are UNRELIABLE. Use at your own discretion")
    
    se = Array{Float64}(undef,size(data,1),size(Xdc,2))
    for ch = 1:size(data,1)
        residualVar = var(data[ch,:] .- Xdc*beta[ch,:]);
        @assert(!isnan(residualVar),"residual Variance was NaN")
        hat = hat_prime .* residualVar
        #se = sqrt(diag(cfg.contrast(:,:)*hat*cfg.contrast(:,:)'));
        se[ch,:] = sqrt.(diag(hat))
    end
    return se
end
function calculate_stderror(X,  data::AbstractArray{T,3},beta) where {T<:Union{Missing, <:Number}}  
#function calculate_stderror(Xdc,data::AbstractArray{T,2},beta) where {T<:Union{Missing, <:Number}}  

    # Hat matrix
    hat_prime = inv(Matrix( X'*X))
    # Calculate residual variance
    @warn("Autocorrelation was NOT taken into account. Therefore SE are UNRELIABLE. Use at your own discretion")

    se = Array{Float64}(undef,size(data,1),size(data,2),size(X,2))
    for ch = 1:size(data,1)
        for t = 1:size(data,2)
            ix = .!ismissing.(data[ch,t,:])
            residualVar = var(data[ch,t,ix] .- X[ix,:]*beta[ch,t,:]);
            @assert(!isnan(residualVar),"residual Variance was NaN")
            hat = hat_prime .* residualVar
            #se = sqrt(diag(cfg.contrast(:,:)*hat*cfg.contrast(:,:)'));
            se[ch,t,:] = sqrt.(diag(hat))
        end
    end
    return se
end
function solver_default(X,data::AbstractArray{T,3};stderror=false) where {T<:Union{Missing, <:Number}}
       beta = Array{Union{Missing,Number}}(undef,size(data,1),size(data,2),size(X,2))
       @showprogress .1 for ch in 1:size(data,1)
           for t in 1:size(data,2)
               @debug("$(ndims(data,)),$t,$ch")
               ix = .!ismissing.(data[ch,t,:])
               beta[ch,t,:] = X[ix,:] \ data[ch,t,ix]
           end
       end
   
       minfo = [undef] # no history implemented (yet?)
       if stderror
            stderror = calculate_stderror(X,data,beta)
            minfo = modelinfo(minfo,stderror)
        end
       return beta, minfo
end

solver_b2b(X,data,cross_val_reps) = solver_b2b(X,data,cross_val_reps = cross_val_reps)
function solver_b2b(X,data::AbstractArray{T,3};cross_val_reps = 10) where {T<:Union{Missing, <:Number}}
    
    X,data = dropMissingEpochs(X,data)


    E = zeros(size(data,2),size(X,2),size(X,2))
    W = Array{Float64}(undef,size(data,2),size(X,2),size(data,1))
    
    prog = Progress(size(data,2)*cross_val_reps,.1)
    for t in 1:size(data,2)        

        for m in 1:cross_val_reps
            k_ix = collect(Kfold(size(data,3),2))
            Y1 = data[:,t,k_ix[1]]
            Y2 = data[:,t,k_ix[2]]
            X1 = X[k_ix[1],:]
            X2 = X[k_ix[2],:]
          

            G = (Y1' \ X1)  
            H = X2 \ (Y2'*G)
            
            E[t,:,:] = E[t,:,:]+Diagonal(H[diagind(H)])
            ProgressMeter.next!(prog; showvalues = [(:time,t), (:cross_val_rep,m)])
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





using MLJ
function solver_b2bcv(X,data::AbstractArray{T,3},cross_val_reps = 10;solver=(a,b,c)->ridge(a,b,c)) where {T<:Union{Missing, <:Number}}
    
    #
    X,data = Unfold.dropMissingEpochs(X,data)

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