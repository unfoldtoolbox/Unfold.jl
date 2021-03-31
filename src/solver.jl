using StatsBase: var
function solver_default(X,data::AbstractArray{T,2}) where {T<:Union{Missing, <:Number}}
    # likely much larger matrix, using lsqr
    modelinfo = []
    beta = Array{Float64}(undef,size(data,1),size(X,2))
   
    @showprogress .1 for ch in 1:size(data,1)
        ix = .!ismissing.(data[ch,:])
        beta[ch,:],h = lsmr(X[ix,:],data[ch,ix],log=true)
        push!(modelinfo,h)   
    end
   
    return beta, modelinfo
end
struct lsmr_modelinfo
    history
    stderror
end
function solver_lsmr_se(X,data::AbstractArray{T,2}) where {T<:Union{Missing, <:Number}}
    # likely much larger matrix, using lsqr
    modelinfo = []
    beta = Array{Float64}(undef,size(data,1),size(X,2))
   
    @showprogress .1 for ch in 1:size(data,1)
        ix = .!ismissing.(data[ch,:])
        beta[ch,:],h = lsmr(X[ix,:],data[ch,ix],log=true)
        push!(modelinfo,h)   
    end
    ix = any(.!ismissing.(data),dims=1)[1,:]
    stderror = calculate_stderror(X[ix,:],data[:,ix],beta)
    return beta, lsmr_modelinfo(modelinfo,stderror)
end
function calculate_stderror(Xdc,data,beta)    
    # Hat matrix
    hat_prime = inv(Matrix( Xdc'*Xdc))
    # Calculate residual variance
    @warn("Autocorrelation was NOT taken into account. Therefore SE are UNRELIABLE. Use at your own discretion")
    
    se = Array{Float64}(undef,size(data,1),size(Xdc,2))
    for ch = 1:size(data,1)
        residualVar = var(data[ch,:] .- Xdc*beta[ch,:]);
        @assert(!isnan(residualVar),"residual Variance was NaN")
        hat = hat_prime .* residualVar
        #se = sqrt(diag(cfg.contrast(:,:)*hat*cfg.contrast(:,:)'));
        se[ch,:] = diag(hat)
    end
    
    return se
    
    
end
function solver_default(X,data::AbstractArray{T,3}) where {T<:Union{Missing, <:Number}}
       beta = Array{Union{Missing,Number}}(undef,size(data,1),size(data,2),size(X,2))
       @showprogress .1 for ch in 1:size(data,1)
           for t in 1:size(data,2)
               @debug("$(ndims(data,)),$t,$ch")
               ix = .!ismissing.(data[ch,t,:])
               beta[ch,t,:] = X[ix,:] \ data[ch,t,ix]
           end
       end
   
       modelinfo = [undef] # no history implemented (yet?)
       return beta, modelinfo
end

function solver_b2b(X,data::AbstractArray{T,3},cross_val_reps = 10) where {T<:Union{Missing, <:Number}}
    
    X,data = dropMissingEpochs(X,data)


    E = zeros(size(data,2),size(X,2),size(X,2))
    W = Array{Float64}(undef,size(data,2),size(X,2),size(data,1))
    println("n = samples = $(size(X,1)) = $(size(data,3))")
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

