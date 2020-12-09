function solver_default(X,data::AbstractArray{T,2}) where {T<:Union{Missing, <:Number}}
    # likely much larger matrix, using lsqr
    @info lsmr
    modelinfo = []
    sizehint!(modelinfo, size(data,1))
    beta = Array{Float64}(undef,size(data,1),size(X,2))

    @showprogress .1 for ch in 1:size(data,1)
        dd = view(data, ch, :)
        ix = @. !ismissing(dd)
        beta[ch,:],h = lsmr!(view(beta, ch, :), view(X, ix, :), view(data, ch, ix),log=true)
        push!(modelinfo,h)
    end

    return beta, modelinfo
end

function solver_default(X,data::AbstractArray{T,3}) where {T<:Union{Missing, <:Number}}
       beta = Array{Union{Missing,Number}}(undef,size(data,1),size(data,2),size(X,2))
       @showprogress .1 for ch in 1:size(data,1)
           for t in 1:size(data,2)
               @debug("$(ndims(data,)),$t,$ch")
               dd = view(data, ch, t,:)
               ix = @. !ismissing(dd)
               # sparse qr! isn't implemented, so views don't save us allocations
               beta[ch,t,:] = X[ix, :] \ data[ch,t,ix]
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

