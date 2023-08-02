using UnfoldSim
using Krylov,CUDA
using Unfold
using Random
using DataFrames
using StableRNGs
using BSplineKit

#--- Generate Data
sfreq = 100
data,evts = UnfoldSim.predef_eeg(StableRNG(1);n_repeats = 200,sfreq=sfreq)
evts = evts[1:end-2,:]
data = reshape(data,1,:)
evts.type = rand(StableRNG(1),[0,1],nrow(evts))

ba1 = firbasis(τ=(0,1),sfreq = sfreq,name="evts1")
ba2 = firbasis(τ=(0,1),sfreq = sfreq,name="evts2")

f1  = @formula 0~1+spl(continuous,5)+condition
f2  = @formula 0~1+spl(continuous,5)+condition

dict_lin = Dict(0 => (f1,ba1),1 => (f2,ba2))

uf = fit(UnfoldModel,dict_lin,evts,data,eventcolumn="type",solver=(x,y)->nothing);

#----
X = modelmatrix(uf)

data_one = data[1:1,1:size(X,1)] # cute the data to have same length
data20 = repeat(data_one,20)
data20 .= data20 .+ rand(StableRNG(1),size(data20)...)

#---
y = data_one
y = data20
@time  Unfold.solver_default(X,y;multithreading=false);
@time  Unfold.solver_default(X,y;multithreading=true);

@time  Unfold.solver_krylov(X,y);
@time  Unfold.solver_krylov(X,y;GPU = true);



#----
#using SuiteSparseGraphBLAS
#@time b7 = solver_gb(X,y);





# function solver_gb(
#     X,
#     data::AbstractArray{T,2};
#     stderror = false,
# ) where {T<:Union{Missing,<:Number}}
#     minfo = []
#     sizehint!(minfo, size(data,1))
#     beta = zeros(size(data,1),size(X,2)) # had issues with undef
    
#     X_l = GBMatrix(disallowmissing(X))
#     for ch = 1:size(data, 1)
# 	    # use the previous channel as a starting point
#         #ch == 1 || copyto!(view(beta, ch, :), view(beta, ch-1, :))
 
# 		beta[ch,:],h = lsmr!(@view(beta[ch, :]), X_l, @view(data[ch, ix]),log=true)

#         push!(minfo, h)
#     end

#     if stderror
#         stderror = calculate_stderror(X, data, beta) 
#         modelfit = Unfold.LinearModelFit(beta, ["lsmr",minfo], stderror)
#     else
#         modelfit = Unfold.LinearModelFit(beta, ["lsmr",minfo])
#     end
#     return modelfit
# end
