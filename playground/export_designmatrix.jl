using Test, StatsModels
using DataFrames

using unfold
# From folder "test" in unfold
include("test_utilities.jl")

tau = (-.5,1.5)

data,evts = loadtestdata("testCase5") #
data = reshape(data,(1,:))
data_e,times = unfold.epoch(data=data,tbl=evts,τ=tau,sfreq=1000) # cut the data into epochs

f  = @formula 0~1#+conditionA#+conditionB # 4
bf1 = firbasis(τ=tau,sfreq=1000,name="ev1")


fi,r1 =  fit(UnfoldLinearModel,f,evts,data,bf1)
fi2,r2 =  fit(UnfoldLinearModel,f,evts,data_e,times)

if 1==0
    # plot it if needed
    # commented because Makie takes long to load
    using Makie
    lines(r1.colname_basis,r1.estimate)
    lines!(r2.colname_basis,r2.estimate,  color=:red)
end

#---
using SparseArrays
I, J, V = findnz(fi.X.Xs)
df = DataFrame([Int,Int,Int], [:I, :J,:V])
for k = 1:length(I)
    if V[k] == 1
        push!(df,(I[k],J[k],V[k]))
    end
end

CSV.write("Xs_large_y~1.csv", df)

#CSV.write("Xs_large_y~1.csv",DataFrame(fi.X.Xs))

# also export a slightly more complex designmatrix
f  = @formula 0~1+conditionA#+conditionB # 4
fi3,r3 =  fit(UnfoldLinearModel,f,evts,data,bf1)
CSV.write("Xs_y~1+conditionA.csv",DataFrame(fi3.X.Xs))

