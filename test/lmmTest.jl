#Pkg.activate("c:")
using Revise
using CSV
using Missings
using Plots
using StatsModels
##
import unfold
testCase = "testCaseMultisubject"
data = CSV.read("test\\"*testCase*"_data.csv", header=0)
data = convert(Matrix,data)
evts = CSV.read("test\\"*testCase*"_events.csv")

basis = unfold.firbasis(τ=(-0.5,1),sfreq=10)

form  = @formula y~1+condA*condB+(1|subject)




ufdesign          = unfold.DesignMatrix(form,evts,basis)
beta,history = unfold.fit(ufdesign,data)

plot(reshape(beta,length(basis.times),Int64(length(beta)/length(basis.times))))
