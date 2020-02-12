#Pkg.activate("c:")
using Revise
using CSV
using Missings
using Plots
using StatsModels
##
import unfold
testCase = "testCase1"
#testCase = "testCaseMultisubject"
data = CSV.read("test\\"*testCase*"_data.csv", header=0)
data = convert(Matrix,data)
evts = CSV.read("test\\"*testCase*"_events.csv")

# we ultimately want to get there:
# designmatrix(:stimulus1,@formula ,evts)

evts_subset = evts[evts.type.=="stimulus2",:]

# Is it possible to pass Missing Type through? That would be fantastic
evts_subset.continuousA = Missings.disallowmissing(evts_subset.continuousA)
evts_subset.conditionA  = Missings.disallowmissing(evts_subset.conditionA)


basis = unfold.firbasis(τ=(-0.5,1),sfreq=10)


form  = @formula y~1+conditionA*continuousA+(1|subject)

ufdesign          = unfold.DesignMatrix(form,evts_subset,basis)
beta,history = unfold.fit(ufdesign.Xdc,data)

plot(reshape(beta,length(basis.times),Int64(length(beta)/length(basis.times))))
