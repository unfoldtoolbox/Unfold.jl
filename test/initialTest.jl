using CSV
using Missings
using Revise
using StatsModels
import unfold
data = CSV.read("test\\testCase1_data.csv", header=0)
data = convert(Matrix,data)
evts = CSV.read("test\\testCase1_events.csv")

# we ultimately want to get there:
# designmatrix(:stimulus1,@formula ,evts)
basis = unfold.firbasis(τ=(-0.5,1),sfreq=10)
form = @formula y~1+conditionA*continuousA

evts_subset = evts[evts.type.=="stimulus2",:]

# Is it possible to pass Missing Type through? That would be fantastic
evts_subset.continuousA = Missings.disallowmissing(evts_subset.continuousA)
evts_subset.conditionA = Missings.disallowmissing(evts_subset.conditionA)
Xdc = unfold.DesignMatrix(form,evts_subset,basis)[5]

beta,history = unfold.fit(Xdc,data)


plot(reshape(beta,length(basis.times),Int64(length(beta)/length(basis.times))))
