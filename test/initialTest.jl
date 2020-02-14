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

##
basis = unfold.firbasis(τ=(-0.5,1),sfreq=10)
form  = @formula y~1+conditionA*continuousA#+(1|subject)

ufdesign          = unfold.DesignMatrix(form,evts_subset,basis)
beta,history = unfold.fit(ufdesign,data)

##
f = @formula 0~1+conditionA*continuousA
tbl = evts
formula  = apply_schema(f,schema(f,tbl))#+(1|subject)
formula = unfold.TimeExpandedTerm(formula.rhs,basisfunction)
Xdc = modelcols(form2,evts_subset)


basisfunction = unfold.firbasis(τ=(-0.5,1),sfreq=10)

m = unfold.fit(unfold.UnfoldLinearModel,f,tbl,dropdims(data,dims=1),basisfunction)


plot(reshape(m.beta,length(m.basisfunction.times),Int64(length(m.beta)/length(m.basisfunction.times))))
