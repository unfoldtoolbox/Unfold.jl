##
using Test, StatsModels
using DataFrames
import Logging
import unfold
include("test_utilities.jl")

Logging.global_logger(Logging.SimpleLogger(stdout, Logging.Debug))


data,evts = loadtestdata("testCase3") #
append!(data,zeros(1000))
data = reshape(data,1,:)
data = vcat(data,data)
data = data.+ 1*randn(size(data)) # we have to add minimal noise, else mixed models crashes.
data_missing = Array{Union{Missing,Number}}(undef,size(data))
data_missing .= data
data_missing[4500:4600] .= missing

categorical!(evts,:subject)
f  = @formula 0~1+condA+condB + (1+condA+condB|subject)
#f  = @formula 0~1 + (1|subject)



# cut the data into epochs
# TODO This ignores subject bounds
data_e,times = unfold.epoch(data=data,tbl=evts,τ=(-1.,1.9),sfreq=10)
evts_e,data_e = unfold.dropMissingEpochs(evts,data_e)



basisfunction = unfold.firbasis(τ=(-0.1,.3),sfreq=10)

f  = @formula 0~1+condA+condB # 1
Xs = unfold.designmatrix(unfold.UnfoldLinearModel,f,evts_e)
ufModel_A = unfold.fit!(unfold.UnfoldLinearModel,Xs,data_e)

basisfunction = unfold.firbasis(τ=(-0.1,.3),sfreq=10)
Xs = unfold.designmatrix(unfold.UnfoldLinearModel,f,evts,basisfunction)
basisfunction = unfold.firbasis(τ=(-0.1,.5),sfreq=10)
Xs2 = unfold.designmatrix(unfold.UnfoldLinearModel,f,evts,basisfunction)
Xs = Xs+Xs2

ufModel_B = unfold.fit!(unfold.UnfoldLinearModel,Xs,data)


f  = @formula 0~1+condA+condB + (1+condA|subject)
Xs = unfold.designmatrix(unfold.UnfoldLinearMixedModel,f,evts_e)

ufModel_C = unfold.fit!(unfold.UnfoldLinearMixedModel,Xs,data_e)

basisfunction = unfold.firbasis(τ=(-0.1,.3),sfreq=10)
Xs = unfold.designmatrix(unfold.UnfoldLinearMixedModel,f,evts,basisfunction)

ufModel_D = unfold.fit!(unfold.UnfoldLinearMixedModel,Xs,data)


ufA = unfold.condense_long(ufModel_A,times)
ufB = unfold.condense_long(ufModel_B)
ufC = unfold.condense_long(ufModel_C,times)
ufD = unfold.condense_long(ufModel_D)

plot(ufA.colnames_basis,ufA.estimate)
plot(ufC.colnames_basis,ufC.estimate)
plot(ufB.colnames_basis,ufB.estimate)
plot(ufD.colnames_basis,ufD.estimate)
