using Revise
using CSV
using StatsModels
using MixedModels
using DataFrames
using DataFramesMeta
import Plots
import unfold
using Gadfly
using DataFramesMeta


testCase = "testCaseMultisubject2"
data = CSV.read("test\\"*testCase*"_data.csv", header=0)
data = convert(Matrix,data)
evts = CSV.read("test\\"*testCase*"_events.csv")
categorical!(evts,:subject);


f  = @formula 0~1+condA*condB+(1+condA*condB|subject);
f  = @formula 0~1+condA*condB+(1|subject);
basisfunction = unfold.firbasis(τ=(-0.05,.375),sfreq=40)
@time mm = unfold.fit(unfold.UnfoldLinearMixedModel,f,evts,data,basisfunction)
