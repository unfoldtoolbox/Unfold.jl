#Pkg.activate("c:")
using Revise
using CSV

using StatsModels
using MixedModels
using DataFrames
##
import unfold
testCase = "testCaseMultisubject_small"
data = CSV.read("test\\"*testCase*"_data.csv", header=0)
data = convert(Matrix,data)
evts = CSV.read("test\\"*testCase*"_events.csv")
categorical!(evts,:subject)
#evts = evts[evts.subject.>=20,:]

basisfunction = unfold.firbasis(τ=(-0.1,.5),sfreq=10)

##
f  = @formula 0~1+condA*condB+(1+condA*condB|subject)
f = @formula 0 ~ 1 + condA + (1 |subject)
#formula  = apply_schema(f,schema(f,evts),LinearMixedModel)#+(1|subject)

#form = FormulaTerm(formula.lhs, unfold.TimeExpandedTerm.(formula.rhs,Ref(basisfunction)))
#_,Xs = modelcols(form, evts)


@time mm = unfold.fit(unfold.UnfoldLinearMixedModel,f,evts,data,basisfunction)

using Gadfly
Gadfly.push_theme(:dark)
using DataFramesMeta
d = @linq mm.results |> where(:group.=="fixed")
plot(d,x=:time,y=:estimate,color=:term,Geom.LineGeometry)
d = @linq mm.results |> where(:group.=="subject")
plot(d,x=:time,y=:estimate,color=:term,Geom.LineGeometry)
