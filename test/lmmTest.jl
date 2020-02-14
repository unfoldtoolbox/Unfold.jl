#Pkg.activate("c:")
using Revise
using CSV
using Missings
using Plots
using StatsModels
using MixedModels
using DataFrames
##
import unfold
testCase = "testCaseMultisubject"
data = CSV.read("test\\"*testCase*"_data.csv", header=0)
data = convert(Matrix,data)
evts = CSV.read("test\\"*testCase*"_events.csv")
categorical!(evts,:subject)
#evts = evts[evts.subject.>=20,:]

basisfunction = unfold.firbasis(τ=(-0.5,1),sfreq=10)

##
f  = @formula 0~1+condA+(1|subject)
formula  = apply_schema(f,schema(f,evts),LinearMixedModel)#+(1|subject)

form = FormulaTerm(formula.lhs, unfold.TimeExpandedTerm.(formula.rhs,Ref(basisfunction)))
_,Xs = modelcols(form, evts)
#formula2 = FormulaTerm(formula.lhs, unfold.TimeExpandedTerm.(formula.rhs,Ref(basisfunction)))
re = Xs[2]
fe1 = MixedModels.FeMat(Xs[1],[""])
fe2 = MixedModels.FeMat( SparseMatrixCSC(reshape(float(sparse(data[1:size(Xs[1],1)])), (:, 1))),[""])
#X = modelcols(formula2,evts)

unfold.fit(unfold.UnfoldLinearMixedModel,f,evts,sparse(data),basisfunction)

plot(reshape(beta,length(basis.times),Int64(length(beta)/length(basis.times))))



###
t = formula.term[2]
d =
lhs = t.lhs
z = Matrix(transpose(modelcols(lhs, d)))
cnames = coefnames(lhs)
T = eltype(z)
S = size(z, 1)
grp = t.rhs
m = reshape(1:abs2(S), (S, S))
inds = sizehint!(Int[], (S * (S + 1)) >> 1)
for j = 1:S, i = j:S
    push!(inds, m[i, j])
end
refs, levels = _ranef_refs(grp, d)

ReMat{T,S}(
    grp,
    refs,
    levels,
    isa(cnames, String) ? [cnames] : collect(cnames),
    z,
    z,
    LowerTriangular(Matrix{T}(I, S, S)),
    inds,
    adjA(refs, z),
    Matrix{T}(undef, (S, length(levels))),
)

##
