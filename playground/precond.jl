using unfold
using IncompleteLU
#using IterativeSolver
using StatsModels

include("../test/test_utilities.jl")

data,evts = loadtestdata("testCase4","test/") #
f  = @formula 0~1+conditionA+conditionB # 1

# Timexpanded Univariate Linear
basisfunction = firbasis(τ=(-1,1),sfreq=100,name="A")
X = designmatrix(UnfoldLinearModel,f,evts,basisfunction)

##
#Crashes
X2 = IncompleteLU.ilu(X.Xs,τ=0.2)
#m,results = fit(UnfoldLinearModel,f,evts,data,basisfunction)

#X = m.X.Xs # extract designmatrix
X,data_r = unfold.zeropad(X,data')
@time beta,h = lsmr(X,data_r[1,:],log=true)

