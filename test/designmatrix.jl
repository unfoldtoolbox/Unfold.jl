##
using Test,DataFrames,StatsModels
using unfold
tbl = DataFrame([1 4]',[:latency])
X = ones(size(tbl))
shouldBeNeg = zeros(4,4)
shouldBeNeg[1,:] = [1,0,0,1]
shouldBeNeg[2,:] = [0,1,0,0]
shouldBeNeg[3,:] = [0,0,1,0]
shouldBeNeg[4,:] = [0,0,0,1]

shouldBePos = zeros(4,4)
shouldBePos[1,:] = [0,0,0,0]
shouldBePos[2,:] = [1,0,0,0]
shouldBePos[3,:] = [0,1,0,0]
shouldBePos[4,:] = [0,0,1,0]

## test negative
basisfunction = unfold.firbasis(τ=(-3,0),sfreq = 1)
term =  unfold.TimeExpandedTerm(Term,basisfunction,:latency );
Xdc = unfold.time_expand(X,term,tbl)
@test all(isapprox.(Matrix(Xdc)[1:4,1:4], shouldBeNeg,atol=1e-15))

## Test Positive only
basisfunction = unfold.firbasis(τ=(1,4),sfreq = 1)
term =  unfold.TimeExpandedTerm(Term,basisfunction,:latency );
Xdc = unfold.time_expand(X,term,tbl)
println(Matrix(Xdc))

@test all(isapprox.(Matrix(Xdc)[1:4,1:4], shouldBePos,atol=1e-15))


## combining designmatrices
tbl = DataFrame([1 4]',[:latency])
X = ones(size(tbl))
basisfunction1 = unfold.firbasis(τ=(0,1),sfreq = 10,eventname="basis1")
basisfunction2 = unfold.firbasis(τ=(0,0.5),sfreq = 10,eventname="basis2")
f = @formula 0~1
Xdc1          = unfold.unfoldDesignmatrix(unfold.UnfoldLinearModel,f,tbl,basisfunction1)
Xdc2          = unfold.unfoldDesignmatrix(unfold.UnfoldLinearModel,f,tbl.+1,basisfunction2)

Xdc = Xdc1+Xdc2
@test size(Xdc.Xs,2) == size(Xdc1.Xs,2) + size(Xdc2.Xs,2)

# XXX concatenating of UnfoldLinearMixedModel designmatrices! Especially FeMat going to be more interesting....
df = unfold.unfoldFit(unfold.UnfoldLinearModel,Xdc,rand(size(Xdc.Xs,1)))
@test size(df.beta,1) == 17

ufresult = unfold.condense(df,tbl)

## Speedtest
if 1==0
    x = collect(range(10,stop=10000,step=10))
    tbl = DataFrame(reshape(x,1000,1),[:latency])
    X = ones(size(tbl,1),3).*[1,2,3]'


    basisfunction = unfold.firbasis(τ=(0,1),sfreq = 100,exact=true)
    term =  unfold.TimeExpandedTerm(Term,basisfunction,:latency );
    @time Xdc = Matrix(unfold.time_expand(X,term,tbl))
    @time Xdc = Matrix(unfold.time_expand2(X,term,tbl))
end
