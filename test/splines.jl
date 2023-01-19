
include("test_utilities.jl")
include("setup.jl")
#include("../src/circsplinepredictors.jl")

data, evts = loadtestdata("test_case_3a") #
#
##
f_circspl = @formula 0 ~ 1 + conditionA + circspl(continuousA, 4, 0.1, 10) # 1
f_spl = @formula 0 ~ 1 + conditionA + spl(continuousA, 4) # 1
f = @formula 0 ~ 1 + conditionA + continuousA # 1
data_r = reshape(data, (1, :))
data_e, times = Unfold.epoch(data = data_r, tbl = evts, τ = (-1.0, 1.0), sfreq = 10) # cut the data into epochs

m_mul = coeftable(fit(UnfoldModel, f, evts, data_e, times))
# to experiment with circular splines. TODO: remove them once test cases got implemented
#test_circspl = fit(UnfoldModel, f_circspl, evts, data_e, times)
#test_spl = fit(UnfoldModel, f_spl, evts, data_e, times)
m_mul_spl = coeftable(fit(UnfoldModel, f_spl, evts, data_e, times))
m_mul_circspl = coeftable(fit(UnfoldModel, f_circspl, evts, data_e, times))
# asking for 4 splines should generate 4 splines 
@test length(unique(m_mul_spl.coefname)) == 6 # XXX check back with Unfold whether this is the same! could be n-1 splines in Unfold. We should keep that comparable I guess

# test safe prediction
m = fit(UnfoldModel, f_spl, evts, data_e, times)
r = predict(m,DataFrame(conditionA=[0,0],continuousA=[0.9,1.9]))
@test all(ismissing.(r.yhat[r.continuousA.==1.9]))
@test !any(ismissing.(r.yhat[r.continuousA.==0.9]))


# test time expanded
basisfunction = firbasis(τ = (-1, 1), sfreq = 10, name = "A")
m_tul = coeftable(fit(UnfoldModel, f, evts, data_r, basisfunction))
m_tul_spl = coeftable(fit(UnfoldModel, f_spl, evts, data_r, basisfunction))

# test safe predict
m = fit(UnfoldModel, f_spl, evts, data_r, basisfunction)
@test_broken predict(m,DataFrame(conditionA=[0,0,0],continuousA=[0.9,0.9,1.9]))
#@test_broken all(ismissing.)

#evts_grid = gridexpand() 
# results from timeexpanded and non should be equal
#yhat_tul  = predict(m_tul_spl,evts_grid)
#yhat_mul  = predict(m_mul_spl,evts_grid)
if 1 == 0
    using AlgebraOfGraphics
    yhat_mul.conditionA = categorical(yhat_mul.conditionA)
    yhat_mul.continuousA = categorical(yhat_mul.continuousA)
    m = mapping(:times, :yhat, color = :continuousA, linestyle = :conditionA)
    df = yhat_mul
    AlgebraOfGraphics.data(df) * visual(Lines) * m |> draw
end

# test much higher number of splines
f_spl_many = @formula 0 ~ 1 + spl(continuousA, 131) # 1
m_mul_spl_many = coeftable(fit(UnfoldModel, f_spl_many, evts, data_e, times))
@test length(unique(m_mul_spl_many.coefname)) == 132


## some circular test_spl
m = fit(UnfoldModel, f_circspl, evts, data_e, times)
f_evaluated = Unfold.formula(m)
m
effValues = [0.1,0.11,0.12,0.5,1,3,5,8,9.98,9.99,10]
effSingle = effects(Dict(:continuousA => effValues), m)
effMeaned = combine(groupby(effSingle, [:continuousA]), [:yhat] .=> mean .=> [:yhat])
@test length(unique(m_mul_circspl.coefname)) == 6
@test size(f_evaluated.rhs.terms[3].fun([1])) == (1,5)
@test size(modelcols(f_evaluated.rhs.terms[3],DataFrame(conditionA=[0],continuousA=[1]))) == (1,4)
# first and last value must approximately be the same
#'check whether there is something like basisname which makes first and last value be of different categories'
@test effMeaned.yhat[1] ≈ effMeaned.yhat[length(effMeaned.yhat)]
#'TODO: run these two tests to see if they work'
# check for first and second derivative whether they are similar in the border cases
# both formulae taken from https://mathformeremortals.wordpress.com/2013/01/12/a-numerical-second-derivative-from-three-points/
@testset "first and second derivative" begin    
    x0 = effMeaned.continuousA[length(effMeaned.continuousA)-2]
    x1 = effMeaned.continuousA[length(effMeaned.continuousA)-1]
    x2_1 = effMeaned.continuousA[length(effMeaned.continuousA)]
    x2_2 = effMeaned.continuousA[1]
    x3 = effMeaned.continuousA[2]
    x4 = effMeaned.continuousA[3]
    y0 = effMeaned.yhat[length(effMeaned.yhat)-2]
    y1 = effMeaned.yhat[length(effMeaned.yhat)-1]
    y2_1 = effMeaned.yhat[length(effMeaned.yhat)]
    y2_2 = effMeaned.yhat[1]
    y3 = effMeaned.yhat[2]
    y4 = effMeaned.yhat[3]
    @test (y2_1 - y1) / (x2_1 - x1) ≈ (y3 - y2_2) / (x3 - x2_2)
    @test hcat(2 / ((x1 - x0)*(x2_1 - x0)), -2 / ((x2_1 - x1)*(x1 - x0)), 2 / ((x2_1 - x1)*(x2_1 - x0))) * vcat(
        y0, y1, y2_1) ≈ hcat(2 / ((x3 - x2_2)*(x4 - x2_2)), -2 / ((x4 - x3)*(x3 - x2_2)), 2 / ((x4 - x3)*(x4 - x2_2))
         ) * vcat(y2_2, y3, y4)
end

@test mapValues([3,8,100,100,11], 0, 10) == [3,8,0,10,1]

# visual test
# inspired by https://unfoldtoolbox.github.io/Unfold.jl/dev/_literate/explanations/nonlinear_effects/
for i in collect(1:10)
    for nspl in [3,5,10]
        rng = MersenneTwister(i) # make repeatable
        n = 20 # datapoints
        evts = DataFrame(:x=>rand(rng,n))
        signal = (3*(evts.x .-0.5)).^2 .+ 0.5 .* rand(rng,n) .* 6
        # making a sinus-like curve by concatenating a parabola with itself but mirrored
        fullsignal = vcat(signal, signal .* -1 .+ 2 * maximum(signal))
        fullsignal = reshape(fullsignal,length(fullsignal),1,1)
        fullsignal = permutedims(fullsignal,[3,2,1])
        size(fullsignal)
        design = nothing
        if(nspl == 3)
            design = Dict(Any=>(@formula(0~1+circspl(x, 3, 0, 2)),[0]))
        elseif(nspl == 5)
            design = Dict(Any=>(@formula(0~1+circspl(x, 5, 0, 2)),[0]))
        elseif(nspl == 10)
            design = Dict(Any=>(@formula(0~1+circspl(x, 10, 0, 2)),[0]))
        else
            error("the visual test currently only works for 3, 5, and 10 splines. Received value: $(nspl)")
        end
        uf_spl = fit(UnfoldModel,design,DataFrame(:x=>vcat(evts.x, evts.x .+ 1)),fullsignal);

        p_spl = Unfold.effects(Dict(:x => range(0,stop=2,length=200)),uf_spl);

        pl = plot(DataFrame(:x=>vcat(evts.x, evts.x .+ 1)).x,fullsignal[1,1,:], xtickfontsize=30,ytickfontsize=30)
        lines!(p_spl.x,p_spl.yhat)
        display(pl)
        save("plots/visualTests_circspl/vistest$(i)_$(nspl)spl.svg", pl)
    end
end

