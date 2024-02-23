using UnfoldSim
using Unfold
using BSplineKit
using Krylov, CUDA
data, evts = UnfoldSim.predef_eeg()
data = repeat(data, 1, 30)'
data = data .+ rand(size(data)...)
bf = firbasis([-1, 1], 100)
evts.continuousB .= rand(size(evts, 1))

f = @formula(0 ~ 1 + condition + spl(continuous, 5) + spl(continuousB, 20))

# 1.3s vs 0.3s in matlab Oo
@time uf = designmatrix!(UnfoldLinearModelContinuousTime(Dict(Any => (f, bf))), evts);

@time fit(UnfoldModel, Dict(Any => (f, bf)), evts, data; solver = (x, y) -> return nothing);

@time fit(UnfoldModel, Dict(Any => (f, bf)), evts, data);


@time fit(
    UnfoldModel,
    Dict(Any => (f, bf)),
    evts,
    data;
    solver = (x, y) -> Unfold.solver_krylov(x, y; GPU = true),
);
