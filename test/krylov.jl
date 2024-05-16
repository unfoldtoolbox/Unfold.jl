# This is not included in the testset because we cant use GPU on github action :(
@testset "krylov with missings" begin
    using Krylov, CUDA
    data, evts = UnfoldSim.predef_eeg()
    data = allowmissing(data)
    data[500:600] .= missing

    f = @formula 0 ~ 1 + condition + continuous
    # generate ModelStruct

    m = fit(
        UnfoldModel,
        [Any => (f, firbasis([-0.111, 0.2312], 100))],
        evts,
        data;
        solver = (x, y) -> Unfold.solver_krylov(x, y; GPU = false),
    )


end