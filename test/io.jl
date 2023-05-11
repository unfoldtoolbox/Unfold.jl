using UnfoldSim

save_path = tempdir()

data1, evts1 = UnfoldSim.predef_eeg(; n_repeats=100, noiselevel=0.8);
evts1[!, :type] = repeat(["event_A", "event_B"], nrow(evts1) ÷ 2);

bf1_A = firbasis(τ=[-0.1, 1], sfreq=100, name="event_A");
bf1_B = firbasis(τ=[-0.1, 1], sfreq=100, name="event_B");

f1_A = @formula 0 ~ 1;
f1_B = @formula 0 ~ 1 + condition + spl(continuous, 3);

bfDict1 = Dict("event_A" => (f1_A, bf1_A), "event_B" => (f1_B, bf1_B));

m1 = Unfold.fit(UnfoldModel, bfDict1, evts1, data1, eventcolumn="type");

@testset "SingleSubjectDesign with two event types and splines" begin
    # save the model to a compressed .jld2 file and load it again
    save(joinpath(save_path, "m1_compressed.jld2"), m1; compress=true)
    m1_loaded = load(joinpath(save_path, "m1_compressed.jld2"), UnfoldModel, generate_Xs=true)

    @test m1.modelfit.estimate == m1_loaded.modelfit.estimate
    @test m1.designmatrix.events == m1_loaded.designmatrix.events
    @test m1.designmatrix.Xs == m1_loaded.designmatrix.Xs

    eff1 = effects(Dict(:condition=>["car","face"],:continuous=>-5:1),m1)
    eff1_loaded = effects(Dict(:condition=>["car","face"],:continuous=>-5:1),m1_loaded)
    @test eff1 == eff1_loaded
end

#----
data2, evts2 = UnfoldSim.predef_2x2(; n_subjects=5, return_epoched=true);

f2 = @formula(0 ~ 1 + A * B + (A * B | subject));
τ2 = [-0.1, 1];
sfreq2 = 100;
times2 = range(τ2[1], length=size(data2, 1), step=1 ./ sfreq2);

m2 = Unfold.fit(UnfoldModel, Dict(Any => (f2, times2)), evts2, reshape(data2, 1, size(data2)...));

save(joinpath(save_path, "m2_compressed.jld2"), m2; compress=true)
m2_loaded = load(joinpath(save_path, "m2_compressed.jld2"), UnfoldModel, generate_Xs=true)

@testset "2x2 MultiSubjectDesign Mixed-effects model" begin
    # save the model to a compressed .jld2 file and load it again
    save(joinpath(save_path, "m2_compressed.jld2"), m2; compress=true)
    m2_loaded = load(joinpath(save_path, "m2_compressed.jld2"), UnfoldModel, generate_Xs=true)

    @test m2.modelfit.fits == m2_loaded.modelfit.fits
    @test m2.designmatrix.events == m2_loaded.designmatrix.events

    @test_broken m2.designmatrix.Xs == m2_loaded.designmatrix.Xs
end