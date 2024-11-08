using UnfoldSim
using DataFrames
save_path = mktempdir(; cleanup = false)#tempdir()

## 1. Test data set:
# - Generate a P1/N1/P3 complex for one subject (using UnfoldSim)
# - Use a Continuous Time Unfold model with two event types
# - Use splines to model the continuous predictor variable



data1, evts1 = UnfoldSim.predef_eeg(; n_repeats = 100, noiselevel = 0.8);
# Assume that the data is from two different events
evts1[!, :type] = repeat(["event_A", "event_B"], nrow(evts1) ÷ 2);

bf1_A = firbasis(τ = [-0.1, 1], sfreq = 100, name = "event_A");
bf1_B = firbasis(τ = [-0.1, 1], sfreq = 100, name = "event_B");

f1_A = @formula 0 ~ 1;
f1_B = @formula 0 ~ 1 + condition + spl(continuous, 4);

bfDict1 = ["event_A" => (f1_A, bf1_A), "event_B" => (f1_B, bf1_B)];

data1_e, times = Unfold.epoch(data1, evts1, [-0.1, 1], 100)
bfDict1_e = ["event_A" => (f1_A, times), "event_B" => (f1_B, times)];


for deconv in [false, true]

    if deconv

        m1 = Unfold.fit(UnfoldModel, bfDict1, evts1, data1, eventcolumn = "type")
    else
        m1 = Unfold.fit(UnfoldLinearModel, bfDict1_e, evts1, data1_e; eventcolumn = "type")
    end
    @testset "SingleSubjectDesign with two event types and splines" begin
        # save the model to a compressed .jld2 file and load it again

        save(joinpath(save_path, "m1_compressed2.jld2"), m1; compress = true)
        m1_loaded = load(
            joinpath(save_path, "m1_compressed2.jld2"),
            UnfoldModel,
            generate_Xs = true,
        )

        @test isempty(Unfold.modelmatrices(designmatrix(m1_loaded))[1]) == false
        @test typeof(m1) == typeof(m1_loaded)
        @test coeftable(m1) == coeftable(m1_loaded)
        @test m1.modelfit.estimate == m1_loaded.modelfit.estimate
        @test m1.designmatrix[end].events == m1_loaded.designmatrix[end].events

        # In the loaded version one gets two matrices instead of one.
        # The dimensions do also not match.
        # Probably the designmatrix reconstruction in the load function needs to be changed.
        @test m1.designmatrix[end].modelmatrix == m1_loaded.designmatrix[end].modelmatrix

        # Test whether the effects function works with the loaded model
        # and the results match the ones of the original model
        eff1 = effects(Dict(:condition => ["car", "face"], :continuous => -5:1), m1)
        eff1_loaded =
            effects(Dict(:condition => ["car", "face"], :continuous => -5:1), m1_loaded)
        @test eff1 == eff1_loaded

        # load the model without reconstructing the designmatrix
        m1_loaded_without_dm = load(
            joinpath(save_path, "m1_compressed2.jld2"),
            UnfoldModel,
            generate_Xs = false,
        )


        # ismissing should only be true fr the deconv case
        @test isempty(modelmatrix(designmatrix(m1_loaded_without_dm)[2])) ==
              (deconv == true)

    end
end
#----
## 2. Test data set:
# - Generate a 2x2 design with Hanning window for multiple subjects (using UnfoldSim)
# - Use a Mixed-effects Unfold model

data2, evts2 = UnfoldSim.predef_2x2(; n_subjects = 5, return_epoched = true);
data2 = reshape(data2, size(data2, 1), :)


# Define a model formula with interaction term and random effects for subjects
f2 = @formula(0 ~ 1 + A * B + (1 | subject));
τ2 = [-0.1, 1];
sfreq2 = 100;
times2 = range(τ2[1], length = size(data2, 1), step = 1 ./ sfreq2);

m2 = Unfold.fit(
    UnfoldModel,
    Dict(Any => (f2, times2)),
    evts2,
    reshape(data2, 1, size(data2)...),
);

save(joinpath(save_path, "m2_compressed2.jld2"), m2; compress = true)
m2_loaded =
    load(joinpath(save_path, "m2_compressed2.jld2"), UnfoldModel, generate_Xs = true)

@testset "2x2 MultiSubjectDesign Mixed-effects model" begin
    # save the model to a compressed .jld2 file and load it again
    save(joinpath(save_path, "m2_compressed2.jld2"), m2; compress = true)
    m2_loaded =
        load(joinpath(save_path, "m2_compressed2.jld2"), UnfoldModel, generate_Xs = true)


    @test isempty(Unfold.modelmatrices(designmatrix(m2_loaded))[1]) == false

    @test typeof(m2) == typeof(m2_loaded)
    @test coeftable(m2) == coeftable(m2_loaded)
    @test modelfit(m2).fits == modelfit(m2_loaded).fits
    @test Unfold.events(m2) == Unfold.events(m2_loaded)
    @test modelmatrix(m2) == modelmatrix(m2_loaded)

    # Test whether the effects function works with the loaded models
    # and the results match the ones of the original model
    conditions = Dict(:A => levels(evts2.A), :B => levels(evts2.B))

    # The effects function is currently not defined for UnfoldLinearMixedModel
    #eff2 = effects(conditions, m2)
    #eff2_loaded = effects(conditions, m2_loaded)

    @test_broken eff2 == eff2_loaded

    # load the model without reconstructing the designmatrix
    m2_loaded_without_dm =
        load(joinpath(save_path, "m2_compressed2.jld2"), UnfoldModel, generate_Xs = false)

    @test isempty(modelmatrix(designmatrix(m2_loaded_without_dm))[1]) == true
end
