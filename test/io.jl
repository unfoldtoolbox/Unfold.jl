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
