using Unfold
using UnfoldSim
using CairoMakie
using DataFrames
using UnfoldMakie

data, evts = UnfoldSim.predef_eeg(; n_repeats = 5, noiselevel = 0.8)
uf = fit(
    UnfoldModel,
    [
        "car" => (@formula(0 ~ 1 + spl(continuous, 4)), firbasis((-0.1, 1), 100)),
        "face" => (@formula(0 ~ 1 + continuous^2), firbasis((-0.2, 0.4), 100)),
    ],
    evts,
    data;
    eventcolumn = "condition",
    show_progress = false,
)
#----
ix = 500
lines(predict(uf)[1, 1:ix])

# predict(m,overlap=false)
#---
f = Figure()
heatmap(
    f[1, 1],
    predict(uf; keep_basis = ["car"], epoch_to = "car", eventcolumn = :condition)[1, :, :],
)
heatmap(
    f[1, 2],
    predict(uf; keep_basis = ["car"], epoch_to = "face", eventcolumn = :condition)[1, :, :],
)
heatmap(
    f[2, 2],
    predict(uf; keep_basis = ["face"], epoch_to = "face", eventcolumn = :condition)[
        1,
        :,
        :,
    ],
)
heatmap(
    f[2, 1],
    predict(uf; keep_basis = ["face"], epoch_to = "car", eventcolumn = :condition)[1, :, :],
)
f
#---

predict(
    uf,
    DataFrame(:continuous => [1, 5], :condition => ["car", "face"], :latency => [1, 40]),
)[1][
    1,
    :,
    :,
]' |> series


#---
eff = effects(Dict(:continuous => [0, 1, 2]), uf)

plot_erp(eff; mapping = (; color = :continuous, col = :eventname))

