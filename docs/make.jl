using Documenter
using Unfold
using AlgebraOfGraphics # can be removed with UnfoldMakie 0.3.0
using DocStringExtensions

using Literate
using Glob

GENERATED = joinpath(@__DIR__, "src", "generated")
SOURCE = joinpath(@__DIR__, "literate")
for subfolder ∈ ["explanations", "HowTo", "tutorials", "references"]
    local SOURCE_FILES = Glob.glob(subfolder * "/*.jl", SOURCE)
    foreach(fn -> Literate.markdown(fn, GENERATED * "/" * subfolder), SOURCE_FILES)

end



makedocs(
    sitename = "Unfold.jl Timeseries Analysis & Deconvolution",
    #root = joinpath(dirname(pathof(Unfold)), "..", "docs"),
    #prettyurls = get(ENV, "CI", nothing) == "true",
    repo = Documenter.Remotes.GitHub("unfoldtoolbox", "Unfold.jl"),
    pages = [
        "index.md",
        "Installing Julia + Unfold.jl" => "installation.md",
        "Tutorials" => [
            "rERP (mass univariate)" => "tutorials/lm_mu.md",
            "rERP (overlap correction)" => "tutorials/lm_overlap.md",
        ],
        "HowTo" => [
            "Multiple events" => "HowTo/multiple_events.md",
            "Change contrasts / coding schema" => "generated/HowTo/contrasts.md",
            "Standard errors" => "HowTo/standarderrors.md",
            "Marginal effects (must read re: splines)" => "generated/HowTo/effects.md",
            "Alternative Solvers (Robust, GPU, B2B)" => "HowTo/custom_solvers.md",
            #"Time domain basis functions"=>"generated/HowTo/timesplines.md",
            "Save and load Unfold models" => "generated/HowTo/unfold_io.md",
            "Duration-scaled basisfunctions (Hassall-style)" => "generated/HowTo/FIRduration.md",
            "🐍 Import EEG with PyMNE.jl" => "HowTo/pymne.md",
            "🐍 Calling Unfold.jl directly from Python" => "generated/HowTo/juliacall_unfold.md",
        ],
        "Explanations" => [
            "Non-Linear effects" => "./generated/explanations/nonlinear_effects.md",
            "Basisfunctions" => "./explanations/basisfunctions.md",
            "Predictions" => "./generated/explanations/predict.md",
            "Window Length Effect" => "./generated/explanations/window_length.md",
        ],
        "Reference" => [
            "Overview of package extensions" => "references/extensions.md",
            #"Development environment" => "explanations/development.md",
            "Solver/optimizer implementations" => "./generated/references/solver.md",
            "Solver benchmarks" => "./references/benchmarks.md",
            "API: Types" => "references/types.md",
            "API: Functions" => "references/functions.md",
        ],
        "Contributing" => ["90-contributing.md"],
        "Developer Guide" => ["91-developer.md"],
    ],
)

deploydocs(;
    repo = "github.com/unfoldtoolbox/Unfold.jl",
    push_preview = true,
    devbranch = "main",
)
