using Documenter
using Unfold
using AlgebraOfGraphics # can be removed with UnfoldMakie 0.3.0
using DocStringExtensions

using Literate
using Glob

GENERATED = joinpath(@__DIR__, "src", "generated")
SOURCE = joinpath(@__DIR__, "literate")
for subfolder ∈ ["explanations", "HowTo", "tutorials"]
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
            "Mass univariate LM" => "tutorials/lm_mu.md",
            "LM overlap correction" => "tutorials/lm_overlap.md",
            "Mass univariate Mixed Model" => "tutorials/lmm_mu.md",
            "LMM + overlap correction" => "tutorials/lmm_overlap.md",
        ],
        "HowTo" => [
            "Overlap: Multiple events" => "HowTo/multiple_events.md",
            "Import EEG with 🐍 PyMNE.jl" => "HowTo/pymne.md",
            "Standard errors" => "HowTo/standarderrors.md",
            "Alternative Solvers (Robust, GPU, B2B)" => "HowTo/custom_solvers.md",
            "🐍 Calling Unfold.jl directly from Python" => "generated/HowTo/juliacall_unfold.md",
            "P-values for mixedModels" => "HowTo/lmm_pvalues.md",
            "Marginal effects (focus on non-linear predictors)" => "generated/HowTo/effects.md",
            #"Time domain basis functions"=>"generated/HowTo/timesplines.md",
            "Save and load Unfold models" => "generated/HowTo/unfold_io.md",
            "Change contrasts / coding schema" => "generated/HowTo/contrasts.md",
        ],
        "Explanations" => [
            "About basisfunctions" => "./explanations/basisfunctions.md",
            "Non-Linear effects" => "./generated/explanations/nonlinear_effects.md",
            "Window Length Effect" => "./generated/explanations/window_length.md",
        ],
        "Reference" => [
            "Overview of package extensions" => "references/extensions.md",
            "Development environment" => "explanations/development.md",
            "API: Types" => "references/types.md",
            "API: Functions" => "references/functions.md",
        ],
    ],
)

deploydocs(;
    repo = "github.com/unfoldtoolbox/Unfold.jl",
    push_preview = true,
    devbranch = "main",
)
