using Documenter
using Unfold
using AlgebraOfGraphics # can be removed with UnfoldMakie 0.3.0
using DocStringExtensions

using Literate
using Glob

GENERATED = joinpath(@__DIR__, "src", "generated")
SOURCE = joinpath(@__DIR__,"literate")
for subfolder ∈ ["explanations","HowTo","tutorials"]
    local SOURCE_FILES = Glob.glob(subfolder*"/*.jl", SOURCE)
    foreach(fn -> Literate.markdown(fn, GENERATED*"/"*subfolder), SOURCE_FILES)

end



makedocs(sitename="Unfold.jl Timeseries Analysis & Deconvolution",
        #root = joinpath(dirname(pathof(Unfold)), "..", "docs"),
        #prettyurls = get(ENV, "CI", nothing) == "true",
        repo = Documenter.Remotes.GitHub("unfoldtoolbox", "Unfold.jl"),
        pages = [
            "index.md",
            "Installing Julia + Unfold.jl" => "installation.md",
            "Tutorials"=>[

                "Mass Univariate" =>"tutorials/lm_mu.md",
                "LM Overlap correction" =>"tutorials/lm_overlap.md",
                "Mass Univariate Mixed Model" =>"tutorials/lmm_mu.md",
                "LMM + Overlap correction" =>"tutorials/lmm_overlap.md",
            ],
            "HowTo"=>[
		            "Overlap: Different events"=>"HowTo/multiple_events.md",
                    "Load Existing Dataset with PyMNE.jl"=>"HowTo/pymne.md" ,
                    "Custom Solvers / StandardErrors / B2B"=>"HowTo/custom_solvers.md",
                    "🐍 Calling Unfold.jl directly from Python" => "generated/HowTo/juliacall_unfold.md",
                    "P-values in Mass Univariate MixedModels" => "HowTo/lmm_pvalues.md",
		            "Marginal effects (what to do with non-linear predictors)" =>"generated/HowTo/effects.md",
		            "Time Basis Functions"=>"generated/HowTo/timesplines.md",
                      ],
            "Explanations"=>[
		            "Temporal Basisfunctions" => "./explanations/basisfunctions.md",
		            "Non-Linear Effects" => "./generated/explanations/nonlinear_effects.md",],
            "Reference"=>[
                    "Overview of Package Extensions" => "references/extensions.md",
                    "Development Environment" => "explanations/development.md",
                    "Types" => "references/types.md",
                    "Functions" => "references/functions.md"],
            
        ])

deploydocs(; repo = "github.com/unfoldtoolbox/Unfold.jl", push_preview = true,        devbranch = "main")
