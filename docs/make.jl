using Documenter
using Unfold
using DocStringExtensions
using Plots
using Literate
using Glob

GENERATED = joinpath(@__DIR__, "src", "_literate")
SOURCE_FILES = Glob.glob("*/*.jl", GENERATED)
foreach(fn -> Literate.markdown(fn, GENERATED), SOURCE_FILES)


gr() # plots - can be removed after updating basisfunction
#unicodeplots()

makedocs(sitename="Unfold.jl Timeseries Analysis & Deconvolution",
        #root = joinpath(dirname(pathof(Unfold)), "..", "docs"),
        #prettyurls = get(ENV, "CI", nothing) == "true",
        pages = [
            "index.md",
            "Tutorials"=>[
                "Running these tutorials" => "tutorials/installation.md",
                "Mass Univariate" =>"tutorials/lm_mu.md",
                "LM Overlap correction" =>"tutorials/lm_overlap.md",
                "Mass Univariate Mixed Model" =>"tutorials/lmm_mu.md",
                "LMM + Overlap correction" =>"tutorials/lmm_overlap.md",
            ],
            "HowTo"=>[
		            "Overlap: Different events"=>"HowTo/multiple_events.md",
                    "Load Existing Dataset with PyMNE.jl"=>"HowTo/pymne.md" ,
                    "Custom Solvers / StandardErrors / B2B"=>"HowTo/custom_solvers.md",
                    "Unfold.jl directly from Python" => "_literate/pyjulia_unfold.md",
                    "P-values in Mass Univariate MixedModels" => "HowTo/lmm_pvalues.md",
		    "marginal effects (what to do with non-linear predictors)" =>"_literate/effects.md",
		    "Time Basis Functions"=>"_literate/timesplines.md",
                      ],
            "Explanations"=>[
		"Temporal Basisfunctions" => "./explanations/basisfunctions.md",
		"Non-Linear Effects" => "./_literate/nonlinear_effects.jl",],
            "Reference"=>["Types" => "references/types.md",
            "Functions" => "references/functions.md"],
            
        ])

deploydocs(; repo = "github.com/unfoldtoolbox/Unfold.jl", push_preview = true,        devbranch = "main")
