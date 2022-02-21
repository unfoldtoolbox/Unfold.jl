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

makedocs(sitename="Unfold.jl",
        #root = joinpath(dirname(pathof(Unfold)), "..", "docs"),
        #prettyurls = get(ENV, "CI", nothing) == "true",
        sitename="Unfold.jl Timeseries Analysis & Deconvolution",
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
                    #"Unfold.jl directly from Python" => "HowTo/pyjulia_unfold.ipynb", # no idea if documenter allows that or not
                    "P-values in Mass Univariate MixedModels" => "HowTo/lmm_pvalues.md",
                      ],
            "Explanations"=>["Temporal Basisfunctions" => "./explanations/basisfunctions.md"],
            "Reference"=>["Types" => "references/types.md",
            "Functions" => "references/functions.md"],
            
        ])

deploydocs(; repo = "github.com/unfoldtoolbox/Unfold.jl", push_preview = true,        devbranch = "main")
