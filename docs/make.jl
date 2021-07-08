using Documenter
using Unfold
using DocStringExtensions
using Plots
gr()
#unicodeplots()

makedocs(sitename="Unfold.jl",
        #root = joinpath(dirname(pathof(Unfold)), "..", "docs"),
        #prettyurls = get(ENV, "CI", nothing) == "true",
        pages = [
            "index.md",
            "Tutorials"=>[
                "Running these tutorials" => "tutorials/installation.md",
                "Mass Univariate" =>"tutorials/lm_mu.md",
                "Overlap correction" =>"tutorials/lm_overlap.md",
                "(to be overhauled) Mass Univariate Mixed Model" =>"tutorials/lmm_mu.md",
                "(to be overhauled) Mixed Model with Overlap" =>"tutorials/lmm_overlap.md",
            ],
            "HowTo"=>[],
            "Explanations"=>["Temporal Basisfunctions" => "./explanations/basisfunctions.md"],
            "Reference"=>[],
            
        ])

deploydocs(; repo = "github.com/unfoldtoolbox/Unfold.jl", push_preview = true,        devbranch = "main")
