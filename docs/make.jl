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
                "Mass Univariate Mixed Model" =>"tutorials/lmm_mu.md",
                "LM Overlap correction" =>"tutorials/lm_overlap.md",
                "(t.b.d.) LMM + Overlap correction" =>"tutorials/lmm_overlap.md",
            ],
            "HowTo"=>["Overlap: Different events"multiple_events.md"],
            "Explanations"=>["Temporal Basisfunctions" => "./explanations/basisfunctions.md"],
            "Reference"=>["references/datastructures.md"],
            
        ])

deploydocs(; repo = "github.com/unfoldtoolbox/Unfold.jl", push_preview = true,        devbranch = "main")
