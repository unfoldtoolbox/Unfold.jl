using Documenter
using Unfold
using DocStringExtensions
using Plots
gr()
#unicodeplots()

makedocs(sitename="Unfold.jl",
        root = joinpath(dirname(pathof(Unfold)), "..", "docs"),
        prettyurls = get(ENV, "CI", nothing) == "true",
        pages = [
            "index.md",
            "Tutorials",[
                "Mass Univariate Tutorial" =>"tutorials/lm_mu.md",
                "Overlap correction Tutorial" =>"tutorials/lm_overlap.md",
                "(to be overhauled) LMM Tutorial" =>"tutorials/lmm_tutorial.md",],
            "HowTo",[],
            "Explanations",["Temporal Basisfunctions" => "explanations/basisfunctions.md"],
            "Reference",[],
            
        ])

deploydocs(; repo = "github.com/unfoldtoolbox/Unfold.jl", push_preview = true,        devbranch = "main")
