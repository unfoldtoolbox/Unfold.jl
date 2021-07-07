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
            "Tutorials" => [ "LM Tutorial" =>"lm_tutorial.md",
                             "LMM Tutorial" =>"lmm_tutorial.md"],
            "How-To" =>[],
            "Explanations" =>[],
            "Reference" => [],

        ])

deploydocs(; repo = "github.com/unfoldtoolbox/Unfold.jl", push_preview = true,        devbranch = "main")
