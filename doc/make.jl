using Documenter
using unfold
using DocStringExtensions
makedocs(sitename="Unfold.jl",
        prettyurls = get(ENV, "CI", nothing) == "true",
        pages = [
            "index.md",
            "LM Tutorial" =>"lm_tutorial.md",
            "LMM Tutorial" =>"lmm_tutorial.md",
        ])
