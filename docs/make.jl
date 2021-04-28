using Documenter
using unfold
using DocStringExtensions
makedocs(sitename="Unfold.jl",
        root = joinpath(dirname(pathof(unfold)), "..", "docs"),
        prettyurls = get(ENV, "CI", nothing) == "true",
        pages = [
            "index.md",
            "LM Tutorial" =>"lm_tutorial.md",
           # "LMM Tutorial" =>"lmm_tutorial.md",
        ])

deploydocs(; repo = "github.com/unfoldtoolbox/unfold.jl", push_preview = true)