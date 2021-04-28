using Documenter
using unfold
using DocStringExtensions
makedocs(sitename="Unfold.jl",
        root = joinpath(dirname(pathof(unfold)), "..", "doc"),
        prettyurls = get(ENV, "CI", nothing) == "true",
        pages = [
            "index.md",
            "LM Tutorial" =>"lm_tutorial.md",
           # "LMM Tutorial" =>"lmm_tutorial.md",
        ])
