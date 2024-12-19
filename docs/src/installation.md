
# [Installation](@id install_instruct)

## Installing Julia

The easiest way to install julia is using [`juliaup`](https://github.com/JuliaLang/juliaup)

TLDR;

- Windows: `winget install julia -s msstore`
- Mac/Linux: `curl -fsSL https://install.julialang.org | sh`

We further recommend to use VSCode. Make sure to install the Julia-Plugin, and install Revise.jl - [a tutorial with screenshots can be found here](http://www.simtech-summerschool.de/installation/julia.html)

## Installing Unfold.jl

You can enter the package manager (similar to conda) using `]` in the REPL ("julia-commandline").

This should result in `(currentFolder) pkg>` (with `currentFolder` being the project you currently work in)

!!! hint
    if you see `(@v1.9) pkg>` instead, you still have to activate your environment. This can be done using:

    `cd("/path/to/your/project")`
     and `]activate .`

     or alternatively `]activate /path/to/your/project/`

Now you can do
`pkg> add Unfold`

and after some installation:

`julia> using Unfold` in the REPL
