
## Installation

### getting into the right environment
First we have to install some packages. in julia you would do this either by putting a `]` in the REPL ("julia-commandline").

This should result in 

`(currentFolder) pkg>` 

(with `currentFolder` being the project you currently work in) - but if you see 

`(@v1.6) pkg>`

instead, you still have to activate your environment (using `cd("/path/to/your/project")` and `]activate .` or `]activate /path/to/your/project/`)

!!! note 
    You should have done this already to install Unfold in the first place. have a look at the Readme.md - there we use the Pkg.add("") syntax, which is equivalent to the `]` package manager.

### installing packages for the documentation
Now we are ready to add packages:

`add StatsModels,MixedModels,DataFrames,DSP.conv,Plots`

Next we have to make sure to be in the `Unfold/docs` folder, else the tutorial will not be able to find the data. Thus `cd("./docs")` in case you cd'ed already to the Unfold project.

