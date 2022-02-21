
## [Installation]((@id install_instruct))

### getting into the right environment
First we have to install some packages. in julia you would do this either by putting a `]` in the REPL ("julia-commandline").

This should result in `(currentFolder) pkg>` (with `currentFolder` being the project you currently work in) - but if you see 

`(@v1.7) pkg>` instead, you still have to activate your environment (using `cd("/path/to/your/project")` and `]activate .` or `]activate /path/to/your/project/`)

!!! note 
    You should have done this already to install Unfold in the first place. have a look at the Readme.md - there we use the Pkg.add("") syntax, which is equivalent to the `]` package manager.

### Install a dev-version of Unfold
In order to see and change the tutorials, you have to install a local dev-version of Unfold via:
`]dev --local Unfold` - which installs it in `./dev/Unfold`

### Instantiating the documentation environment
- Now we have to add the packages required for the documentation.
- Next we have to make sure to be in the `Unfold/docs` folder, else the tutorial will not be able to find the data. Thus `cd("./docs")` in case you cd'ed already to the Unfold project. 
- And the `]activate .` to activate the docs-environment.
- Finally run `]instantiate` to install the required packages. Now you are ready to run the tutorials locally



