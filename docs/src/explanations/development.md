
## Install a dev-version of Unfold
In order to see and change the tutorials, you have to install a local dev-version of Unfold via:

`]dev --local Unfold` 

This clones the `git#main` into `./dev/Unfold`

### Instantiating the documentation environment
To generate documentation, we recommend to install LiveServer.jl - then you can do:
```julia
using LiveServer
servedocs(skip_dir=joinpath("docs","src","generated"),literate_dir=joinpath("docs","literate"))
```

- Next we have to make sure to be in the `Unfold/docs` folder, else the tutorial will not be able to find the data. Thus `cd("./docs")` in case you cd'ed already to the Unfold project. 
- And the `]activate .` to activate the docs-environment.
- Finally run `]instantiate` to install the required packages. Now you are ready to run the tutorials locally

