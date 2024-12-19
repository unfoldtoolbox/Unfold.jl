
## Install a dev-version of Unfold

In order to see and change the tutorials, you have to install a local dev-version of Unfold via:

`]dev --local Unfold`

This clones the `git#main` into `./dev/Unfold`

### Instantiating the documentation environment

To generate documentation, we recommend to install LiveServer.jl - then you can do:

```julia
using LiveServer
servedocs(skip_dirs=joinpath("docs","src","generated"),literate_dir=joinpath("docs","literate"))
```

If you prefer a one-off:

- activate the  `./docs` folder (be sure to `]instantiate` the first time!)
- run `include("docs/make.jl")`
