# Benchmarks

We ran benchmarks on 2024-11-07 as described in `./benchmark/cuda/solver_comparison.jl`. Given that some were run on a GPU, we cannot run them on continuous-integration online.

!!! important
    - Allocations are only CPU allocations - GPU allocations were not counted.
    - Solvers other than `default_multi` are currently *NOT* multi-threaded
    - Solvers other than `default_multi` and `krylov_gpu` solve $X'Xb = X'y$ instead of $Xb=y$ directly. They are likely less accurate, but should be faster for multi-channel data, as we can precalulate cholesky, qr or similar & the to-be-inverted matrix is much smaller.

## Small Model

```julia
n_channels = 1,
sfreq = 10,
n_splines = 4,
n_repeats = 10;
```

| **gpu** | **method**     | **el\_type** | **time**    | **GB**      | **percent\_X\_filled** | **sizeDesign** | **n\_channels** | **overlap** | **comment**         |
|--------:|---------------:|-------------:|------------:|------------:|-----------------------:|---------------:|----------------:|------------:|--------------------:|
| true    | cholesky       | Float64      |             |             | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  | PosDefException(-1) |
| false   | cholesky       | Float64      | **0.00056** | 0.00069     | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  |                     |
| true    | intern         | Float64      | 0.00088     | 0.00017     | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  |                     |
| false   | intern         | Float64      | 0.0011      | 0.00069     | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  |                     |
| true    | qr             | Float64      | 0.0013      | 0.00019     | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  |                     |
| false   | cg             | Float64      | 0.0015      | 0.00057     | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  |                     |
| false   | default\_multi | Float64      | 0.0017      | **0.00016** | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  |                     |
| false   | qr             | Float64      | 0.002       | 0.00076     | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  |                     |
| true    | cg             | Float64      | 0.0054      | 0.00056     | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  |                     |
| true    | pinv           | Float64      | 0.0054      | 0.00032     | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  |                     |
| false   | pinv           | Float64      | 0.016       | 0.0016      | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  |                     |
| true    | krylov\_gpu    | Float64      | 0.032       | 0.0013      | 0.068                  | (1190, 130)    | 1               | (0.2, 0.2)  |                     |

## small-to-midsize: multi-channel

```julia
n_channels = 128,
sfreq = 100,
n_splines = 4,
n_repeats = 200;
```

### Float64

| **gpu** | **method**     | **el\_type** | **time** | **GB**   | **percent\_X\_filled** | **sizeDesign** | **n\_channels** | **overlap** | **comment**         |
|--------:|---------------:|-------------:|---------:|---------:|-----------------------:|---------------:|----------------:|------------:|--------------------:|
| true    | cholesky       | Float64      |          |          | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  | PosDefException(-1) |
| true    | qr             | Float64      | **0.38** | 0.25     | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| true    | pinv           | Float64      | 0.42     | 0.26     | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| true    | intern         | Float64      | 0.7      | **0.25** | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| true    | cg             | Float64      | 1.2      | 0.32     | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | cholesky       | Float64      | 1.5      | 0.31     | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | qr             | Float64      | 1.7      | 0.31     | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | cg             | Float64      | 2.0      | 0.3      | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | pinv           | Float64      | 2.1      | 0.38     | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| true    | krylov\_gpu    | Float64      | 5.9      | 0.4      | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | default\_multi | Float64      | 13.0     | 1.2      | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | intern         | Float64      | 13.0     | 1.7      | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |

### Float32

| **gpu** | **method**     | **el\_type** | **time**   | **GB**     | **percent\_X\_filled** | **sizeDesign** | **n\_channels** | **overlap** | **comment**         |
|--------:|---------------:|-------------:|-----------:|-----------:|-----------------------:|---------------:|----------------:|------------:|--------------------:|
| true    | cholesky       | Float32      |            |            | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  | PosDefException(-1) |
| true    | krylov\_gpu    | Float32      |            |            | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| true    | pinv           | Float32      | **0.39**   | 0.25       | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| true    | qr             | Float32      | 0.62       | 0.24       | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| true    | intern         | Float32      | 0.69       | 0.24       | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| true    | cg             | Float32      | 1.2        | 0.31       | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | cholesky       | Float32      | 1.2        | 0.17       | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | cg             | Float32      | 1.3        |**0.16**    | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | qr             | Float32      | 1.4        | 0.17       | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | pinv           | Float32      | 1.6        | 0.21       | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | intern         | Float32      | 13.0       | 0.86       | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |
| false   | default\_multi | Float32      | 13.0       | 0.97       | 0.0068                 | (239522, 1210) | 128             | (0.2, 0.2)  |                     |

## large, realistic model

```julia
    n_channels = 128,
    sfreq = 500,
    n_splines = (4, 4),
    n_repeats = 500,
```

| **gpu** | **method**     | **el\_type** | **time** | **GB**  | **percent\_X\_filled** | **sizeDesign**  | **n\_channels** | **overlap** | **comment**             |
|--------:|---------------:|-------------:|---------:|--------:|-----------------------:|----------------:|----------------:|------------:|------------------------:|
| true    | cholesky       | Float64      |          |         | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  | PosDefException(-1)     |
| false   | cholesky       | Float64      |          |         | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  | PosDefException(2760)   |
| false   | intern         | Float64      |          |         | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  | SingularException(9599) |
| true    | cg             | Float64      | **9.3**  | 3.6     | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  |                         |
| true    | qr             | Float64      | 11.0     | 3.5     | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  |                         |
| true    | intern         | Float64      | 13.0     | **3.5** | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  |                         |
| false   | qr             | Float64      | 80.0     | 6.3     | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  |                         |
| true    | pinv           | Float64      | 80.0     | 4.2     | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  |                         |
| true    | krylov\_gpu    | Float64      | 107.0    | 3.9     | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  |                         |
| false   | default\_multi | Float64      | 500.0    | 15.0    | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  |                         |
| false   | pinv           | Float64      | 520.0    | 11.0    | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  |                         |
| false   | cg             | Float64      | 939.0    | 5.7     | 0.0015                 | (3001479, 9616) | 128             | (0.2, 0.2)  |                         |
