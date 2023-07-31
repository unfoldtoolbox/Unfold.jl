module UnfoldMixedModelsExt

using Unfold

import Unfold: isMixedModelFormula,make_estimate
using MixedModels
#import MixedModels.FeMat # extended for sparse femats, type piracy => issue on MixedModels.jl github
using StaticArrays # for MixedModels extraction of parametrs (inherited from MixedModels.jl, not strictly needed )
import MixedModels: likelihoodratiotest,ranef
using StatsModels
import StatsModels: fit!,coef,coefnames,modelcols
using SparseArrays
using DocStringExtensions
using LinearAlgebra # LowerTriangular

include("condense.jl")
include("designmatrix.jl")
include("fit.jl")
include("statistics.jl")
include("timeexpandedterm.jl")
include("typedefinitions.jl")

end