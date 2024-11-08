module UnfoldMixedModelsExt

using Unfold

import Unfold: isa_lmm_formula, make_estimate, modelmatrices
using MixedModels
#import MixedModels.FeMat # extended for sparse femats, type piracy => issue on MixedModels.jl github
using StaticArrays # for MixedModels extraction of parametrs (inherited from MixedModels.jl, not strictly needed )
import MixedModels: likelihoodratiotest, ranef
using StatsModels
import StatsModels: fit!, coef, coefnames, modelcols, modelmatrix
using SparseArrays
using DocStringExtensions
using LinearAlgebra # LowerTriangular
using DataFrames
using ProgressMeter
using SimpleTraits

include("typedefinitions.jl")
include("condense.jl")
include("designmatrix.jl")
include("fit.jl")
include("statistics.jl")
include("timeexpandedterm.jl")
include("effects.jl")

end
