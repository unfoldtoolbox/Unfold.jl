module UnfoldBSplineKitExt
using Unfold
    using BSplineKit
    using StatsModels
import StatsModels: termvars,width,coefnames,modelcols,apply_schema
import Base: show
using Effects

abstract type AbstractSplineTerm <:AbstractTerm end

    include("basisfunctions.jl")
    include("splinepredictors.jl")


## Effects
Effects._trmequal(t1::Unfold.AbstractSplineTerm,t2::AbstractTerm) = _symequal(t1.term,t2)
Effects._trmequal(t1::Unfold.AbstractSplineTerm,t2::Unfold.AbstractSplineTerm) = _symequal(t1.term,t2.term)

Effects._trmequal(t1::AbstractTerm,t2::Unfold.AbstractSplineTerm) = _symequal(t1,t2.term)
Effects._symequal(t1::AbstractTerm,t2::Unfold.AbstractSplineTerm) = _symequal(t1,t2.term)
end