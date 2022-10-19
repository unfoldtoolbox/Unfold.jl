using StatsBase
using StatsModels
import Statistics.quantile
#import Splines2.bs
import Base.show
using Random
using LinearAlgebra

function genCircSplFunction(x, df, bounds)
    x = mapValues(x, bounds[1], bounds[2])
    p = range(0.0, length = df + 2, stop = 1.0)
    knots = quantile(x, p)
    basis = genCircSplBasis(x, knots)
    return y -> circSplFunction(mapValues(y,bounds[1],bounds[2]), basis, knots)
end

function mapValues(x, lowerBound, upperBound)
    if lowerBound >= upperBound
        error("the lower bound has to be smaller than the upper bound.")
    end


    x[x .> upperBound] .= lowerBound .+ (x[x .> upperBound] .- upperBound) .% (upperBound - lowerBound)
    x[x .< lowerBound] .= upperBound .- (lowerBound .- x[x .< lowerBound]) .% (upperBound - lowerBound)
    #for (index, value) in enumerate(x)
    #    if value > upperBound
    #        x[index] = lowerBound + (value - upperBound)
    #    elseif value < lowerBound
    #        x[index] = upperBound - (lowerBound - value)
    #    end
    #end
    return x
end

function genCircSplBasis(x, knots)
    j = findknotslowerbounds(x, knots)
    h = knots[2:length(knots)] - knots[1:length(knots)-1]
    hj = h[j]
    xj1_x = knots[j.+1] - x
    x_xj = x - knots[j]

    ajm = xj1_x ./ hj
    ajp = x_xj ./ hj

    cjm_3 = xj1_x .* xj1_x .* xj1_x ./ (6. .* hj)
    cjm_3[x .> max(knots...)] .= 0.
    cjm_1 = hj .* xj1_x ./ 6
    cjm = cjm_3 - cjm_1

    cjp_3 = x_xj .* x_xj .* x_xj ./ (6. .* hj)
    cjp_3[x .< min(knots...)] .= 0
    cjp_1 = hj .* x_xj ./ 6
    cjp = cjp_3 - cjp_1

    return Dict("ajm" => ajm, "ajp" => ajp, "cjm" => cjm, "cjp" => cjp, "j" => j)
end

function findknotslowerbounds(x, knots)
    lb = zeros(Int64, length(x))

    for xi in eachindex(x)
        tmp = findfirst(knots .> x[xi])
        lb[xi] = isnothing(tmp) ? length(knots) : tmp - 1
    end

    lb[lb .== 0] .= 1
    lb[lb .== length(knots)] .= length(knots) - 1

    return lb
end

function circSplFunction(x, basis, knots)
    n = length(knots) - 1
    j1 = basis["j"] .+ 1
    j1[j1 .== n+1] .= 1
    i = Matrix{Int}(I, n, n)
    f = getCyclicF(knots)
    return basis["ajm"] .* i[basis["j"],:] + basis["ajp"] .* i[j1,:] + basis["cjm"] .* f[basis["j"],:] + basis["cjp"] .* f[j1,:]
end

function getCyclicF(knots)
    h = knots[2:length(knots)] - knots[1:length(knots)-1]
    n = length(knots) - 1
    b = zeros(n, n)
    d = zeros(n, n)

    b[1, 1] = (h[n - 1] + h[1]) ./ 3.
    b[1, n] = h[n - 1] ./ 6.
    b[n, 1] = h[n - 1] ./ 6.

    d[1, 1] = -1. ./ h[1] - 1. ./ h[n - 1]
    d[1, n] = 1. ./ h[n - 1]
    d[n, 1] = 1. ./ h[n - 1]

    for i in 2:n
        b[i, i] = (h[i - 1] + h[i]) ./ 3.
        b[i, i - 1] = h[i - 1] ./ 6.
        b[i - 1, i] = h[i - 1] ./ 6.
        
        d[i, i] = -1. ./ h[i - 1] - 1. ./ h[i]
        d[i, i - 1] = 1. ./ h[i - 1]
        d[i - 1, i] = 1. ./ h[i - 1]
    end
    return (b\d)
end

# type of model where syntax applies: here this applies to any model type
const CIRCSPL_CONTEXT = Any

# syntax: best practice to define a _new_ function
circspl(x, df, lbound, ubound) = 1

# make a nice call if the function is called via REPL
circspl(t::Symbol, d::Int, l::Int, u::Int) = uf_circsplineTerm(term(t), term(d), term(l), term(u))

# struct for behavior
mutable struct uf_circSplTerm{T,D,L,U} <: AbstractTerm
    term::T
    deg::D
    lb::L
    ub::U
    fun::Any # function handle
end

function uf_circSplineTerm(term, df, lb, ub)
    uf_circSplTerm(term, df, lb, ub, nothing)
end

Base.show(io::IO, p::uf_circSplTerm) = print(io, "circspl($(p.term), $(p.deg), $(p.lb), $(p.ub))")

# for `circspl` use at run-time (outside @formula), return a schema-less uf_circSplTerm
#circspl(t::Symbol, d::Int, b::Vector{Int64}) = uf_circSplTerm(term(t), term(d), term(b))

# for `circspl` use inside @formula: create a schemaless uf_circSplTerm and apply_schema
function StatsModels.apply_schema(
    t::FunctionTerm{typeof(circspl)},
    sch::StatsModels.Schema,
    Mod::Type{<:CIRCSPL_CONTEXT}
)
    apply_schema(uf_circSplineTerm(t.args_parsed...), sch, Mod)
end

# apply_schema to internal Terms and check for proper types
function StatsModels.apply_schema(
    t::uf_circSplTerm,
    sch::StatsModels.Schema,
    Mod::Type{<:CIRCSPL_CONTEXT}
)
    term = apply_schema(t.term, sch, Mod)
    isa(term, ContinuousTerm) ||
        throw(ArgumentError("uf_circSplineTerm only works with continuous terms (got $term)"))
    isa(t.deg, ConstantTerm) ||
        throw(ArgumentError("uf_circSplineTerm df must be a number (got $t.deg)"))
    isa(t.lb, ConstantTerm) ||
        throw(ArgumentError("uf_circSplineTerm lb must be a number (got $t.lb)"))
    isa(t.ub, ConstantTerm) ||
        throw(ArgumentError("uf_circSplineTerm ub must be a number (got $t.ub)"))
    
    uf_circSplineTerm(term, t.deg.n, t.lb.n, t.ub.n)
end

function StatsModels.modelcols(p::uf_circSplTerm, d::NamedTuple)
    col = modelcols(p.term, d)
    if isnothing(p.fun)
        p.fun = genCircSplFunction(col, p.deg, [p.lb, p.ub])#Splines2.bs_(col,df=p.df+1,intercept=true)
    end
    #X = Splines2.bs(col, df=p.df+1,intercept=true)
    X = p.fun(col)

    # remove middle X to negate intercept = true, generating a pseudo effect code 
    X[:, Not(Int(ceil(end / 2)))]
end

# the basic terms contained within a uf_circSplTerm (for schema extraction)
#StatsModels.terms(p::uf_circSplTerm) = terms(p.term)
# names variables from the data that a uf_circSplTerm relies on
StatsModels.termvars(p::uf_circSplTerm) = StatsModels.termvars(p.term)
# number of columns in the matrix this term produces
# QUESTION: this has to be deg - 1 for circulat splines, no?
StatsModels.width(p::uf_circSplTerm) = p.deg
#StatsBase.coefnames(p::uf_circSplTerm) = "circspl(" .* coefnames(p.term) .* "," .* string.(1:p.deg) .* "," .* string.(p.lb) .* "," .* string.(p.ub) .* ")"
StatsBase.coefnames(p::uf_circSplTerm) = coefnames(p.term) .* "^" .* string.(1:p.deg)