
const bsPLINE_CONTEXT = Any

mutable struct BSplineTerm{T,D} <: AbstractSplineTerm
    term::T
    df::D
    order::Int
    breakpoints::Vector
end
mutable struct PeriodicBSplineTerm{T,D} <: AbstractSplineTerm
    term::T
    df::D
    order::Int
    low::Real
    high::Real
    breakpoints::Vector
end

"""
generate a cubic spline basis function set, with df-1 breakpoints fixed to the quantiles of x
"""
function genSpl_breakpoints(x, df)
    p = range(0.0, length = df - 1, stop = 1.0)
    breakpoints = quantile(x, p)
    return breakpoints
end

"""
In the circular case, we do not use quantiles, (circular quantiles are difficult)
"""
function genSpl_breakpoints(x,df,lo,hi)
    # periodic case - 
    return range(lo,hi,length=df-1)
end

"""
evaluate a spline basisset `basis` at `x`

returns `Missing` if x is outside of the basis set
"""
function splFunction(x, bs)
    df = length(bs)
    
    large = zeros(Union{Missing,Float64},length(x), df)
    
    #bs_eval = evaluate_all.(basis, x)
    bnds = boundaries(bs)
    for k = 1:df
        #@show k
        large[:,k] .= bs[k](x)
    end

    ix = x .< bnds[1] .|| x .>bnds[2]
    
    if sum(ix) != 0
        @warn("spline prediction outside of possible range, putting those values to missing")
        large[ix,:] .= missing
    end

    
    return large
end

function splFunction(x,spl::PeriodicBSplineTerm)
    basis = PeriodicBSplineBasis(BSplineOrder(spl.order),spl.breakpoints)
    splFunction(x,basis)
end

function splFunction(x,spl::BSplineTerm)
    basis = BSplineKit.BSplineBasis(BSplineOrder(spl.order),spl.breakpoints)
    splFunction(x,basis)
end
#spl(x,df) = Splines2.bs(x,df=df,intercept=true) # assumes intercept
spl(x, df) = 1 # fallback

# make a nice call if the function is called via REPL
spl(t::Symbol, d::Int) = BSplineTerm(term(t), term(d),BSplineBasis,4)
circspl(t::Symbol, d::Int,low,high) = PeriodicBSplineTerm(term(t), term(d),4,low,high)

"""
Construct a BSplineTerm, if breakpoints/basis are not defined yet, put to `nothing`
"""
function BSplineTerm(term, df,order=4)
    BSplineTerm(term, df,order,[])
end

function PeriodicBSplineTerm(term, df,low,high)
    PeriodicBSplineTerm(term, df,4,low,high,[])
end
function PeriodicBSplineTerm(term, df,order,low,high)
    PeriodicBSplineTerm(term, df,order,low,high,[])
end

Base.show(io::IO, p::BSplineTerm) = print(io, "spl($(p.term), $(p.df))")
Base.show(io::IO, p::PeriodicBSplineTerm) = print(io, "circspl($(p.term), $(p.df),$(p.low):$(p.high))")

function StatsModels.apply_schema(
    t::FunctionTerm{typeof(spl)},
    sch::StatsModels.Schema,
    Mod::Type{<:bsPLINE_CONTEXT},
)
    apply_schema(BSplineTerm(t.args_parsed...), sch, Mod)
end

function StatsModels.apply_schema(
    t::FunctionTerm{typeof(circspl)},
    sch::StatsModels.Schema,
    Mod::Type{<:bsPLINE_CONTEXT},
)
    apply_schema(PeriodicBSplineTerm(t.args_parsed...), sch, Mod)
end
function StatsModels.apply_schema(
    t::AbstractSplineTerm,
    sch::StatsModels.Schema,
    Mod::Type{<:bsPLINE_CONTEXT},
)
    term = apply_schema(t.term, sch, Mod)
    isa(term, ContinuousTerm) ||
        throw(ArgumentError("BSplineTerm only works with continuous terms (got $term)"))

    isa(t.df, ConstantTerm) ||
        throw(ArgumentError("BSplineTerm df must be a number (got $(t.df))"))
    return construct_spline(t,term)
    end
construct_spline(t::BSplineTerm,term)=BSplineTerm(term, t.df.n,t.order)
construct_spline(t::PeriodicBSplineTerm,term)=PeriodicBSplineTerm(term, t.df.n,t.order,t.low,t.high)

function StatsModels.modelcols(p::AbstractSplineTerm, d::NamedTuple)

    col = modelcols(p.term, d)

    if isempty(p.breakpoints)
        p.breakpoints = genSpl_breakpoints(col,p.df)
    end
    
    #basis = genSpl_basis(pp.breakpoints,p.order)#Splines2.bs_(col,df=p.df+1,intercept=true)
    
    #X = Splines2.bs(col, df=p.df+1,intercept=true)
    X = splFunction(col,p)

    # remove middle X to negate intercept = true, generating a pseudo effect code 
    return X[:, Not(Int(ceil(end / 2)))]
end

#StatsModels.terms(p::BSplineTerm) = terms(p.term)
StatsModels.termvars(p::AbstractSplineTerm) = StatsModels.termvars(p.term)
StatsModels.width(p::AbstractSplineTerm) = p.df
StatsModels.coefnames(p::BSplineTerm) =
    "spl(" .* coefnames(p.term) .* "," .* string.(1:p.df) .* ")"
StatsModels.coefnames(p::PeriodicBSplineTerm) =
    "circspl(" .* coefnames(p.term) .* "," .* string.(1:p.df) .* ",$(p.low):$(p.high))"
