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

minus one due to the intercept.

Note: that due to the boundary condition (`natural`) spline, we repeat the boundary knots to each side `order` times, enforcing smoothness there - this is done within BSplineKit

"""
function genSpl_breakpoints(p::AbstractSplineTerm, x)
    p = range(0.0, length = p.df - 2, stop = 1.0)
    breakpoints = quantile(x, p)
    return breakpoints
end

"""
In the circular case, we do not use quantiles, (circular quantiles are difficult)
"""
function genSpl_breakpoints(p::PeriodicBSplineTerm, x)
    # periodic case - 
    return range(p.low, p.high, length = p.df + 2)
end

"""
function that fills in an Matrix `large` according to the evaluated values in `x` with the designmatrix values for the spline.

    Two separate functions are needed here, as the periodicBSplineBasis implemented in BSplineKits is a bit weird, that it requires to evalute "negative" knots + knots above the top-boundary and fold them down
"""
function spl_fillMat!(bs::PeriodicBSplineBasis, large::Matrix, x::AbstractVector)
    # wrap values around the boundaries
    bnds = boundaries(bs)
    x = deepcopy(x)
    x = mod.(x .- bnds[1], period(bs)) .+ bnds[1]

    for k = -1:length(bs)+2
        ix = basis_to_array_index(bs, axes(large, 2), k)
        large[:, ix] .+= bs[k](x)
    end
end
function spl_fillMat!(bs::BSplineBasis, large::Matrix, x::AbstractVector)
    for k = 1:length(bs)

        large[:, k] .+= bs[k](x)
    end

    bnds = boundaries(bs)
    ix = x .< bnds[1] .|| x .> bnds[2]

    if sum(ix) != 0
        @warn(
            "spline prediction outside of possible range  putting those values to missing.\n `findfirst(Out-Of-Bound-value)` is x=$(x[findfirst(ix)]), with bounds: $bnds"
        )
        large = allowmissing(large)
        large[ix, :] .= missing
    end

end

"""
    splFunction(x::AbstractVector,bs::AbstractBSplineTerm)
    _splFunction(x::AbstractVector,bs::AbstractBSplineTerm)
evaluate a spline basisset `basis` at `x`. Automatically promotes `x` to Float64 if not AbstractFloat

returns `Missing` if x is outside of the basis set
"""
function _splFunction(x::AbstractVector{T}, bs) where {T<:AbstractFloat}
    @debug "spl" typeof(x)
    # init array
    large = zeros(T, length(x), length(bs))

    # fill it with spline values
    spl_fillMat!(bs, large, x)

    return large
end
_splFunction(x::AbstractVector, bs) = _splFunction(Float64.(x), bs)

splFunction(x, bs) = _splFunction(x, bs)

function splFunction(x::AbstractVector, spl::PeriodicBSplineTerm)
    basis = PeriodicBSplineBasis(BSplineOrder(spl.order), deepcopy(spl.breakpoints))
    _splFunction(x, basis)
end

function splFunction(x::AbstractVector, spl::BSplineTerm)
    basis = BSplineKit.BSplineBasis(BSplineOrder(spl.order), deepcopy(spl.breakpoints))
    _splFunction(x, basis)
end
#spl(x,df) = Splines2.bs(x,df=df,intercept=true) # assumes intercept
Unfold.spl(x, df) = 0 # fallback

# make a nice call if the function is called via REPL
Unfold.spl(t::Symbol, d::Int) = BSplineTerm(term(t), d, 4, [])
Unfold.circspl(t::Symbol, d::Int, low, high) =
    PeriodicBSplineTerm(term(t), term(d), 4, low, high)

"""
Construct a BSplineTerm, if breakpoints/basis are not defined yet, put to `nothing`
"""
function BSplineTerm(term, df, order = 4)
    @assert df > 3 "Minimal degrees of freedom has to be 4"
    BSplineTerm(term, df, order, [])
end

function BSplineTerm(term, df::ConstantTerm, order = 4)
    BSplineTerm(term, df.n, order, [])
end


function PeriodicBSplineTerm(term, df, low, high)
    PeriodicBSplineTerm(term, df, 4, low, high)
end
function PeriodicBSplineTerm(
    term::AbstractTerm,
    df::ConstantTerm,
    order,
    low::ConstantTerm,
    high::ConstantTerm,
    breakvec,
)
    PeriodicBSplineTerm(term, df.n, order, low.n, high.n, breakvec)
end
function PeriodicBSplineTerm(term, df, order, low, high)
    PeriodicBSplineTerm(term, df, order, low, high, [])
end

Base.show(io::IO, p::BSplineTerm) = print(io, "spl($(p.term), $(p.df))")
Base.show(io::IO, p::PeriodicBSplineTerm) =
    print(io, "circspl($(p.term), $(p.df),$(p.low):$(p.high))")

function StatsModels.apply_schema(
    t::FunctionTerm{typeof(Unfold.spl)},
    sch::StatsModels.Schema,
    Mod::Type{<:bsPLINE_CONTEXT},
)
    @debug "BSpline spl Schema"
    ar = nothing
    try
        ar = t.args
    catch
        ar = t.args_parsed # statsmodels < 0.7
    end
    apply_schema(BSplineTerm(ar...), sch, Mod)
end

function StatsModels.apply_schema(
    t::FunctionTerm{typeof(Unfold.circspl)},
    sch::StatsModels.Schema,
    Mod::Type{<:bsPLINE_CONTEXT},
)
    ar = nothing
    try
        ar = t.args
    catch
        ar = t.args_parsed # statsmodels < 0.7
    end
    apply_schema(PeriodicBSplineTerm(ar...), sch, Mod)
end
function StatsModels.apply_schema(
    t::AbstractSplineTerm,
    sch::StatsModels.Schema,
    Mod::Type{<:bsPLINE_CONTEXT},
)
    @debug "BSpline Inner schema"
    term = apply_schema(t.term, sch, Mod)
    isa(term, ContinuousTerm) ||
        throw(ArgumentError("BSplineTerm only works with continuous terms (got $term)"))

    if isa(t.df, ConstantTerm)
        try
            # in case of ConstantTerm of Èffects.jl``
            t.df.n
        catch
            throw(ArgumentError("BSplineTerm df must be a number (got $(t.df))"))
        end
    end
    return construct_spline(t, term)
end
construct_spline(t::BSplineTerm, term) = BSplineTerm(term, t.df, t.order)
construct_spline(t::PeriodicBSplineTerm, term) =
    PeriodicBSplineTerm(term, t.df, t.order, t.low, t.high)

function StatsModels.modelcols(p::AbstractSplineTerm, d::NamedTuple)

    col = modelcols(p.term, d)

    if isempty(p.breakpoints)
        p.breakpoints = genSpl_breakpoints(p, col)
    end

    #basis = genSpl_basis(pp.breakpoints,p.order)#Splines2.bs_(col,df=p.df+1,intercept=true)

    #X = Splines2.bs(col, df=p.df+1,intercept=true)
    X = splFunction(col, p)

    # remove middle X to negate intercept = true, generating a pseudo effect code 
    return X[:, Not(Int(ceil(end / 2)))]
end

#StatsModels.terms(p::BSplineTerm) = terms(p.term)
StatsModels.termvars(p::AbstractSplineTerm) = StatsModels.termvars(p.term)
StatsModels.width(p::AbstractSplineTerm) = p.df - 1
StatsModels.coefnames(p::BSplineTerm) =
    "spl(" .* coefnames(p.term) .* "," .* string.(1:p.df-1) .* ")"
StatsModels.coefnames(p::PeriodicBSplineTerm) =
    "circspl(" .* coefnames(p.term) .* "," .* string.(1:p.df-1) .* ",$(p.low):$(p.high))"
