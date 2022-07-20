using StatsModels
import Statistics.quantile

import Base.show
using Random
using BSplines
const bsPLINE_CONTEXT = Any



function offsetArrayToZeros!(oneRow, spl)
    oneRow[spl.offsets[1]+1:spl.offsets[1]+length(spl)] = parent(spl)
end

function genSplFunction(x, df)
    p = range(0.0, length = df - 1, stop = 1.0)
    breakpoints = quantile(x, p)
    basis = BSplineBasis(4, breakpoints) # 4 = cubic
    return x -> splFunction(x, basis)
end
function splFunction(x, basis)
    df = length(basis)
    
    large = zeros(Union{Missing,Float64},length(x), df)

    
    bs_eval = bsplines.(Ref(basis), x)
    
    
    for k = 1:length(bs_eval)
        if isnothing(bs_eval[k]) 
            @warn("spline prediction outside of possible range, putting to missing")
            large[k,:] .= missing
        else
            offsetArrayToZeros!(view(large, k, :), bs_eval[k])
        end
    end
    
    return large
end


#spl(x,df) = Splines2.bs(x,df=df,intercept=true) # assumes intercept
spl(x, df) = 1 # fallback

# make a nice call if the function is called via REPL
spl(t::Symbol, d::Int) = BSplineTerm(term(t), term(d))

mutable struct BSplineTerm{T,D} <: AbstractTerm
    term::T
    df::D
    fun::Any # function handle
end
function BSplineTerm(term, df)
    BSplineTerm(term, df, nothing)
end


Base.show(io::IO, p::BSplineTerm) = print(io, "spl($(p.term), $(p.df))")

function StatsModels.apply_schema(
    t::FunctionTerm{typeof(spl)},
    sch::StatsModels.Schema,
    Mod::Type{<:bsPLINE_CONTEXT},
)
    apply_schema(BSplineTerm(t.args_parsed...), sch, Mod)
end
function StatsModels.apply_schema(
    t::BSplineTerm,
    sch::StatsModels.Schema,
    Mod::Type{<:bsPLINE_CONTEXT},
)
    term = apply_schema(t.term, sch, Mod)
    isa(term, ContinuousTerm) ||
        throw(ArgumentError("BSplineTerm only works with continuous terms (got $term)"))

    isa(t.df, ConstantTerm) ||
        throw(ArgumentError("BSplineTerm df must be a number (got $(t.df))"))
    BSplineTerm(term, t.df.n)
end
function StatsModels.modelcols(p::BSplineTerm, d::NamedTuple)

    col = modelcols(p.term, d)


    if isnothing(p.fun)
        p.fun = genSplFunction(col, p.df)#Splines2.bs_(col,df=p.df+1,intercept=true)
    end
    #X = Splines2.bs(col, df=p.df+1,intercept=true)
    X = p.fun(col)

    # remove middle X to negate intercept = true, generating a pseudo effect code 
    return X[:, Not(Int(ceil(end / 2)))]
end

#StatsModels.terms(p::BSplineTerm) = terms(p.term)
StatsModels.termvars(p::BSplineTerm) = StatsModels.termvars(p.term)
StatsModels.width(p::BSplineTerm) = p.df
StatsModels.coefnames(p::BSplineTerm) =
    "spl(" .* coefnames(p.term) .* "," .* string.(1:p.df) .* ")"
