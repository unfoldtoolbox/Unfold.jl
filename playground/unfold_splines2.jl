using StatsModels
import Splines2.bs
import Base.show
const bsPLINE_CONTEXT = Any
bs2(x,df) = Splines2.bs(x,df=df,intercept=true) # assumes intercept
mutable struct uf_bsplineTerm2{T,D} <: AbstractTerm
    term::T
    df::D
    fun # function handle
end
function uf_bsplineTerm2(term,df)
    uf_bsplineTerm2(term,df,nothing)
end


Base.show(io::IO, p::uf_bsplineTerm2) = print(io, "bs2($(p.term), $(p.df))")

function StatsModels.apply_schema(t::FunctionTerm{typeof(bs2)},
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:bsPLINE_CONTEXT})
    apply_schema(uf_bsplineTerm2(t.args_parsed...), sch, Mod)
end
function StatsModels.apply_schema(t::uf_bsplineTerm2,
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:bsPLINE_CONTEXT})
    term = apply_schema(t.term, sch, Mod)
    isa(term, ContinuousTerm) ||
        throw(ArgumentError("uf_bsplineTerm2 only works with continuous terms (got $term)"))
    isa(t.df, ConstantTerm) ||
        throw(ArgumentError("uf_bsplineTerm2 df must be a number (got $(t.df))"))
    uf_bsplineTerm2(term, t.df.n)
end
function StatsModels.modelcols(p::uf_bsplineTerm2, d::NamedTuple)
    col = modelcols(p.term, d)
    if isnothing(p.fun)
        p.fun = Splines2.bs_(col,df=p.df+1,intercept=true)
    end
        #X = Splines2.bs(col, df=p.df+1,intercept=true)
    X = p.fun(col)
    # remove middle X to negate intercept = true, generating a pseudo effect code 
    X[:,Not(Int(ceil(end/2)))]
end
StatsModels.terms(p::uf_bsplineTerm2) = terms(p.term)
StatsModels.termvars(p::uf_bsplineTerm2) = StatsModels.termvars(p.term)
StatsModels.width(p::uf_bsplineTerm2) = 1
StatsModels.coefnames(p::uf_bsplineTerm2) = "bs2(" .* coefnames(p.term) .* "," .* string.(1:p.df) .* ")"