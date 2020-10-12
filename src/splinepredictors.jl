using StatsModels
import Splines2.bs
import Base.show
const bsPLINE_CONTEXT = Any
spl(x,df) = Splines2.bs(x,df=df,intercept=true) # assumes intercept
mutable struct uf_bsplineTerm{T,D} <: AbstractTerm
    term::T
    df::D
    fun # function handle
end
function uf_bsplineTerm(term,df)
    uf_bsplineTerm(term,df,nothing)
end


Base.show(io::IO, p::uf_bsplineTerm) = print(io, "spl($(p.term), $(p.df))")

function StatsModels.apply_schema(t::FunctionTerm{typeof(spl)},
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:bsPLINE_CONTEXT})
    apply_schema(uf_bsplineTerm(t.args_parsed...), sch, Mod)
end
function StatsModels.apply_schema(t::uf_bsplineTerm,
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:bsPLINE_CONTEXT})
    term = apply_schema(t.term, sch, Mod)
    isa(term, ContinuousTerm) ||
        throw(ArgumentError("uf_bsplineTerm only works with continuous terms (got $term)"))

    isa(t.df, ConstantTerm) ||
            throw(ArgumentError("uf_bsplineTerm df must be a number (got $(t.df))"))
    uf_bsplineTerm(term, t.df.n)
end
function StatsModels.modelcols(p::uf_bsplineTerm, d::NamedTuple)
    col = modelcols(p.term, d)
    if isnothing(p.fun)
        p.fun = Splines2.bs_(col,df=p.df+1,intercept=true)
    end
        #X = Splines2.bs(col, df=p.df+1,intercept=true)
    X = p.fun(col)
    # remove middle X to negate intercept = true, generating a pseudo effect code 
    X[:,Not(Int(ceil(end/2)))]
end
StatsModels.terms(p::uf_bsplineTerm) = terms(p.term)
StatsModels.termvars(p::uf_bsplineTerm) = StatsModels.termvars(p.term)
StatsModels.width(p::uf_bsplineTerm) = 1
StatsModels.coefnames(p::uf_bsplineTerm) = "spl(" .* coefnames(p.term) .* "," .* string.(1:p.df) .* ")"