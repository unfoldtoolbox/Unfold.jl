using StatsModels
import Splines2.bs
const bsPLINE_CONTEXT = Any
bs(x,df) = Splines2.bs(x,df=df,intercept=true) # assumes intercept
struct bsplineTerm{T,D} <: AbstractTerm
    term::T
    df::D
end
Base.show(io::IO, p::bsplineTerm) = print(io, "bs($(p.term), $(p.df))")
function StatsModels.apply_schema(t::FunctionTerm{typeof(bs)},
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:bsPLINE_CONTEXT})
    apply_schema(bsplineTerm(t.args_parsed...), sch, Mod)
end
function StatsModels.apply_schema(t::bsplineTerm,
                                  sch::StatsModels.Schema,
                                  Mod::Type{<:bsPLINE_CONTEXT})
    term = apply_schema(t.term, sch, Mod)
    isa(term, ContinuousTerm) ||
        throw(ArgumentError("bsplineTerm only works with continuous terms (got $term)"))
    isa(t.df, ConstantTerm) ||
        throw(ArgumentError("bsplineTerm df must be a number (got $(t.df))"))
    bsplineTerm(term, t.df.n)
end
function StatsModels.modelcols(p::bsplineTerm, d::NamedTuple)
    col = modelcols(p.term, d)
    X = Splines2.bs(col, df=p.df+1,intercept=true)
    # remove middle X to negate intercept = true, generating a pseudo effect code 
    X[:,Not(Int(ceil(end/2)))]
end
StatsModels.terms(p::bsplineTerm) = terms(p.term)
StatsModels.termvars(p::bsplineTerm) = StatsModels.termvars(p.term)
StatsModels.width(p::bsplineTerm) = 1
StatsModels.coefnames(p::bsplineTerm) = "bs(" .* coefnames(p.term) .* "," .* string.(1:p.df) .* ")"