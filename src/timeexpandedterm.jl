"""
Object with a *term* and an applicable *BasisFunction* and a eventfield that are later passed to the basisfunction.

$(TYPEDEF)
$(TYPEDSIGNATURES)
$(FIELDS)
# Examples
```julia-repl
julia>  b = TimeExpandedTerm(term,kernel,[:latencyTR,:durationTR])
```
"""
struct TimeExpandedTerm{T<:AbstractTerm} <: AbstractTerm 
    "Term that the basis function is applied to. This is regularly called in other functions to get e.g. term-coefnames and timeexpand those"
    term::T
    "Kernel that determines what should happen to the designmatrix of the term"
    basisfunction::BasisFunction
    "Which fields of the event-table should be passed to the basisfunction.Important: The first entry has to be the event-latency in samples!"
    eventfields::Array{Symbol}
end


function TimeExpandedTerm(term, basisfunction, eventfields::Nothing)
    # if the eventfield is nothing, call default
    TimeExpandedTerm(term, basisfunction)
end

function TimeExpandedTerm(term, basisfunction, eventfields::Symbol)
    # call with symbol, convert to array of symbol
    TimeExpandedTerm(term, basisfunction, [eventfields])
end

function TimeExpandedTerm(term, basisfunction; eventfields = [:latency])
    # default call, use latency
    TimeExpandedTerm(term, basisfunction, eventfields)
end

collabel(term::TimeExpandedTerm) = collabel(term.basisfunction)
StatsModels.width(term::TimeExpandedTerm) = width(term.basisfunction)
StatsModels.terms(t::TimeExpandedTerm) = terms(t.term)

function Base.show(io::IO, p::TimeExpandedTerm)
    print(io, "$(p.basisfunction.name): timeexpand($(p.term)) for times $(times(p.basisfunction))")
    
    #println(io, "$(coefnames(p))")
end

