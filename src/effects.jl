import Effects: effects
import Effects: expand_grid
import Effects: typify
import Effects.typify
import Effects: _symequal
import StatsModels.collect_matrix_terms
import Base.getproperty
using Effects


"""
effects(design::AbstractDict, model::UnfoldModel;typical=mean)

Calculates marginal effects for all term-combinations in `design`.

 Implementation based on Effects Package; likely could repackage in UnfoldEffects; somebody wants to do it? This would make it easier to cross-maintain it to changes/bugfixes in the Effects.jl Package
 `design` is a Dictionary containing those predictors (as keys) with levels (as values), that you want to evaluate. The `typical` refers to the value, that other predictors not in the Dictionary should take on.


For MixedModels, the returned effects are based on the "typical" subject, i.e. all random effects are put to 0.

 # Example
 ```julia-repl
 julia> f = @formula 0 ~ categoricalA + continuousA + continuousB
 julia> uf = fit(UnfoldModel,(Any=>(f,times)),data,events)
 julia> d = Dict(:categoricalA=>["levelA","levelB"],:continuousB=>[-2,0,2])
 julia> effects(d,uf)
```
 will result in 6 predicted values: A/-2, A/0, A/2, B/-2, B/0, B/2.
"""


function effects(design::AbstractDict, model::T; typical = mean) where {T<:UnfoldModel}
    if isempty(design)
        return effects(Dict(:dummy=>[:dummy]),model;typical)
    end
    reference_grid = expand_grid(design)
    form = Unfold.formulas(model) # get formula

    # replace non-specified fields with "constants"
    m = Unfold.modelmatrix(model, false) # get the modelmatrix without timeexpansion
    #@debug "type form[1]", typeof(form[1])

    form_typical = _typify(T, reference_grid, form, m, typical)
    @debug typeof(form_typical) typeof(form_typical[1])

    #form_typical = vec(form_typical)
    reference_grids = repeat([reference_grid], length(form_typical))

    eff = predict(model, form_typical, reference_grids; overlap = false)
    if :latency ∈ unique(vcat(names.(reference_grids)...))
        reference_grids = select.(reference_grids, Ref(DataFrames.Not(:latency)))
    end
    @debug "effects" size(eff[1]) reference_grid size(times(model)[1]) eventnames(model)
    return result_to_table(eff, reference_grids, times(model), eventnames(model))


end

Effects.typify(reference_grid, form::AbstractArray, X; kwargs...) =
    typify.(Ref(reference_grid), form, Ref(X); kwargs...)


# cast single form to a vector
_typify(
    reference_grid,
    form::FormulaTerm{<:InterceptTerm,<:Unfold.TimeExpandedTerm},
    m,
    typical,
) = _typify(reference_grid, [form], [m], typical)

@traitfn function _typify(
    ::Type{UF},
    reference_grid,
    form::Vector{T},
    m::Vector,
    typical,
) where {T<:FormulaTerm,UF<:UnfoldModel;ContinuousTimeTrait{UF}}
    @debug "_typify - stripping away timeexpandedterm"
    form_typical = Array{FormulaTerm}(undef, length(form))
    for f = 1:length(form)

        # strip of basisfunction and put it on afterwards again
        tmpf = deepcopy(form[f])

        # create a Formula without Timeexpansion
        tmpf = FormulaTerm(tmpf.lhs, tmpf.rhs.term)

        # typify that
        tmpf = typify(reference_grid, tmpf, m[f]; typical = typical)

        # regenerate TimeExpansion
        tmpf = Unfold.TimeExpandedTerm(
            tmpf,
            form[f].rhs.basisfunction;
            eventfields = form[f].rhs.eventfields,
        )
        form_typical[f] = FormulaTerm(form[f].lhs, tmpf)
    end
    return form_typical
end
function _typify(reference_grid, form::FormulaTerm, m, typical)

    return [typify(reference_grid, form, m; typical = typical)]

end

@traitfn function _typify(
    ::Type{UF},
    reference_grid,
    form::AbstractArray{T},
    m::Vector,
    typical,
) where {T<:FormulaTerm,UF<:UnfoldModel;!ContinuousTimeTrait{UF}}
    # Mass Univariate with multiple effects
    @debug "_typify going the mass univariate route - $(typeof(form))"
    @debug length(form), length(m), typeof(m)
    out = FormulaTerm[]
    for k = 1:length(form)
        tmpf = typify(reference_grid, form[k], m[k]; typical = typical)
        push!(out, FormulaTerm(form[k].lhs, tmpf))
    end
    @debug :_typify typeof(form)
    return out

end
# mixedModels case - just use the FixEff, ignore the ranefs
Effects.typify(reference_grid, form, m::Tuple; typical) =
    typify(reference_grid, form, m[1]; typical)
