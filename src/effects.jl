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
    reference_grid = expand_grid(design)
    form = Unfold.formulas(model) # get formula

    # replace non-specified fields with "constants"
    m = Unfold.modelmatrix(model, false) # get the modelmatrix without timeexpansion
    #@debug "type form[1]", typeof(form[1])

    form_typical = _typify(T, reference_grid, form, m, typical)

    form_typical = vec(form_typical)
    reference_grids = repeat([reference_grid], length(form_typical))


    eff = predict(model, form_typical, reference_grids; overlap = false)

    return predict_to_table(model, eff, form_typical, reference_grids)


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
    form::AbstractArray{<:FormulaTerm},
    m::Vector,
    typical,
) where {UF <: UnfoldModel; !ContinuousTimeTrait{UF}}
    # Mass Univariate with multiple effects
    @debug "_typify going the mass univariate route - $(typeof(form))"
    @debug length(form), length(m), typeof(m)
    out = []
    for k = 1:length(form)
        push!(out, typify(reference_grid, form[k], m[k]; typical = typical))
    end
    return out

end
# mixedModels case - just use the FixEff, ignore the ranefs
Effects.typify(reference_grid, form, m::Tuple; typical) =
    typify(reference_grid, form, m[1]; typical)

function cast_referenceGrid(r, eff::AbstractArray{T}, times; eventname = nothing) where {T}
    @debug typeof(eff), typeof(r)
    nchan = size(eff, 2) # correct
    neff = size(r, 1) # how many effects requested
    neffCol = size(r, 2) # how many predictors
    ncols = size(eff, 1) ÷ neff # typically ntimes



    # replicate
    # for each predictor in r (reference grid), we need this at the bottom

    if isnothing(eventname)
        nbases = 1
    else
        nbases = length(unique(eventname))
    end
    coefs_rep = Array{Array}(undef, nbases, neffCol)


    for k = 1:neffCol
        # in case we have only a single basis (e.g. mass univariate), we can directly fill in all values
        ixList = []
        if isnothing(eventname)
            ix = ones(ncols) .== 1.0
            append!(ixList, [ix])
        else
            #in case of multiple bases, we have to do it iteratively, because the bases can be different length
            for b in unique(eventname)
                ix = eventname[1:neff:end] .== b
                append!(ixList, [ix])
            end
        end
        for i_ix = 1:length(ixList)

            coefs_rep[i_ix, k] = linearize(
                permutedims(
                    repeat(r[:, k], outer = [1, nchan, sum(ixList[i_ix])]),
                    [2, 3, 1],
                ),
            )
        end
    end

    # often the "times" vector
    if length(times) == neff * ncols
        # in case we have timeexpanded, times is already in long format and doesnt need to be repeated for each coefficient
        colnames_basis_rep = permutedims(repeat(times, 1, nchan, 1), [2 1 3])
    else
        colnames_basis_rep = permutedims(repeat(times, 1, nchan, neff), [2 1 3])
    end

    # for multiple channels
    chan_rep = repeat(1:nchan, 1, ncols, neff)

    # for mass univariate there is no eventname
    if isnothing(eventname)
        eventname = fill(nothing, ncols)
        eventname_rep = permutedims(repeat(eventname, 1, nchan, neff), [2, 1, 3])
    else
        eventname_rep = permutedims(repeat(eventname, 1, nchan, 1), [2, 1, 3])
    end


    result = Dict(
        :yhat => linearize(eff'),
        :time => linearize(colnames_basis_rep),
        :channel => linearize(chan_rep),
        :eventname => linearize(eventname_rep),
    )
    #@debug size(coefs_rep), typeof(coefs_rep), size(coefs_rep[1, 1])
    #@debug coefs_rep[1, 2][1:2]
    for k = 1:neffCol
        #@debug names(r)[k]
        push!(result, Symbol(names(r)[k]) => reduce(vcat, coefs_rep[:, k]))
    end

    return result
end

#Effects._symequal(t1::AbstractTerm,t2::Unfold.TimeExpandedTerm) = _symequal(t1,t2.term)
#function Effects._replace(matrix_term::MatrixTerm{<:Tuple{<:Unfold.TimeExpandedTerm}},typicals::Dict)

#    replaced_term = MatrixTerm((Effects._replace.(matrix_term.terms, Ref(typicals))...,))
#    basisfunctionTerm = getfield(matrix_term,:terms)[1]
#    return TimeExpandedTerm(replaced_term,basisfunctionTerm.basisfunction,basisfunctionTerm.eventfields)

#end
