

function extract_coef_info(coefs, ix)
    # 1 = eventname, 2 = coefname, 3 = colname
    return [c[ix] for c in split.(coefs, " : ")]
end

function get_coefnames(uf::Union{<:UnfoldModel,<:AbstractDesignMatrix})
    coefnames = Unfold.coefnames(formulas(uf))
    coefnames = vcat(coefnames...) # gets rid of an empty Any() XXX not sure where it comes from, only in MixedModel Timexpanded case
end



modelfit(uf::UnfoldModel) = uf.modelfit

StatsModels.coef(uf::UnfoldModel) = coef(modelfit(uf))
StatsModels.coef(mf::LinearModelFit) = mf.estimate


#=
function coeftable2(uf)
    coefs = coef(uf)
    ix = get_basis_indices.(Ref(uf), eventnames(uf))

    # XXX stderror, group, coefname
    coefs_views = [@view(coefs[:, i]) for i in ix]
    XXX = DataFrame() # specify the coefficient dataframe somehow?!
    result_to_table(uf, coefs_views, XXX)
end
=#

@traitfn function StatsModels.coeftable(uf::T) where {T<:UnfoldModel;ContinuousTimeTrait{T}}
    coefsRaw = get_coefnames(uf) |> poolArray
    coefs = extract_coef_info(coefsRaw, 2) |> poolArray
    #colnames_basis_raw = get_basis_colnames(formulas(uf))# this is unconverted basisfunction basis,
    colnames_basis = extract_coef_info(coefsRaw, 3) |> poolArray # this is converted to strings!
    basisnames = extract_coef_info(coefsRaw, 1) |> poolArray

    nchan = size(coef(uf), 1)

    coefs_rep = permutedims(repeat(coefs, 1, nchan), [2, 1])

    if length(colnames_basis) == 1
        colnames_basis = [colnames_basis]
    end


    colnames_basis_rep = permutedims(repeat(colnames_basis, 1, nchan), [2, 1])
    try
        colnames_basis_rep = parse.(Float64, colnames_basis_rep)
    catch
        #do nothing
    end
    chan_rep = repeat((1:nchan) |> poolArray, 1, size(colnames_basis_rep, 2))



    designkeys = collect(first.(design(uf)))
    if design(uf) != [:empty => ()]
        basiskeys = [string(b.name) for b in last.(last.(design(uf)))]

        eventnames =
            Array{Union{eltype(designkeys),eltype(basiskeys)}}(undef, length(basisnames))

        for (b, d) in zip(basiskeys, designkeys)

            eventnames[basisnames .== b] .= d
        end
    else
        @warn "No design found, falling back to basisnames instead of eventnames"
        eventnames = basisnames
    end

    eventnames_rep = permutedims(repeat(eventnames |> poolArray, 1, nchan), [2, 1])
    @debug typeof(coefs_rep) coefs_rep[1]
    return make_long_df(
        uf,
        coefs_rep,
        chan_rep,
        colnames_basis_rep,
        eventnames_rep,
        collabel(uf),
    )
end


@traitfn function StatsModels.coeftable(
    uf::T,
) where {T<:UnfoldModel;!ContinuousTimeTrait{T}}
    # Mass Univariate Case
    coefnames = get_coefnames(uf)

    colnames_basis = collect(last.(last.(design(uf)[1])))

    @debug "coefs: $(size(coefnames)),colnames_basis:$(size(colnames_basis)))"
    @debug "coefs: $coefnames,colnames_basis:$colnames_basis"
    nchan = size(coef(uf), 1)
    ncols = length(colnames_basis)
    ncoefs = length(coefnames)

    # replicate
    coefs_rep = permutedims(repeat(coefnames, outer = [1, nchan, ncols]), [2, 3, 1])
    colnames_basis_rep = permutedims(repeat(colnames_basis, 1, nchan, ncoefs), [2 1 3])
    chan_rep = repeat(1:nchan, 1, ncols, ncoefs)

    designkeys = collect(first.(design(uf)))
    if length(designkeys) == 1
        # in case of 1 event, repeat it by ncoefs
        eventnames = repeat([designkeys[1]], ncoefs)
    else
        eventnames = String[]
        for (ix, evt) in enumerate(designkeys)
            push!(eventnames, repeat([evt], size(modelmatrix(uf)[ix], 2))...)
        end
    end

    eventnames_rep = permutedims(repeat(eventnames, 1, nchan, ncols), [2, 3, 1])
    #
    results =
        make_long_df(uf, coefs_rep, chan_rep, colnames_basis_rep, eventnames_rep, :time)

    return results
end


#---
# Returns a long df given the already matched
function make_long_df(m, coefs, chans, colnames, eventnames, collabel)
    @assert all(size(coefs) .== size(chans)) "coefs, chans and colnames need to have the same size at this point, $(size(coefs)),$(size(chans)),$(size(colnames)), should be $(size(coef(m))))"
    @assert all(size(coefs) .== size(colnames)) "coefs, chans and colnames need to have the same size at this point"
    estimate, stderror, group = make_estimate(m)


    return DataFrame(
        Dict(
            :coefname => String.(linearize(coefs)) |> poolArray,
            :channel => linearize(chans),
            :eventname => linearize(eventnames),
            collabel => linearize(colnames),
            :estimate => linearize(estimate),
            :stderror => linearize(stderror),
            :group => linearize(group),
        ),
    )
end
#---------

stderror(m::UnfoldModel) = stderror(modelfit(m))

function stderror(m::LinearModelFit)
    if isempty(m.standarderror)
        stderror = fill(nothing, size(coef(m)))
    else
        stderror = Float64.(m.standarderror)
    end
end

function make_estimate(uf::UnfoldModel)
    return Float64.(coef(uf)), stderror(uf), fill(nothing, size(coef(uf))) |> poolArray
end

# Return the column names of the basis functions.
function get_basis_colnames(formula::FormulaTerm)
    return get_basis_colnames(formula.rhs)
end

function get_basis_colnames(rhs::Tuple)
    return colnames(rhs[1].basisfunction)
end

function get_basis_colnames(rhs::AbstractTerm)
    return colnames(rhs.basisfunction)
end

"""
    get_basisnames(model::UnfoldModel)

Return the basisnames for all predictor terms as a vector.

The returned vector contains the name of the event type/basis, repeated by their actual coefficient number (after StatsModels.apply_schema / timeexpansion).
If a model has more than one event type (e.g. stimulus and fixation), the vectors are concatenated.
"""
@traitfn function get_basis_names(m::T) where {T<:UnfoldModel;!ContinuousTimeTrait{T}}

    # Extract the event names from the design
    design_keys = first.((Unfold.design(m)))

    # Create a list of the basis names corresponding to each model term
    basisnames = String[]
    for (ix, event) in enumerate(design_keys)
        push!(basisnames, repeat([event], size(modelmatrix(m)[ix], 2))...)
    end
    return basisnames
end


@traitfn get_basis_names(m::T) where {T<:UnfoldModel;ContinuousTimeTrait{T}} =
    get_basis_names.(formulas(m))
function get_basis_names(m::FormulaTerm)
    bf = m.rhs.basisfunction
    #    @debug bf
    return repeat([name(bf)] |> poolArray, width(m.rhs))
end



"""
    get_basis_colnames(m)
    get_basis_colnames(formulas)
returns list of colnames - e.g. times for firbasis.

"""
get_basis_colnames(formulas::AbstractArray{<:FormulaTerm}) = get_basis_colnames.(formulas)



"""
    result_to_table(model<:UnfoldModel, eff::AbstractArray, events::Vector{<:DataFrame})
    result_to_table(
        eff::AbstractArray,
        events::Vector{<:DataFrame},
        times::Vector{<:Vector{<:Number}},
        eventnames::Vector)
    result_to_table(
        eff::Vector{<:AbstractArray},
        events::Vector{<:DataFrame},
        times::Vector,
        eventnames::Vector,
    )
Converts an array-result (prediction or coefficient) together with the events, to a tidy dataframe.

To support multi-event models, we expect everything to be put into `Vectors` - this should be refactored at some point to be compatible with broadcasting, but it is not right now.

## args
`eff`: A vector that contains the array(s) to be converted to a tidy dataframe. Each event in your events dataframe your have it's own array (i.e. the array should have `length(eff) == length(unique(events.event))`). Each array should be 3D, with channel x time x predictor
`events`: A vector of event-dataframes, each need to match `size(eff,3)`. Each individual event (`unique(events.event)`) should have it's own dataframe. E.g. `[[df_event1::DataFrame], [df_event1::DataFrame]]
`times`: A vector of time-vectors with `length(eff)`, each time-vector needs to match `size.(eff,2)`
`eventnames`: A vector of eventnames, either symbols or strings, should be a single entry per event
"""
result_to_table(model, eff, events::Vector{<:DataFrame}) =
    result_to_table(eff, events, times(model), eventnames(model))

# Array directly without in Vector
function result_to_table(
    eff::AbstractArray{<:Union{Number,Missing}},
    events::Vector,
    times::Vector,
    eventnames,
)
    1 == length(times) == length(events)
    result_to_table([eff], events, times, eventnames)
end
function result_to_table(
    eff::Vector{<:AbstractArray},
    events::Vector{<:DataFrame},
    times::Vector,
    eventnames::Vector,
)
    @assert length(eventnames) == length(events) == length(times)
    times_vecs = repeat.(poolArray.(times), size.(events, 1))
    # times_vecs = map((x, n) -> repeat(x, outer = n), poolArray.(times), size.(events, 1))
    # init the meta dataframe
    @debug typeof(times_vecs) size(times_vecs)
    data_list = []
    for (single_events, single_eff, single_times, single_eventname) in
        zip(events, eff, times_vecs, eventnames)
        @debug "sizes" size(single_events) size(single_eff) size(single_times)
        @debug "sizes2" size.(events) size.(times)
        ntimes = size(single_eff, 2)
        evts_tbl = repeat(Table(single_events), inner = (ntimes))
        time_event = Table(
            time = single_times |> poolArray,
            eventname = repeat([single_eventname] |> poolArray, length(single_times)),
        )
        @debug size(time_event) size(evts_tbl)
        @assert size(evts_tbl) == size(time_event)
        metadata = Table(evts_tbl, time_event)
        @debug "metadata" size(metadata)
        single_data = Table(
            yhat = single_eff[:],#vec(reshape(single_eff, :, 1)),
            channel = repeat(
                (1:size(single_eff, 1)) |> poolArray,
                outer = size(single_eff, 2) * size(single_eff, 3),
            ),
        )

        single_yhat = single_data.yhat
        single_channel = single_data.channel

        # get metadata columns as NamedTuple of vectors
        meta_cols = Tables.columntable(metadata)
        nchan = size(single_eff, 1)
        @debug size(metadata) typeof(metadata)
        # expand each metadata column by repeating each row for every channel
        expanded_meta = Dict{Symbol,Any}()
        for (k, v) in pairs(meta_cols)
            #@debug k, v
            key = isa(k, Symbol) ? k : Symbol(string(k))
            expanded_meta[key] = repeat(v, inner = nchan)
        end


        all_cols =
            merge(Dict(:yhat => single_yhat, :channel => single_channel), expanded_meta)

        push!(data_list, DataFrame(all_cols))

    end
    return reduce(vcat, data_list) |> t -> DataFrame(columns(t))


end
