# helper scripts for event handlung
"""
copy field-info from source to closest target

## Arguments
- `evts::DataFrame`
- `from_to::Pair` specifies from which entry to which entry to copy, e.g. "source"=>"target"
- `field::String` name of the column that contains the data to be copied to other events. Can be a `pair` in order to copy to a new column and thereby not replace any entries
## Keyword arguments
- `search_fun::Symbol/Function` can be `:closest` (default), `:forward` or `:backward`` or a custom function following the interface `search_fun(source_latency::Float64,target_latencies::Vector)` returning a single integer index
- `column::String` (default `"event"`) the column where the `from_to` source and target events can be found in

## Example
Copy reaction time values from button press to closest stimulus immediately before button press
julia> copy_eventinfo!(evts,"button"=>"stimulus","reaction_time";search_fun="s after t")
"""
copy_eventinfo(evts, args...; kwargs...) =
    copy_eventinfo!(deepcopy(evts), args...; kwargs...)

function copy_eventinfo!(
    evts,
    from_to::Pair,
    field;
    search_fun = "closest",
    column = "event",
)

    source = first(from_to)
    target = last(from_to)

    if isa(field, Pair)
        source_field = first(field) |> string
        target_field = last(field) |> string
    else
        source_field = field |> string
        target_field = field |> string
    end
    source_ix = findall(isequal(source), evts[:, column])
    target_ix = findall(isequal(target), evts[:, column])

    isempty(source_ix) && error("couldnt find source entries ($source) in evts.$column")
    isempty(target_ix) && error("couldnt find target entries ($target) in evts.$column")

    # for some matching functions we want to find the minimum, but no negative number
    filter_greaterzero = x -> x >= 0 ? x : Inf

    if search_fun == :closest
        search_fun = (a, b) -> argmin(abs.(a .- b)) # find closest before or after
    elseif search_fun == :backward
        search_fun = (s, t) -> argmin(filter_greaterzero.(s .- t)) # source after target
    elseif search_fun == :forward
        search_fun = (s, t) -> argmin(filter_greaterzero.(t .- s)) # source before target
    elseif !isa(search_fun, Function)
        error("bad search_fun input")
    end
    ix = [search_fun(evts[s, :latency], evts[target_ix, :latency]) for s in source_ix]
    ix = target_ix[ix]
    payload = evts[source_ix, source_field]

    # if targetfield does not exist, create one with Union(missing,target_field_datatype)
    if !any(names(evts) .== target_field)
        evts[!, target_field] =
            Array{Union{Missing,typeof(payload[1])}}(missing, nrow(evts))
    end

    # fill it in
    evts[ix, target_field] .= payload
    return evts
end
