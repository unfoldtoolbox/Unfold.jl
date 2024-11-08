# helper scripts for event handlung
"""
copy field-info from source to closest target

`match_fun` can be "closest", "s before t" (eq. "t after s") or "t before s" (eq. "s after t")

## Example
julia> # copy reaction time values from button press to closest stimulus immediately before button press
julia> copy_eventinfo!(evts,"button"=>"stimulus",:reaction_time;match_fun="s after t")
"""
function copy_eventinfo!(evts, fromTo, field; match_fun = "closest")

    source = fromTo.first
    target = fromTo.second

    if isa(field, Pair)
        source_field = field.first
        target_field = field.second
    else
        source_field = field
        target_field = field
    end
    source_ix = findall(isequal(source), evts.trial_type)
    target_ix = findall(isequal(target), evts.trial_type)

    isempty(source_ix) && error("couldnt find source entries ($source) in evts.trial_type")
    isempty(target_ix) && error("couldnt find target entries ($target) in evts.trial_type")

    # for some matching functions we want to find the minimum, but no negative number
    filter_greaterzero = x -> x >= 0 ? x : Inf

    if match_fun == "closest"
        match_fun = (a, b) -> argmin(abs.(a .- b)) # find closest before or after
    elseif match_fun == "s after t" || match_fun == "t before s"
        match_fun = (s, t) -> argmin(filter_greaterzero.(s .- t)) # source after target
    elseif match_fun == "s before t" || match_fun == "t after s"
        match_fun = (s, t) -> argmin(filter_greaterzero.(t .- s)) # source before target
    elseif !callable(match_fun)
        error("bad match_fun input")
    end
    ix = [match_fun(evts[s, :latency], evts[target_ix, :latency]) for s in source_ix]
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
