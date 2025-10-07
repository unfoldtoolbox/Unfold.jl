#  # [Eventhandling](@id eventhandling)
# This tutorial introduces some helpful scripts to copy information between events - something that commonly happens.

# # Setup things
# Setup some packages

using Unfold
using DataFrames

# Let's start with a typical event structure you might get from a stimulus - response paradigm. The condition is only encoded in the stimulus, the reaction time only in the RT-event. Also, there is a nasty "break" event in between, making our task a bit harder
evts = DataFrame(
    :event => ["S", "R", "break", "S", "R", "S", "R"],
    :latency => [1, 3, 4, 6, 7, 10, 12],
    :condition => ["face", missing, missing, "bike", missing, "bike", missing],
    :rt => [missing, 0.3, missing, missing, 0.4, missing, 0.5],
)

# The quest is now to copy some info from `S` to `R` and others from `R` to preceeding `S`

evts_new = copy_eventinfo(evts, "S" => "R", :condition; search_fun = :forward)

# In order to copy the RT, we want to have a "lookback", that is the preceeding "S" shuold be used
copy_eventinfo!(evts_new, "R" => "S", "rt"; search_fun = :backward)


# Other convenient pattern is to copy to a new column - useful to test the copying behavior
copy_eventinfo!(evts_new, "R" => "S", "rt" => "newcolumn"; search_fun = :backward)


# You can also use "standard" DataFrames patterns e.g.

for grp in groupby(evts_new, :event)
    grp.newcolumn .= grp.event[1] .== "S" ? "STIMULUS" : "OTHER"
end
evts_new
