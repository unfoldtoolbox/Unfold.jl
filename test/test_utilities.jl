using CSV
using DelimitedFiles

function loadtestdata(testCase::String, dataPath::String=(@__DIR__)*"/data_new_testcases")
    #println(pwd()) # to debug github action
    data = readdlm(joinpath(dataPath, "$(testCase)_data.csv"), ',', Float64, '\n')
    data = dropdims(data,dims=1) # convert to vector
    evts = CSV.read(joinpath(dataPath, "$(testCase)_events.csv"), DataFrame)
    return data, evts
end


function gridexpand(conditionA  = [0,1.],continuousA = [-1.,0,1.])

tmp = reshape([ [x,y]  for x=conditionA, y=continuousA ],length(conditionA)*length(continuousA))
evts_grid = DataFrame(hcat(tmp...)')
rename!(evts_grid,["conditionA","continuousA"])
return evts_grid
end