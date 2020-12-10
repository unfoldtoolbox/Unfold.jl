using CSV
using DelimitedFiles

function loadtestdata(testCase::String, dataPath::String=@__DIR__)
    #println(pwd()) # to debug github action
    data = readdlm(joinpath(dataPath, "data/$(testCase)_data.csv"), ',', Float64, '\n')
    data = dropdims(data,dims=1) # convert to vector
    evts = CSV.read(joinpath(dataPath, "data/$(testCase)_events.csv"), DataFrame)
    return data, evts
end
