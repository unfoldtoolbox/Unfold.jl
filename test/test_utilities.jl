using CSV
function loadtestdata(testCase::String,dataPath::String="")
    #println(pwd()) # to debug github action
    data = CSV.read(dataPath*"data/"*testCase*"_data.csv", header=0)
    data = convert(Matrix,data)
    data = dropdims(data,dims=1) # convert to vector
    evts = CSV.read(dataPath*"data/"*testCase*"_events.csv");
    return data,evts
end
