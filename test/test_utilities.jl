using CSV
function loadtestdata(testCase::String)
    data = CSV.read("./test/data/"*testCase*"_data.csv", header=0)
    data = convert(Matrix,data)
    data = dropdims(data,dims=1) # convert to vector
    evts = CSV.read("./test/data/"*testCase*"_events.csv");
    return data,evts
end
