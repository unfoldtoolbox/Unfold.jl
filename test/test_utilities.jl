function loadtestdata(
    testCase::String;
    dataPath::String = (@__DIR__) * "/data_new_testcases",
)
    data = readdlm(joinpath(dataPath, "$(testCase)_data.csv"), ',', Float64, '\n')
    data = dropdims(data, dims = 1) # convert to vector
    evts = CSV.read(joinpath(dataPath, "$(testCase)_events.csv"), DataFrame)
    return data, evts
end

function gridexpand(conditionA = [0, 1.0], continuousA = [-1.0, 0, 1.0])
    tmp = reshape(
        [[x, y] for x in conditionA, y in continuousA],
        length(conditionA) * length(continuousA),
    )
    DataFrame(collect(hcat(tmp...)'), ["conditionA", "continuousA"])
end
