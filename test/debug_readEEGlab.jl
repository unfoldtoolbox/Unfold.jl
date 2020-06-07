using unfold
using MAT
file = matopen("../data/sub-01_desc-MSDARKpreprocessed.set")
EEG = read(file, "EEG") # note that this does NOT introduce a variable ``varname`` into scope

srate = EEG["srate"]
data = EEG["data"]
evtsMatlab = EEG["event"]

DataFrame(duration = evtsMatlab["duration"],latency = evtsMatlab["latency"])
