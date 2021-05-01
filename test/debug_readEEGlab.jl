#using Unfold
using MAT
using DataFrames
using DelimitedFiles

function import_eeglab(filename)

    file = MAT.matopen(filename)
    # open file
    EEG = read(file, "EEG")


    function parse_struct(s::Dict)
        return DataFrame(map(x->dropdims(x,dims=1),values(s)),collect(keys(s)))
    end

    evts_df = parse_struct(EEG["event"])
    chanlocs_df = parse_struct(EEG["chanlocs"])
    
#    epoch_df = parse_struct(EEG["epochs"])
    
    srate = EEG["srate"]
    if typeof(EEG["data"]) == String
        datapath = joinpath(splitdir(filename)[1],EEG["data"])
        data = Array{Float32, 3}(undef,Int(EEG["nbchan"]),Int(EEG["pnts"]),Int(EEG["trials"]))
        read!(datapath,data)
    else
        data = EEG["data"]
    end
    
    if (ndims(data)==3) & (size(data,3) == 1)
        data = dropdims(data,dims=3)
    end
return data,srate,evts_df,chanlocs_df,EEG
end