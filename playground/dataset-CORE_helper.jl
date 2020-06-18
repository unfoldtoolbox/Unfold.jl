

function parse_trigger_p3(evts)
    lookup = Dict(202.0=>"right",201.0=>"left",1.0 =>"A",2.0=>"B",3.0=>"C",4.0=>"D",5.0=>"E")
    stimulus =  []
    target=[]
    correct =[]
     button = []
     eventtype = []
     rts = []
    lasttarget = NaN
    laststimulus = NaN
    t = evts.type
    for i = 1:length(t)
        k = t[i]
        try
            
            push!(button,lookup[k]) # will fail except for 202 & 201
            
            push!(target,lasttarget)
            
            push!(stimulus,laststimulus)
            try
             laststim = lookup[(t[i-1])%10]
             rt = evts.latency[i]-evts.latency[i-1]
             push!(correct,laststim==lasttarget)
             push!(rts,rt)
            catch
             push!(correct,NaN)
             push!(rts,NaN)
            end
            push!(eventtype,"button")
            
        catch e
            #print(e)
            laststimulus = lookup[k%10]
            push!(stimulus,laststimulus)

            lasttarget = lookup[floor(k/10)]
            push!(target,lasttarget)
            push!(button,NaN)
            push!(correct,NaN)
            push!(eventtype,"stimulus")
            push!(rts,NaN)
        end
    end

    evts = hcat(evts, DataFrame(stimulus=stimulus,target=target,button=button,correct=correct,eventtype=eventtype,rt=Float64.(rts)))
    evts.stimtype = ["distractor","target"][(evts.stimulus .== evts.target).+1]
    return evts
end

