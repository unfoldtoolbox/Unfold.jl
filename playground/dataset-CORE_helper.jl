

function parse_trigger_p3(evts)
    ##
    evts[:,:eventtype] .= ""
    evts[:,:stimulus]  .= ""
    evts[:,:target]    .= ""
    evts[:,:correct]   .= NaN
    evts[:,:answer]    .= ""
    evts[:,:rt]        .= NaN
    evts[:,:invalidresponse].= false
    
    evts.eventtype = ["stimulus","button"][(evts.type .> 100).+1]
    
    # fill in stimulus + target + answer
    for i = 1:size(evts,1)
        if evts[i,:eventtype] == "stimulus"
            evts[i,:stimulus] = string("ABCDEX"[Int(evts[i,:type]%10)])
            evts[i,:target]   = string("ABCDEX"[Int(floor(evts[i,:type]/10))])
        else evts[i,:eventtype] == "button"
            #201 => distractor
            #202 => Target
            #ANNA: I changed this to be %200 and adjusted the positon of dist & target
            evts[i,:answer] = ["distractor","target"][Int(evts[i,:type]%200)]
        end
    end
    
    # fill in stim+target for button
    
    for i = 2:size(evts,1)
        if evts[i,:eventtype] == "stimulus"
            continue
        end
        evts[i,:stimulus] =evts[i-1,:stimulus] 
        evts[i,:target]   =evts[i-1,:target] 
    end

    # trial type = target trial or distractor trial
    evts.trialtype = ["distractor","target"][(evts.stimulus .== evts.target).+1]

    #correct
    # ANNA: changes this, because I found it easier to understand
    evts.correct = (((evts.trialtype .== "target") .& (evts.answer .=="target")) .|
                    ((evts.trialtype .== "distractor") .& (evts.answer .=="distractor")))

    #evts.correct = (((evts.stimulus .== evts.target) .& (evts.answer .=="target")) .|
    #                ((evts.stimulus .!= evts.target) .& (evts.answer .=="distractor")))
    #evts[evts.eventtype.=="stimulus",:correct] .= NaN
    
    # invalidresponse
    if evts[1,:eventtype] == "button"
        evts[1,:invalidresponse] = true
    end
    # remove events that are double
    evts[vcat(true,evts[1:end-1,:eventtype] .== evts[2:end,:eventtype]),:invalidresponse] .= true
    evts[vcat(evts[1:end-1,:eventtype] .== evts[2:end,:eventtype],true),:invalidresponse] .= true
    
    #reaction time
    evts.rt = vcat(NaN,evts.latency[2:end]-evts.latency[1:end-1]).*1000/srate
    for i = 1:size(evts,1)-1
        if evts[i,:eventtype] == "button"
            continue
        end
        evts[i,:rt] =evts[i+1,:rt]        
    end


    
    return evts
end

