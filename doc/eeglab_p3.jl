
#---
#title: "Overlap Correction with Linear Models (aka unfold.jl)"
#author: "Benedikt Ehinger, with help Dave Kleinschmidt"
#date: 2020-06-07
#options:
##    line_width: 92
#---

#```julia
using unfold, DataFrames
using StatsModels
using Plots,StatsPlots

#```

# Import data
#```julia
include("../test/debug_readEEGlab.jl")
filename = "C:/Users/behinger/Downloads/1_P3_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist.set"
data,srate,evts_df,chanlocs_df = import_eeglab(filename)
## ---


function parse_trigger_p3(evts,srate)
    lookup = Dict(202.0=>"right",201.0=>"left",1.0 =>"A",2.0=>"B",3.0=>"C",4.0=>"D",5.0=>"E")
    stimulus =  []
    target=[]
    correct =[]
     button = []
     eventtype = []
     rts = []
    lasttarget = NaN
    t = evts.type
    for i = 1:length(t)
        k = t[i]
        try
            
            push!(button,lookup[k]) # will fail except for 202 & 201
            
            push!(target,lasttarget)
            
            push!(stimulus,NaN)
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
            push!(stimulus,lookup[k%10])

            lasttarget = lookup[floor(k/10)]
            push!(target,lasttarget)
            push!(button,NaN)
            push!(correct,NaN)
            push!(eventtype,"stimulus")
            push!(rts,NaN)
        end
    end
    return stimulus,target,button,correct,eventtype,rts
end


stimulus,target,button,correct,eventtype,rts= parse_trigger_p3(evts_df,srate)
evts = hcat(evts_df, DataFrame(stimulus=stimulus,target=target,button=button,correct=correct,eventtype=eventtype,rt=rts))
evts.stimtype = ["distractor","target"][(evts.stimulus .== evts.target).+1]
##
#evts.eventtype = ifelse(isnan.(evts.stimulus),"button","stimulus")
#```
#```julia
evts_stim = filter(x->x.eventtype=="stimulus",evts)
data_e,times = unfold.epoch(data=data,tbl=evts_stim,τ=(-0.3,1),sfreq=srate)

##

f = @formula 0~0+stimtype
um,res = fit(UnfoldLinearModel,f, evts_stim,data_e,times)
@df filter(x->x.channel==findfirst(chanlocs.labels.=="Cz"),res) plot(:colnames_basis,:estimate,group=(:term))
#```
##

Xstim = designmatrix(UnfoldLinearModel,@formula(0~0+stimtype),filter(x->x.eventtype=="stimulus",evts),firbasis((-0.3,1),srate))

Xbutt = designmatrix(UnfoldLinearModel,@formula(0~0+correct+rt),filter(x->x.eventtype=="button",evts),firbasis((-0.5,0.3),srate))

m = unfoldfit(UnfoldLinearModel,Xstim+Xbutt,data)

res = condense_long(m)
res.times = parse.(Float64,res.colnames_basis)
res.event = [x[1] for x in split.(res.term,": ")]

@df filter(x->(x.event=="stimtype")&(x.term!="correct: NaN")&(x.channel==findfirst(chanlocs.labels.=="Cz")),res) plot(:times,:estimate,group=(:term))

@df filter(x->(x.event=="correct")&(x.term!="correct: NaN")&(x.channel==findfirst(chanlocs.labels.=="Cz")),res) plot(:times,:estimate,group=(:term))