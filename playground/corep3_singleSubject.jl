using unfold, DataFrames
using StatsModels
using CSV

using Statistics
using DataFramesMeta

include("../test/debug_readEEGlab.jl")
include("dataset-CORE_helper.jl")

## Single Subject

sub = 2
filename = "C:/Users/behinger/Downloads/P3_clean/$(sub)_P3_shifted_ds_reref_ucbip_hpfilt_ica_weighted_clean.set"

data,srate,evts_df,chanlocs_df,EEG = import_eeglab(filename)
data = unfold.clean_data(data,EEG["reject"]["rejmanual"]) # we marked bad segments in matlab/eeglab
evts= parse_trigger_p3(evts_df,srate)
evts = @where(evts, .!:invalidresponse)

## Epoch & Mass Univariate
evts_stim = filter(x->x.eventtype=="stimulus",evts)
data_e,times = unfold.epoch(data=data,tbl=evts_stim,τ=(-0.3,1),sfreq=srate)


##

using Splines2
include("unfold_splines2.jl")
f = @formula 0~0+trialtype+bs2(rt,5)
um,res = fit(UnfoldLinearModel,f, evts_stim,data_e,times)
res.estimate = Float64.(res.estimate)

# plot it in makie
lines(Data(@where(res,:channel.==findfirst(chanlocs_df.labels.=="Cz"))),Group(color=:term),:colname_basis,:estimate)


##

Xstim_rt = designmatrix(UnfoldLinearModel,@formula(0~0+trialtype+bs2(rt,8)),filter(x->(x.eventtype=="stimulus"),evts),firbasis((-0.3,1),srate,"stimulus"))
Xbutt_rt = designmatrix(UnfoldLinearModel,@formula(0~1+trialtype+bs2(rt,8)),filter(x->(x.eventtype=="button")&(typeof(x.target)==String),evts),firbasis((-0.5,0.4),srate,"button"))

Xstim = designmatrix(UnfoldLinearModel,@formula(0~0+trialtype),filter(x->(x.eventtype=="stimulus"),evts),firbasis((-0.3,1),srate,"stimulus"))
Xbutt = designmatrix(UnfoldLinearModel,@formula(0~1+trialtype),filter(x->(x.eventtype=="button")&(typeof(x.target)==String),evts),firbasis((-0.5,0.4),srate,"button"))
#,contrasts=Dict(:rt_s=>EffectsCoding,:correct=>EffectsCoding))

@time m_rt = unfoldfit(UnfoldLinearModel,Xstim_rt+Xbutt_rt,data)
@time m = unfoldfit(UnfoldLinearModel,Xstim+Xbutt,data)

##


##

x = expandgrid(Dict(:trialtype=>["target","distractor"],:rt=>[400,mean(m_rt.X.events[1].rt),500,600.]))
x = expandgrid(Dict(:trialtype=>["target","distractor"],:rt=>[mean(m_rt.X.events[1].rt)]))

yhat_rt = unfold.predict(m_rt,x)
yhat = unfold.predict(m,x)

res_j = join(yhat_rt,yhat,on=names(yhat)[names(yhat).!="yhat"],makeunique=true)

lines(Data(@where(res_j,:channel.==findfirst(chanlocs_df.labels.=="Cz"),:basisname .=="stimulus")),Group(color=:rt,linestyle=:trialtype),:times,:yhat)
lines!(Data(@where(res_j,:channel.==findfirst(chanlocs_df.labels.=="Cz"),:basisname .=="stimulus")),Group(color=:rt,linestyle=:trialtype),:times,:yhat_1)

lines(Data(@where(stack(res_j,[:yhat,:yhat_1]),:channel.==findfirst(chanlocs_df.labels.=="Cz"),:basisname .=="button")),Group(color=:variable,linestyle=:trialtype),:times,:value)
#lines!(Data(@where(res_j,:channel.==findfirst(chanlocs_df.labels.=="Cz"),:basisname .=="stimulus")),Group(color=:rt,linestyle=:trialtype),:times,:yhat_1)


##
import StatsBase
x = expandgrid(Dict(:trialtype=>["target","distractor"],:rt=>StatsBase.percentile(m_rt.X.events[1].rt,[25. 35. 50 65 75])))
yhat_rt = unfold.predict(m_rt,x)
sc = Scene()
a = lines!(sc,Data(@where(yhat_rt,:channel.==findfirst(chanlocs_df.labels.=="Cz"),:basisname .=="stimulus")),Group(color=:rt,linestyle=:trialtype),:times,:yhat)
lgd = legend(a.plots[2].plots,unique(kron(string.(x.rt) .*" x ", x.trialtype)))
p_rt = vbox(a,lgd)

yhat_rt = unfold.predict(m,x)
sc = Scene()
a = lines!(sc,Data(@where(yhat_rt,:channel.==findfirst(chanlocs_df.labels.=="Cz"),:basisname .=="stimulus")),Group(color=:rt,linestyle=:trialtype),:times,:yhat)
lgd = legend(a.plots[2].plots,unique(kron(string.(x.rt) .*" x ", x.trialtype)))
p = vbox(a,lgd)

hbox(p_rt,p)

##
res_rt = condense_long(m_rt)
res = condense_long(m)
res_j.estimate_1 = Float64.(res_j.estimate_1)