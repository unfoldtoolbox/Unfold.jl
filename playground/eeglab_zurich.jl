
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

include("../test/debug_readEEGlab.jl")
filename = "C:/Users/behinger/Downloads/Benedikt.set"
data,srate,evts_df,chanlocs_df = import_eeglab(filename)


Xstim = designmatrix(UnfoldLinearModel,@formula(0~1+cond),filter(x->(x.type=="stimulus_onset"),evts_df),firbasis((-0.3,1),srate,"stimulus"))

Xsacc = designmatrix(UnfoldLinearModel,@formula(0~1+cond+dir),filter(x->(x.type=="saccade"),evts_df),firbasis((-0.5,0.5),srate,"sacc"))

Xblink = designmatrix(UnfoldLinearModel,@formula(0~1),filter(x->(x.type=="L_blink"),evts_df),firbasis((-0.5,0.5),srate,"blink"))

##

@time m = unfoldfit(UnfoldLinearModel,Xstim+Xsacc+Xblink,data);
X = Xstim+Xsacc+Xblink
using IncompleteLU
X2 = ilu(X.Xs)
##
res = condense_long(m)
res.times = res.colname_basis



@df filter(x->(x.basisname=="stimulus")&(x.term!="correct: NaN")&(x.channel==findfirst(chanlocs_df.labels.=="E14")),res) plot(:times,:estimate,group=(:term))

@df filter(x->(x.basisname=="button")&(x.term!="correct: NaN")&(x.channel==findfirst(chanlocs_df.labels.=="Cz")),res) plot(:times,:estimate,group=(:term))