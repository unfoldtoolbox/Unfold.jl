sub = 1

# load data
filename = "C:/Users/behinger/Downloads/p3_clean/$(sub)_P3_shifted_ds_reref_ucbip_hpfilt_ica_weighted_clean.set"
data,srate,evts_df,chanlocs_df,EEG = import_eeglab(filename)
data = data[1:29,:] # remove EOGs
# reref to average ref
data = data .- mean(data,dims=1)
# clean data
data = clean_data(data,EEG["reject"]["rejmanual"]) # we marked bad segments in matlab/eeglab
evts= parse_trigger_p3(evts_df)
evts = @where(evts, .!:invalidresponse)
# Mass Univariate
f_stim = @formula 0~0+trialtype
f_butt = @formula 0~0+trialtype#+answer
contrast = Dict(:trialtype => DummyCoding(base="distractor"),:answer=>DummyCoding(base="distractor"))
evts_stim = filter(x->x.eventtype=="stimulus",evts)
data_e,times = unfold.epoch(data=data,tbl=evts_stim,τ=(-0.8,1.3),sfreq=srate)
##
Xstim = designmatrix(UnfoldLinearModel,f_stim,filter(x->(x.eventtype=="stimulus"),evts))

um,res_tmp = fit(UnfoldLinearModel,f_stim, evts_stim,data_e,times,solver=(a,b)->unfold.solver_b2b(a,b))


#

e,modelinfo = optimizer_b2b(Xstim.Xs[good_ix,:],data_m)