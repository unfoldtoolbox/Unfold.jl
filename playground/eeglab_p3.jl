using unfold, DataFrames
using StatsModels
using CSV

using Statistics
using DataFramesMeta

include("../test/debug_readEEGlab.jl")
include("dataset-CORE_helper.jl")

# Single Subject

sub = 2
filename = "C:/Users/behinger/Downloads/P3_clean/$(sub)_P3_shifted_ds_reref_ucbip_hpfilt_ica_weighted_clean.set"

data,srate,evts_df,chanlocs_df,EEG = import_eeglab(filename)
data = clean_data(data,EEG["reject"]["rejmanual"]) # we marked bad segments in matlab/eeglab
evts= parse_trigger_p3(evts_df)
## Epoch & Mass Univariate
evts_stim = filter(x->x.eventtype=="stimulus",evts)
data_e,times = unfold.epoch(data=data,tbl=evts_stim,τ=(-0.3,1),sfreq=srate)

##
f = @formula 0~0+trialtype
um,res = fit(UnfoldLinearModel,f, evts_stim,data_e,times)
res.estimate = Float64.(res.estimate)

# plot it in makie
lines(Data(@where(res,:channel.==findfirst(chanlocs_df.labels.=="Cz"))),Group(color=:term),:colname_basis,:estimate)

##

Xstim = designmatrix(UnfoldLinearModel,@formula(0~0+trialtype),filter(x->(x.eventtype=="stimulus"),evts),firbasis((-0.3,1),srate,"stimulus"))

#evts.rt_s = evts.rt .- mean(filter(!isnan,evts.rt))
#evts.rt_s = evts.rt_s ./ std(filter(!isnan,evts.rt_s))
Xbutt = designmatrix(UnfoldLinearModel,@formula(0~1+trialtype),filter(x->(x.eventtype=="button")&(typeof(x.target)==String),evts),firbasis((-0.5,0.4),srate,"button"))

#,contrasts=Dict(:rt_s=>EffectsCoding,:correct=>EffectsCoding))

@time m = unfoldfit(UnfoldLinearModel,Xstim+Xbutt,data)


res = condense_long(m)



lines(Data(@where(res,:channel.==findfirst(chanlocs_df.labels.=="Cz"))),Group(color=:term,linestyle=:basisname),:colname_basis,:estimate)


## All Subjects

resAll = DataFrame()


for sub in 1:40
    if sub == 17 || sub == 25 || sub == 30 || sub == 31
        continue
    end
    println("subject $sub")

    # load data
    filename = "C:/Users/behinger/Downloads/p3_clean/$(sub)_P3_shifted_ds_reref_ucbip_hpfilt_ica_weighted_clean.set"
    data,srate,evts_df,chanlocs_df,EEG = import_eeglab(filename)
    data = data[1:29,:] # remove EOGs
    # reref to average ref
    data = data .- mean(data,dims=1)
    # clean data
    data = clean_data(data,EEG["reject"]["rejmanual"]) # we marked bad segments in matlab/eeglab
    evts= parse_trigger_p3(evts_df)
    # save
    if 1 == 1
        CSV.write(filename*".csv",evts)
    end
    evts = @where(evts, .!:invalidresponse)
    # Mass Univariate
    f_stim = @formula 0~0+trialtype
    f_butt = @formula 0~0+trialtype#+answer
    contrast = Dict(:trialtype => DummyCoding(base="distractor"),:answer=>DummyCoding(base="distractor"))
    evts_stim = filter(x->x.eventtype=="stimulus",evts)
    data_e,times = unfold.epoch(data=data,tbl=evts_stim,τ=(-0.8,1.3),sfreq=srate)
    um,res = fit(UnfoldLinearModel,f_stim, evts_stim,data_e,times,contrasts=contrast)
    res.basisname .= "stimulus"
    
    evts_butt = @where(evts,:eventtype.=="button",.!isnan.(:correct))
    data_e,times = unfold.epoch(data=data,tbl=evts_butt,τ=(-0.5,0.4),sfreq=srate)
    um,res_tmp = fit(UnfoldLinearModel,f_butt, evts_butt,data_e,times,contrasts=contrast)
    res_tmp.basisname .= "button"

    append!(res,res_tmp)
    res.group .="no-dc"

    # Deconv
    Xstim = designmatrix(UnfoldLinearModel,f_stim,filter(x->(x.eventtype=="stimulus"),evts),firbasis((-0.8,1.3),srate,"stimulus"),contrasts=contrast)
    Xbutt = designmatrix(UnfoldLinearModel,f_butt, filter(x->(x.eventtype=="button")&(typeof(x.target)==String),evts),firbasis((-0.5,0.4),srate,"button"),contrasts=contrast)

    m = unfoldfit(UnfoldLinearModel,Xstim+Xbutt,data)
    res_dc = condense_long(m)
    res_dc.group .="dc"
    append!(res,res_dc)
    res[:,:subject] .= sub
    append!(resAll,res)
end

##




using Makie
using StatsMakie
using MakieLayout
sub = 1
resAll = dropmissing(resAll,:estimate)
resAll.estimate = Float64.(resAll.estimate)
resAll.term = String.(resAll.term)



x = @linq resAll |>
where(:basisname.=="stimulus", :channel .==findfirst(chanlocs_df.labels.=="Cz"))

outer_padding = 30
scene, layout = layoutscene(outer_padding, resolution = (1800, 1200))
    #backgroundcolor = RGBf0(0.98, 0.98, 0.98))


for group in ["no-dc","dc"]
    
scene, layout = layoutscene(outer_padding, resolution = (1800, 1200))
    
axlist = []
linelist = []
for (i,sub) in enumerate(unique(x.subject))
    # keep the middle 4 plots free for the GA
    if i>20
        i = i+2
    end
    if i>26
        i = i+2
    end

    # calc col/row position
    col = (i-1)%6+1
    row = Int(ceil(i/6))
    
    println("$i,$sub,$row,$col")
    ax = layout[row, col] = LAxis(scene)
    append!(axlist,[ax]) # for later linking
    d = @linq(x |>where(:subject.==sub,:group.==group)) 
    l = lines!(ax,Data(d),Group(color=:term),:colname_basis,:estimate)
    #lines!(ax,[0.,0.],[-20.,30.],linestyle=:dot)
    tightlimits!(ax)
    ylims!(ax,[-15,15])
    
    
    #text!(ax,"Sub: $sub",position=(1.1,10.),textsize=4)
    append!(linelist,[l])

end

# plot average over subjects
x_avg = @linq x|>
    where(:group.==group)|>
    groupby([:term,:group,:colname_basis,:channel,:basisname]) |>
    based_on(estimate = mean(:estimate),sd = std(:estimate))

ax = layout[4:5, 3:4] = LAxis(scene, title = "Grand Average P3 $group")
l = lines!(ax,Data(x_avg),Group(color=:term),:colname_basis,:estimate)

tightlimits!(ax)
ylims!(ax,[-15,15])

linkaxes!(axlist...)
hidedecorations!.(axlist[1:end])
#axlist[end].xticks=[-0.7,0,1.0]
#axlist.[:frames][:linecolor].=nothing
scene




Makie.save("$group.png",scene)
end
