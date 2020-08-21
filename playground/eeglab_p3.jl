using unfold, DataFrames
using StatsModels
using CSV

using Statistics
using DataFramesMeta
using RCall

using Makie
using StatsMakie
using Formatting
using MakieLayout

include("../test/debug_readEEGlab.jl")
include("dataset-CORE_helper.jl")

## All Subjects

resAll = DataFrame()
resEvts = DataFrame()
global chanlocs_df = [] # keep the loaded chanlocs

# Calculate Unfold with/without
for sub in 1:40
    global chanlocs_df
    # These subjects have some issues e.g. only target trials (at least one of them, not sure what the probem was with the others)
    if sub == 17 || sub == 25 || sub == 30 || sub == 31
        continue
    end
    ##
    
    println("subject $sub")

    # load data
    filename = "C:/Users/behinger/Downloads/p3_clean/$(sub)_P3_shifted_ds_reref_ucbip_hpfilt_ica_weighted_clean.set"
    data,srate,evts_df,chanlocs_df,EEG = import_eeglab(filename)
    data = data[1:29,:] # remove EOGs
    # reref to average ref
    data = data .- mean(data,dims=1)
    # clean data
    data = unfold.clean_data(data,EEG["reject"]["rejmanual"]) # we marked bad segments in matlab/eeglab
    evts= parse_trigger_p3(evts_df,srate)

    # save all evts for e.g. RT analysis
    evts[:,:subject] .= sub
    append!(resEvts,evts)

    #continue # for speed up dont fit the models ;)

    # save as csv?
    if 1 == 0
        CSV.write(filename*".csv",evts)
    end
    evts = @where(evts, .!:invalidresponse)

    # Mass Univariate
    f_stim = @formula 0~0+trialtype
    f_butt = @formula 0~0+trialtype#+answer
    contrast = Dict(:trialtype => DummyCoding(base="distractor"),:answer=>DummyCoding(base="distractor"))
    evts_stim = filter(x->x.eventtype=="stimulus",evts)
    data_e,times = unfold.epoch(data=data,tbl=evts_stim,τ=(-0.8,1.3),sfreq=srate)
    um,res = unfold.fit(UnfoldLinearModel,f_stim, evts_stim,data_e,times,contrasts=contrast)
    res.basisname .= "stimulus"
    
    evts_butt = @where(evts,:eventtype.=="button",.!isnan.(:correct))
    data_e,times = unfold.epoch(data=data,tbl=evts_butt,τ=(-0.5,0.4),sfreq=srate)
    um,res_tmp = unfold.fit(UnfoldLinearModel,f_butt, evts_butt,data_e,times,contrasts=contrast)
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


resAll = dropmissing(resAll,:estimate)
resAll.estimate = Float64.(resAll.estimate)
resAll.term = String.(resAll.term)

subjectwiseRT = @linq resEvts |>where(:invalidresponse .== 0) |> where(:eventtype .=="button") |> groupby([:subject]) |> based_on(m=mean(:rt),std=std(:rt))
datRT = @linq resEvts |> where(:invalidresponse .== 0,:eventtype .=="button") |> groupby(:subject)


## Calculate shifted lognormal on reaction times
@rlibrary brms

R"library('brms')"
brmfit = brm(R"rt  ~ answer", datRT[1], family=shifted_lognormal())
@rput(brmfit)


fitList = []

for g in datRT
    @rput(g)
    brmfit = R"posterior_summary(update(brmfit, newdata=g,recompile=FALSE))"
    append!(fitList,[brmfit])
end

lognorm_param = hcat([rcopy(f)[:,1] for f in fitList]...)
## Plot the effects per subject
d = @linq resAll |>
    where(:basisname.=="stimulus", :channel .==findfirst(chanlocs_df.labels.=="Cz"))
d = unstack(select(d,Not(:stderror)),:term,:estimate)
#d.diff = d.dc - d[:,"no-dc"]
d.subject = Float64.(d.subject)
d.diff = d[:,end] - d[:,end-1]

# Plot difference wave according to sort
axlist = []

scene,  layout = layoutscene(30, resolution = (1800, 1200))
#ix = sortperm(lognorm_param[1,:])
minrt = (@linq datRT|>based_on(min=minimum(:rt))).min
ix = sortperm(minrt)
for (i,sub) in enumerate(unique(d.subject)[ix])
    ncols = 5
    row = Int(ceil(i/ncols))
    col = Int((i%ncols)+1)
    println("$row,$col")
    #subtitle = "$sub,med:"*sprintf1("%.1f",exp(lognorm_param[1,ix[i]]))*", ndt:"*sprintf1("%.0f",(lognorm_param[4,ix[i]]))*" minrt:"*sprintf1("%.1f",minrt[ix[i]])
    subtitle = "$sub, minrt:"*sprintf1("%.1f",minrt[ix[i]])
    ax = layout[row, col] = LAxis(scene,title=subtitle)
    append!(axlist,[ax]) # for later linking

    d_sub = @linq d |>where(:subject.==sub)

    lines!(ax,Data(d_sub),Group(:group),:colname_basis,:diff)
end
linkaxes!(axlist...)
scene

##
x = @linq resAll |>
where(:basisname.=="stimulus", :channel .==findfirst(chanlocs_df.labels.=="Cz"))

outer_padding = 30
scene, layout = layoutscene(outer_padding, resolution = (1800, 1200))
    #backgroundcolor = RGBf0(0.98, 0.98, 0.98))


for group in ["no-dc","dc"]
#group = ["no-dc"]#
    
scene,  layout = layoutscene(outer_padding, resolution = (1800, 1200))
    
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
    ax = layout[row, col] = LAxis(scene,title="STD:"*sprintf1("%.1f",(@linq subjectwiseRT|> where(:subject .== sub)).std[1]))
    append!(axlist,[ax]) # for later linking
    d = @linq(x |>where(:subject.==sub,:group.==group)) 
    l = lines!(ax,Data(d),Group(color=:term),:colname_basis,:estimate)

    # text is absolutely fucked in makie. Don't even bother!

    #annotations!(ax, ["$((@linq subjectwiseRT|> where(:subject .== sub)).std)"],[Point(0.5,-13)],textsize=3) 
    #text!(ax, "$((@linq subjectwiseRT|> where(:subject .== sub)).std)",position=(0.5,-8),textsize=20px) 
    
    #lines!(ax,[0.,0.],[-20.,30.],linestyle=:dot)
    tightlimits!(ax)
    ylims!(ax,[-15,15])
    xlims!(ax,[-0.7,1.25])
    
    
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
#[ax.ticks]
#hidedecorations!.(axlist[1:end])

#axlist[end].xticks=[-0.7,0,1.0]
#axlist.[:frames][:linecolor].=nothing
scene   




Makie.save("$group.png",scene)
end
