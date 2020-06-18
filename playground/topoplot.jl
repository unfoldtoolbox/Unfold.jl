using Dierckx
using Makie
using ColorSchemes

# Load eeglab dataset
include("test/debug_readEEGlab.jl")
filename = "C:/Users/behinger/Downloads/1_P3_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist.set"
data,srate,evts_df,chanlocs_df = import_eeglab(filename)
##
X= Float64.(chanlocs_df.X[1:end-7]) # cast to Float, remove some irrelevant channels
Y= -Float64.(chanlocs_df.Y[1:end-7]) # Reverse Y, not sure why
D = data[1:end-7,10000] #choose random timepoint

# calculate 2D spline (Dierckx)
spl = Spline2D(Y,X,D,kx=3,ky=3,s=800) # s = smoothing parameter, kx/ky = spline order

# get extrema and extend by 20%
xlim = extrema(X) .+ abs.(extrema(X)).*[-.2, .2]
ylim = extrema(Y).+ abs.(extrema(Y)).*[-.2, .2]

# generate and evaluate a grid
xg = range(xlim[1],stop=xlim[2],  step=0.5)
yg = range(ylim[1],stop=ylim[2], step=0.5)
v = evalgrid(spl,yg,xg) # evaluate the spline at the grid locs

# remove everything in a circle (should be ellipse at some point)
ix = sqrt.([i.^2+j.^2 for i in yg, j in xg]).>90
v[ix] .= NaN
# nicer colormaps
cmap = ColorSchemes.vik;
# generate scene+subscene to keep aspect ratio
scene = Scene(resolution = (1000, 1000))
sub = Scene(scene, IRect(0, 0, 1000, 1000))


heatmap!(sub,yg,xg,v,colormap=cmap, show_axis = false)
contour!(sub,yg,xg,Float64.(v),linewidth=3,colormap=cmap)
scatter!(sub,Y,X,markersize=3,color="white",strokewidth=3,strokecolor="black") # add electrodes
[text!(sub,String(chanlocs_df.labels[i]),position=(Y[i],X[i]),textsize=7) for i in 1:length(X)]; # add labels
