using Makie
using MakieLayout
using Random




data = rand(10,50000) # 10 channels, 50000 samples
nchan = size(data,1)
nsamp = size(data,2)
srate = 1000
#x_zoom = slider(LinRange(0, 100, 100), raw = true, camera = campixel!, start = 10)
#y_zoom = slider(LinRange(0, 100, 100), raw = true, camera = campixel!, start = 10)
#x_pos = slider(LinRange(0, nsamp, 1000), raw = true, camera = campixel!, start = 0)

##

scene, layout = layoutscene(30, resolution = (1800, 1200));



ax     = layout[1, 1] = LAxis(scene);
x_zoom = layout[2, 1] = LSlider(scene, range = 1:0.1:100, startvalue = 100,tellheight=true);
x_pos  = layout[4, 1] = LSlider(scene, range = 0:1:nsamp, startvalue = 1,tellheight = true);
y_zoom = layout[:, 2] = LSlider(scene, range = 0:0.1:10, horizontal = false,
         tellwidth = true, height = nothing, width = Auto(),startvalue=3);
         
##

data_zoom = []
for ch = 1:size(data,1)
    ix_range = lift((from,delta)->range(Int(floor(from+1)),stop=Int(floor(1+from+delta))),x_pos.value,x_zoom.value)
    d = lift((ix,zoom)->ch .+ zoom.*data[ch,ix],ix_range,y_zoom.value)
    append!(data_zoom,[d])
    lines!(ax,@lift($ix_range.*1/srate),data_zoom[ch])
end

#curr_xlim =    xlims!(x_zoom.value,x_zoom.value:)
#ylims!([-y_zoom[end][:value].*y_zoom_fac,y_zoom[end][:value].*y_zoom_fac])

#p = heatmap()
scene
#final = hbox(p, vbox(x_pos,x_zoom,y_zoom))#, parent = Scene(resolution = (1500, 1500)))

