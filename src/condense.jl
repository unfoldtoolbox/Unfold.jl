
times = [parse(Float64,k[2]) for k in split.(string.(names(d_tmp)[4:end]),'_')]
pred = [k[1] for k in split.(string.(names(d_tmp)[4:end]),'_')]
#times = parse.(Float64,string.(names(d_tmp)[3:end]))
p_df = DataFrame(times = times,beta = fixef(mm)[2:end],cond = pred,se=stderror(mm)[2:end])

vc = VarCorr(mm).σρ

ranef_sub =[k for k in vc[1].σ][2:end]
ranef_stim = []
groups = [v[1] for v in pairs(vc)]
content = [v[2].σ for v in pairs(vc)]

σρ = vc
nmvec = string.([keys(σρ)...])
cnmvec = string.(foldl(vcat, [keys(sig)...] for sig in getproperty.(values(σρ), :σ)))
σvec = vcat(collect.(values.(getproperty.(values(σρ), :σ)))...)
nvec = length.(keys.(getproperty.(values(vc),:σ)))

nmvec2 = []
for n in zip(nmvec,nvec)
    append!(nmvec2,repeat([n[1]],n[2]))
end

p_ranef = DataFrame(grouping = [k[1] for k in split.(nmvec2,'_')],σ=σvec,term=cnmvec)
p_ranef = p_ranef[2:end,:] # remove the lonely intercept
p_ranef.times = [parse(Float64,k[2]) for k in split.(p_ranef.term,'_')]
p_ranef.term = [k[1] for k in split.(p_ranef.term,'_')]

plot(p_ranef,x=:times,y=:σ,color=:term,linestyle=:grouping,xgroup=:grouping,Geom.subplot_grid(Geom.LineGeometry))

p_df.ymin = p_df.beta-p_df.se
p_df.ymax = p_df.beta+p_df.se
