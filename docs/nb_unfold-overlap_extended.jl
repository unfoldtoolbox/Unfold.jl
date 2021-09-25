### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 65edde2e-d03e-11eb-300b-897bdc66d875
begin
import Pkg;
	Pkg.activate(".")
	Pkg.add(url="https://github.com/unfoldtoolbox/Unfold.jl")
	Pkg.add("PlutoUI")
	Pkg.add("StatsPlots")
	Pkg.add("DataFrames")
	Pkg.add("StatsModels")
	Pkg.add("SignalAnalysis")
	Pkg.add("DSP")
end

# ╔═╡ 0a821c04-1aff-4fca-aab8-07a4d450838c
let
	using Revise
	using Unfold
	using StatsModels,PlutoUI,StatsPlots,DataFrames,Random
	using SignalAnalysis,DSP
end

# ╔═╡ bec07615-8ca3-4251-9e44-bf5044828274
include("../test/test_utilities.jl");

# ╔═╡ 72bffaed-d3d9-44d8-9a74-0e0a0d3b0882
md"""# Simple Unfold.jl Example with overlap correction"""

# ╔═╡ ee305b1e-0c62-4a9f-ab80-d937b1401347
md""" ### Setup
Setup might take a while because the libraries are very large. It might be better to install some of the packages to your julia-environment prior to starting Pluto
"""

# ╔═╡ 2ec25bc8-1154-4850-ba02-62ed634cb558
md"""### Loading Data & Setting up Unfold"""

# ╔═╡ 76afd758-ac50-4f3c-a484-7f268ab059ef
data,evts = loadtestdata("test_case_3b");

# ╔═╡ c7fb216e-d177-4623-a135-db4cd856d4ad
f = @formula 0~1+conditionA+continuousA

# ╔═╡ a55058c9-3b56-447b-8f1a-f7319a29f478
md"""### Fitting the models"""

# ╔═╡ 8bfb1b5c-1c6a-472c-9087-dd03b0884922
md"""### Plotting the results"""

# ╔═╡ 9a4c715a-7a51-411d-a354-f79ec2369cb7
let	
	md"""change window size τ = (-0.3, $(@bind τ2 Slider(0:0.1:5,default=2, show_value=true))) 
	"""
end

# ╔═╡ ba2f0c48-f5f8-40eb-a45b-441b8134ffb1
τ = (-0.3,τ2)

# ╔═╡ f99efde3-bd34-46e2-ad1c-1a97ad866b79
fir_basis = firbasis(τ=τ,sfreq=20,name="basisA");

# ╔═╡ dcbdef13-675e-47f7-b484-d340f0e55934
f_dict = Dict("eventA"=>(f,fir_basis)) # combine eventname, formula & basis

# ╔═╡ ed61f099-343c-47d4-9af4-c4486c93fd4a
let	
	md"""change noise: σ = $(@bind σ Slider(0:0.1:5,default=0, show_value=true))
	"""
end

# ╔═╡ 09705728-c6e5-41a5-9859-2c7e28828d5b
# add some noise
#data_noise = data .+ σ .* randn(size(data));
data_noise = data .+ σ .* rand(PinkGaussian(size(data,1)));

# ╔═╡ 1f3419f0-02e3-4346-8ac4-026cb76fa0fe
	md"""filter: $(@bind low_cutoff Slider(0.01:0.01:2,default=0.01, show_value=true))
	"""

# ╔═╡ 3e6f413d-f0b3-434f-a2bd-e002592247f8
hpf = fir(501,low_cutoff; fs=20);

# ╔═╡ b2893e0a-5712-4af2-923d-e80078f7277c
data_filter = filtfilt(hpf, data_noise);

# ╔═╡ b31a80a5-cca2-4c69-82f5-2b477f5b62be
begin
	plot(data_filter[1:400],title="Simulated Data with overlap")
	plot!(data_noise[1:400],title="Simulated Data with overlap")
	vline!(evts.latency[1:10],legend=false)
end

# ╔═╡ 3f4dd059-e275-483d-af42-a09a2795fca4
data_e,times = Unfold.epoch(data=data_filter,tbl=evts,τ=τ,sfreq=20);

# ╔═╡ 9aa42667-2629-48e1-b7fe-979dabc28f96
# non-overlapping (mass univariate)
m_e,res_e = fit(UnfoldLinearModel,f,evts,data_e,times);

# ╔═╡ bce42119-33d7-4e10-ac5e-b37028381b09
# overlap-corrected
m,res = fit(UnfoldLinearModel,f_dict,evts,data_filter,eventcolumn="type");

# ╔═╡ de267143-c029-4103-9fd7-ada153588ab6
begin
	
	
p1 = @df res_e plot(:colname_basis,:estimate,group=:coefname,title="Without ...",ylims=(-3,5),legend=false,xlims=(-0.3,5))
	p2 = @df res plot(:colname_basis,:estimate,group=:coefname,title="... with overlap correction",ylims=(-3,5),xlims=(-0.3,5),legend=:bottomright)
	plot(p1,p2,layout=(1,2))
end

# ╔═╡ Cell order:
# ╟─72bffaed-d3d9-44d8-9a74-0e0a0d3b0882
# ╟─ee305b1e-0c62-4a9f-ab80-d937b1401347
# ╠═65edde2e-d03e-11eb-300b-897bdc66d875
# ╠═0a821c04-1aff-4fca-aab8-07a4d450838c
# ╠═bec07615-8ca3-4251-9e44-bf5044828274
# ╟─2ec25bc8-1154-4850-ba02-62ed634cb558
# ╠═76afd758-ac50-4f3c-a484-7f268ab059ef
# ╠═09705728-c6e5-41a5-9859-2c7e28828d5b
# ╠═3e6f413d-f0b3-434f-a2bd-e002592247f8
# ╠═b2893e0a-5712-4af2-923d-e80078f7277c
# ╠═b31a80a5-cca2-4c69-82f5-2b477f5b62be
# ╠═3f4dd059-e275-483d-af42-a09a2795fca4
# ╠═c7fb216e-d177-4623-a135-db4cd856d4ad
# ╠═f99efde3-bd34-46e2-ad1c-1a97ad866b79
# ╠═dcbdef13-675e-47f7-b484-d340f0e55934
# ╟─a55058c9-3b56-447b-8f1a-f7319a29f478
# ╠═bce42119-33d7-4e10-ac5e-b37028381b09
# ╠═9aa42667-2629-48e1-b7fe-979dabc28f96
# ╟─8bfb1b5c-1c6a-472c-9087-dd03b0884922
# ╠═ba2f0c48-f5f8-40eb-a45b-441b8134ffb1
# ╟─9a4c715a-7a51-411d-a354-f79ec2369cb7
# ╟─ed61f099-343c-47d4-9af4-c4486c93fd4a
# ╟─1f3419f0-02e3-4346-8ac4-026cb76fa0fe
# ╠═de267143-c029-4103-9fd7-ada153588ab6
