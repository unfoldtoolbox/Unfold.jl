### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 65edde2e-d03e-11eb-300b-897bdc66d875
begin
import Pkg;
	Pkg.activate(".")
	Pkg.add(url="https://github.com/unfoldtoolbox/Unfold.jl")
	Pkg.add("PlutoUI")
	Pkg.add("StatsPlots")
	Pkg.add("DataFrames")
	Pkg.add("StatsModels")
end

# ╔═╡ 0a821c04-1aff-4fca-aab8-07a4d450838c
let
	using Revise
	using Unfold
	using StatsModels,PlutoUI,StatsPlots,DataFrames
end

# ╔═╡ bec07615-8ca3-4251-9e44-bf5044828274
include("../test/test_utilities.jl");

# ╔═╡ 72bffaed-d3d9-44d8-9a74-0e0a0d3b0882
md"""# Simple Unfold.jl Example with overlap correction"""

# ╔═╡ ee305b1e-0c62-4a9f-ab80-d937b1401347
md""" ### Setup"""

# ╔═╡ 2ec25bc8-1154-4850-ba02-62ed634cb558
md"""### Loading Data & Setting up Unfold"""

# ╔═╡ 76afd758-ac50-4f3c-a484-7f268ab059ef
data,evts = loadtestdata("test_case_3b");

# ╔═╡ b31a80a5-cca2-4c69-82f5-2b477f5b62be
begin
	plot(data[1:600],title="Simulated Data with overlap")
	vline!(evts.latency[1:20],legend=false)
end

# ╔═╡ 3f4dd059-e275-483d-af42-a09a2795fca4
data_e,times = Unfold.epoch(data=data,tbl=evts,τ=(-1.,1.9),sfreq=10);

# ╔═╡ c7fb216e-d177-4623-a135-db4cd856d4ad
f = @formula 0~1+conditionA+continuousA

# ╔═╡ f99efde3-bd34-46e2-ad1c-1a97ad866b79
fir = firbasis(τ=(-0.5,1.2),sfreq=20,name="basisA");

# ╔═╡ dcbdef13-675e-47f7-b484-d340f0e55934
f_dict = Dict("eventA"=>(f,fir)) # combine eventname, formula & basis

# ╔═╡ a55058c9-3b56-447b-8f1a-f7319a29f478
md"""### Fitting the models"""

# ╔═╡ bce42119-33d7-4e10-ac5e-b37028381b09
# overlap-corrected
m,res = fit(UnfoldLinearModel,f_dict,evts,data,eventcolumn="type");

# ╔═╡ 9aa42667-2629-48e1-b7fe-979dabc28f96
# non-overlapping (mass univariate)
m_e,res_e = fit(UnfoldLinearModel,f,evts,data_e,times);

# ╔═╡ 8bfb1b5c-1c6a-472c-9087-dd03b0884922
md"""### Plotting the results"""

# ╔═╡ b2561baf-c4c0-4e44-b184-16de912aa2ab
@df res_e plot(:colname_basis,:estimate,group=:term,title="Mass-Univariate (no Overlap Correction")

# ╔═╡ b488d24e-fd3b-4f71-91a4-97ce55ae2c70
@df res plot(:colname_basis,:estimate,group=:term,title="Overlap Corrected")

# ╔═╡ Cell order:
# ╠═72bffaed-d3d9-44d8-9a74-0e0a0d3b0882
# ╟─ee305b1e-0c62-4a9f-ab80-d937b1401347
# ╠═65edde2e-d03e-11eb-300b-897bdc66d875
# ╠═0a821c04-1aff-4fca-aab8-07a4d450838c
# ╠═bec07615-8ca3-4251-9e44-bf5044828274
# ╟─2ec25bc8-1154-4850-ba02-62ed634cb558
# ╠═76afd758-ac50-4f3c-a484-7f268ab059ef
# ╠═b31a80a5-cca2-4c69-82f5-2b477f5b62be
# ╠═3f4dd059-e275-483d-af42-a09a2795fca4
# ╠═c7fb216e-d177-4623-a135-db4cd856d4ad
# ╠═f99efde3-bd34-46e2-ad1c-1a97ad866b79
# ╠═dcbdef13-675e-47f7-b484-d340f0e55934
# ╟─a55058c9-3b56-447b-8f1a-f7319a29f478
# ╠═bce42119-33d7-4e10-ac5e-b37028381b09
# ╠═9aa42667-2629-48e1-b7fe-979dabc28f96
# ╟─8bfb1b5c-1c6a-472c-9087-dd03b0884922
# ╠═b2561baf-c4c0-4e44-b184-16de912aa2ab
# ╠═b488d24e-fd3b-4f71-91a4-97ce55ae2c70
