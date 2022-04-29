### A Pluto.jl notebook ###
# v0.19.2

using Markdown
using InteractiveUtils

# ╔═╡ 25d21156-d006-11eb-3655-61be76e7db93
begin
import Pkg;
	Pkg.activate("../../../")
end

# ╔═╡ 34aefa9a-e034-4293-b53d-cc1dc349ecc7
let
	using Revise
	using Unfold
	using DataFramesMeta
end

# ╔═╡ c73bf3ed-4af3-4e5a-997e-6399c9b85bbe
using StatsModels,PlutoUI,DataFrames,StatsPlots

# ╔═╡ d1e1de70-7e37-46a4-8ce2-3ca591bb84cf
include(pathof(Unfold)*"../../../test/test_utilities.jl")

# ╔═╡ 2e9ec6e0-03c6-41f6-bd3b-a44eaed5b553
data,evts = loadtestdata("test_case_4a");

# ╔═╡ 0b27e0bb-7aee-4672-8501-6096539c7ad7
evts.continuousA = rand(size(evts,1))

# ╔═╡ 97443eb8-ac67-4f8e-9e8a-eb09176c8516


# ╔═╡ 2a5d5493-b88d-473f-8d81-9e9240b3ffe0
#f  = @formula 0~1+conditionA+continuousA # 1
f  = @formula 0~1+continuousA # 1

# ╔═╡ 5d23b74c-988d-48bc-9c30-c2af0dbabdb4


# ╔═╡ 3f2b7a92-a064-4be6-8a39-cabaa2453777
basisfunction = Unfold.splinebasis(τ=(-1,1),sfreq=20,nsplines=10,name="basisA")

# ╔═╡ 8f32dbf3-8958-4e67-8e25-5f24d48058e1
Unfold.splinebasis

# ╔═╡ 7ef24e5a-63df-4e99-a021-7119d7e8d3bf
fir = firbasis(τ=(-1,1),sfreq=20,name="basisB")

# ╔═╡ 12c617c1-1994-4efe-b8bc-23fd2325c2ca
fir.kernel(2)

# ╔═╡ 6f7bc2f6-28de-40a5-83ee-c15cd6f74f6c
basisfunction.kernel(1)

# ╔═╡ e3d7ac15-e7a2-490b-abbc-aea3d7b4f4c7
size(Unfold.times(basisfunction))

# ╔═╡ 73022b4d-4edd-44e2-8cfa-26ded44a8921
size(Unfold.colnames(basisfunction))

# ╔═╡ d38f9cc9-409f-4c52-b0d2-d09312cc695a
size(basisfunction.colnames)

# ╔═╡ 786d0114-e10e-4e29-b6d7-bbe635c55f08
evts

# ╔═╡ c4227637-81b4-4ea7-a7b7-741860fef8c3
m = fit(UnfoldModel,Dict("eventA"=>(f,fir),"eventB"=>(f,basisfunction)),evts,data,eventcolumn="type");

# ╔═╡ cd19784c-686e-41c5-bd32-cdf793cecf03
res = coeftable(m)

# ╔═╡ 29e875d3-45b9-4a6a-aaa6-c6eaf35244c5
@df res plot(:time,:estimate,group=(:basisname,:coefname))

# ╔═╡ b669f55b-9121-4d66-9c82-f89bbfe0b2df
Unfold.coef(m)

# ╔═╡ 64cf7e78-2677-4551-9491-82c9f560ed61
	res

# ╔═╡ 02b8ff44-c514-492e-bf94-0cbf44431097
x = [1]

# ╔═╡ d964a5ca-7e63-41d5-bf97-526772d259e9
#[kernel(formTerm.basisfunction)(0) for formTerm in getfield.(designmatrix(m).formulas,:rhs)]

# ╔═╡ 1ee404e1-3039-410d-99de-ce713e751280
let
	basisnames = Unfold.get_basis_name(m)
	betaOut = []	
	for (k,b) = enumerate(unique(basisnames))
		
		ix = basisnames .== b
		
		betaOut = vcat(betaOut,designmatrix(m).formulas[k].rhs.basisfunction.kernel(0)	* m.modelfit.estimate[:,ix]')
	end
	plot(Float64.(betaOut))
end

# ╔═╡ 5466ec0b-d28a-4be7-b686-cb5ad7fea92b
let
	basisnames = Unfold.get_basis_name(m)
	k  = 2
	b = unique(basisnames)[k]
	

	bf = designmatrix(m).formulas[k].rhs.basisfunction
	ix = basisnames .== b
	
	hcat(bf.kernel.([0,0])	...) * m.modelfit.estimate[:,ix]'
end

# ╔═╡ 033cf2e0-a3ca-40e3-8f27-e88e017ca8de
size(Unfold.times(basisfunction))

# ╔═╡ a8238bd7-0198-45d2-89bd-d1434c537b2a
	designmatrix(m).formulas[2].rhs.basisfunction.kernel(0)


# ╔═╡ fe5a826f-aede-4e91-b6cd-ee7fca33ce31
zeros()

# ╔═╡ 8b9d8e32-d2ca-49db-986f-851ffabc3a3e
function undoBasis(d,m)
	bname = unique(d.basisname) # the current coefficient basename
	@assert(length(bname)==1)
	@assert(length(unique(d.channel))==1)
	bnames = unique(Unfold.get_basis_name(m)) # all model basenames
	
	# find out which formula / basisfunction we have
	k = findfirst(bname .== bnames) 
	bf = designmatrix(m).formulas[k].rhs.basisfunction
		
	
	estimate= bf.kernel.(0) *d.estimate
	@show size(bf.kernel.(0))
	@show size(Unfold.times(bf))
	return DataFrame(:time=>Unfold.times(bf),:estimate=>estimate)

end

# ╔═╡ 350a7f4b-f019-4bf6-8878-2cbb0cf82005
begin
	gd = groupby(res,[:basisname,:channel,:group,:coefname])
	a = combine(gd) do d
		return undoBasis(d,m)
	end
	a
end

# ╔═╡ a876f554-a034-461e-9522-490ca4bab5f8
@df a plot(:time,:estimate,group=(:basisname,:coefname))

# ╔═╡ 1d5b8b05-cdc6-4058-9bc1-bd52e6e3b444
Unfold.get_basis_name(m)

# ╔═╡ 111189ec-be04-4687-91ff-53b540a7b4db
Unfold.get_basis_name(m)

# ╔═╡ 4b572283-14cd-4a75-9ff2-f5a18f40228d
res

# ╔═╡ Cell order:
# ╠═25d21156-d006-11eb-3655-61be76e7db93
# ╠═34aefa9a-e034-4293-b53d-cc1dc349ecc7
# ╠═c73bf3ed-4af3-4e5a-997e-6399c9b85bbe
# ╠═d1e1de70-7e37-46a4-8ce2-3ca591bb84cf
# ╠═2e9ec6e0-03c6-41f6-bd3b-a44eaed5b553
# ╠═0b27e0bb-7aee-4672-8501-6096539c7ad7
# ╠═97443eb8-ac67-4f8e-9e8a-eb09176c8516
# ╠═2a5d5493-b88d-473f-8d81-9e9240b3ffe0
# ╠═5d23b74c-988d-48bc-9c30-c2af0dbabdb4
# ╠═3f2b7a92-a064-4be6-8a39-cabaa2453777
# ╠═8f32dbf3-8958-4e67-8e25-5f24d48058e1
# ╠═7ef24e5a-63df-4e99-a021-7119d7e8d3bf
# ╠═12c617c1-1994-4efe-b8bc-23fd2325c2ca
# ╠═6f7bc2f6-28de-40a5-83ee-c15cd6f74f6c
# ╠═e3d7ac15-e7a2-490b-abbc-aea3d7b4f4c7
# ╠═73022b4d-4edd-44e2-8cfa-26ded44a8921
# ╠═d38f9cc9-409f-4c52-b0d2-d09312cc695a
# ╠═786d0114-e10e-4e29-b6d7-bbe635c55f08
# ╠═c4227637-81b4-4ea7-a7b7-741860fef8c3
# ╠═cd19784c-686e-41c5-bd32-cdf793cecf03
# ╠═29e875d3-45b9-4a6a-aaa6-c6eaf35244c5
# ╠═b669f55b-9121-4d66-9c82-f89bbfe0b2df
# ╠═64cf7e78-2677-4551-9491-82c9f560ed61
# ╠═02b8ff44-c514-492e-bf94-0cbf44431097
# ╠═d964a5ca-7e63-41d5-bf97-526772d259e9
# ╠═1ee404e1-3039-410d-99de-ce713e751280
# ╠═5466ec0b-d28a-4be7-b686-cb5ad7fea92b
# ╠═033cf2e0-a3ca-40e3-8f27-e88e017ca8de
# ╠═a8238bd7-0198-45d2-89bd-d1434c537b2a
# ╠═350a7f4b-f019-4bf6-8878-2cbb0cf82005
# ╠═a876f554-a034-461e-9522-490ca4bab5f8
# ╠═fe5a826f-aede-4e91-b6cd-ee7fca33ce31
# ╠═8b9d8e32-d2ca-49db-986f-851ffabc3a3e
# ╠═1d5b8b05-cdc6-4058-9bc1-bd52e6e3b444
# ╠═111189ec-be04-4687-91ff-53b540a7b4db
# ╠═4b572283-14cd-4a75-9ff2-f5a18f40228d
