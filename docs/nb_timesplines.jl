### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 25d21156-d006-11eb-3655-61be76e7db93
begin
import Pkg;
	Pkg.activate("../")
end

# ╔═╡ 34aefa9a-e034-4293-b53d-cc1dc349ecc7
let
	using Revise
	using Unfold
end

# ╔═╡ c73bf3ed-4af3-4e5a-997e-6399c9b85bbe
using StatsModels,PlutoUI,StatsPlots,DataFrames,BlockDiagonals

# ╔═╡ d1e1de70-7e37-46a4-8ce2-3ca591bb84cf
include("../dev/Unfold/test/test_utilities.jl")

# ╔═╡ 2e9ec6e0-03c6-41f6-bd3b-a44eaed5b553
data,evts = loadtestdata("test_case_4a");

# ╔═╡ 97443eb8-ac67-4f8e-9e8a-eb09176c8516


# ╔═╡ 2a5d5493-b88d-473f-8d81-9e9240b3ffe0
#f  = @formula 0~1+conditionA+continuousA # 1
f  = @formula 0~1 # 1

# ╔═╡ 5d23b74c-988d-48bc-9c30-c2af0dbabdb4


# ╔═╡ 3f2b7a92-a064-4be6-8a39-cabaa2453777
basisfunction = Unfold.splinebasis((-1,1),20,10,"basisA")

# ╔═╡ 7ef24e5a-63df-4e99-a021-7119d7e8d3bf
fir = firbasis(τ=(-1,1),sfreq=20,name="basisB")

# ╔═╡ 12c617c1-1994-4efe-b8bc-23fd2325c2ca
fir.kernel(2)

# ╔═╡ 6f7bc2f6-28de-40a5-83ee-c15cd6f74f6c
basisfunction.kernel(1)

# ╔═╡ 786d0114-e10e-4e29-b6d7-bbe635c55f08
evts

# ╔═╡ c4227637-81b4-4ea7-a7b7-741860fef8c3
m,res = fit(UnfoldLinearModel,Dict("eventA"=>(f,fir),"eventB"=>(f,basisfunction)),evts,data,eventcolumn="type");

# ╔═╡ cd19784c-686e-41c5-bd32-cdf793cecf03
res

# ╔═╡ 29e875d3-45b9-4a6a-aaa6-c6eaf35244c5
@df res plot(:colname_basis,:estimate,group=:basisname)

# ╔═╡ 02b8ff44-c514-492e-bf94-0cbf44431097
x = [1]

# ╔═╡ d964a5ca-7e63-41d5-bf97-526772d259e9
[formTerm.basisfunction.kernel(0) for formTerm in getfield.(m.X.formulas,:rhs)]

# ╔═╡ 54862fe8-be53-47c1-a249-e47bcd492ce4
formTerm.basisfunction

# ╔═╡ 1ee404e1-3039-410d-99de-ce713e751280
let
	basisnames = Unfold.get_basis_name(m)
	betaOut = []	
	for (k,b) = enumerate(unique(basisnames))
		
		ix = basisnames .== b
		
		betaOut = vcat(betaOut,m.X.formulas[k].rhs.basisfunction.kernel(0)	* m.beta[:,ix]')
	end
	plot(Float64.(betaOut))
end

# ╔═╡ d27167fb-9cf4-4e0d-8c65-fde19577be0f
size(m.X.formulas[2].rhs.basisfunction.colnames)

# ╔═╡ 5efce43e-7cb5-4af4-84ef-4dc2c82090a7
res

# ╔═╡ cf1f63c7-8858-462a-809e-139d3872c989


# ╔═╡ c4ebfcbb-4230-46d2-854f-f7f7c3e63a98


# ╔═╡ be2758e1-83f6-4555-a0bf-2941687c3097


# ╔═╡ 07077e0e-87aa-437d-893e-9668e62f5d01


# ╔═╡ ce39efac-c362-4d9f-b6c4-38e13251e68f


# ╔═╡ 454d82fd-4c6d-4f08-a516-4bb020bad421


# ╔═╡ ba67232a-b1a1-4de3-9e38-dacb2ee81fa2


# ╔═╡ 826f0033-1add-4fc3-920a-7d86ef62d31e


# ╔═╡ 706be005-8330-4812-a596-6469f7a57201


# ╔═╡ 64735d5c-eff5-47cb-a843-6d0a941975c0


# ╔═╡ 6e6fa3e3-7c9b-4a32-a1a5-fb260f143f3d


# ╔═╡ 26f58723-c663-468d-b578-a59eff0d9acf


# ╔═╡ fa514e1d-e251-4431-97f3-780da7fdfb95


# ╔═╡ a009593b-f1f3-47ab-af5a-60612dfd37a9


# ╔═╡ b5e5f89b-158b-44b0-8945-5d2f0bd79ef7


# ╔═╡ 07e724d3-4180-471d-8a08-9a0f0b869093


# ╔═╡ ded80200-9a8c-47fe-9e10-c1a9e5ad08ac


# ╔═╡ 671696e3-68a0-45cd-a08d-aa6c7e400230


# ╔═╡ 98c4ffea-0228-4e2b-850e-92c64391bf61


# ╔═╡ 29ebdb59-57c2-4000-b037-2483e8c69d01


# ╔═╡ dc4e4671-a744-4003-adb5-010a3e67628a


# ╔═╡ 70d14446-b3f0-4bf3-8b2b-e5836bd46633
	

# ╔═╡ 282ac7f1-1647-447f-af4c-371251524a37


# ╔═╡ 66121b10-930f-49a7-a222-2d1e5ea81343
k = 1

# ╔═╡ 0fe93070-e434-43d8-ad95-521ef3d42402
m.X.formulas[k].rhs.basisfunction.kernel(0)	* m.beta[:,ix]'

# ╔═╡ 0aa995b2-1c7b-40b8-8913-b17852a290f9
Unfold.get_basis_name(m) * m.beta[:,ix]'

# ╔═╡ 357bc67d-e581-41a0-a378-905072e55289


# ╔═╡ d480d611-d53e-4f15-b2e9-a22ad993b64e


# ╔═╡ 9f4a59b6-a422-457a-a5d9-88beb65f7a78


# ╔═╡ 0181d139-7a38-4ba1-ac6b-0c3d76ab6a47


# ╔═╡ 1bfd0909-05b6-4782-a34d-1ea1ba649115
size(m.beta)

# ╔═╡ 52ec6011-ecf6-4b26-9bf5-0523f1ad035f
size(basisfunction.kernel(0))

# ╔═╡ d9a14d6d-5edf-47e3-8bf3-92b35854db2b
plot(m.beta')

# ╔═╡ 02866544-788a-40cd-a914-285c80bbbf42
modelcols(m.X.formulas.rhs,evts)

# ╔═╡ 5466ec0b-d28a-4be7-b686-cb5ad7fea92b


# ╔═╡ Cell order:
# ╠═25d21156-d006-11eb-3655-61be76e7db93
# ╠═34aefa9a-e034-4293-b53d-cc1dc349ecc7
# ╠═c73bf3ed-4af3-4e5a-997e-6399c9b85bbe
# ╠═d1e1de70-7e37-46a4-8ce2-3ca591bb84cf
# ╠═2e9ec6e0-03c6-41f6-bd3b-a44eaed5b553
# ╠═97443eb8-ac67-4f8e-9e8a-eb09176c8516
# ╠═2a5d5493-b88d-473f-8d81-9e9240b3ffe0
# ╠═5d23b74c-988d-48bc-9c30-c2af0dbabdb4
# ╠═3f2b7a92-a064-4be6-8a39-cabaa2453777
# ╠═7ef24e5a-63df-4e99-a021-7119d7e8d3bf
# ╠═12c617c1-1994-4efe-b8bc-23fd2325c2ca
# ╠═6f7bc2f6-28de-40a5-83ee-c15cd6f74f6c
# ╠═786d0114-e10e-4e29-b6d7-bbe635c55f08
# ╠═c4227637-81b4-4ea7-a7b7-741860fef8c3
# ╠═cd19784c-686e-41c5-bd32-cdf793cecf03
# ╠═29e875d3-45b9-4a6a-aaa6-c6eaf35244c5
# ╠═02b8ff44-c514-492e-bf94-0cbf44431097
# ╠═d964a5ca-7e63-41d5-bf97-526772d259e9
# ╠═54862fe8-be53-47c1-a249-e47bcd492ce4
# ╠═1ee404e1-3039-410d-99de-ce713e751280
# ╠═d27167fb-9cf4-4e0d-8c65-fde19577be0f
# ╠═5efce43e-7cb5-4af4-84ef-4dc2c82090a7
# ╠═cf1f63c7-8858-462a-809e-139d3872c989
# ╠═c4ebfcbb-4230-46d2-854f-f7f7c3e63a98
# ╠═be2758e1-83f6-4555-a0bf-2941687c3097
# ╠═07077e0e-87aa-437d-893e-9668e62f5d01
# ╠═ce39efac-c362-4d9f-b6c4-38e13251e68f
# ╠═454d82fd-4c6d-4f08-a516-4bb020bad421
# ╠═ba67232a-b1a1-4de3-9e38-dacb2ee81fa2
# ╠═826f0033-1add-4fc3-920a-7d86ef62d31e
# ╠═706be005-8330-4812-a596-6469f7a57201
# ╠═64735d5c-eff5-47cb-a843-6d0a941975c0
# ╠═6e6fa3e3-7c9b-4a32-a1a5-fb260f143f3d
# ╠═26f58723-c663-468d-b578-a59eff0d9acf
# ╠═fa514e1d-e251-4431-97f3-780da7fdfb95
# ╠═a009593b-f1f3-47ab-af5a-60612dfd37a9
# ╠═b5e5f89b-158b-44b0-8945-5d2f0bd79ef7
# ╠═07e724d3-4180-471d-8a08-9a0f0b869093
# ╠═ded80200-9a8c-47fe-9e10-c1a9e5ad08ac
# ╠═671696e3-68a0-45cd-a08d-aa6c7e400230
# ╠═98c4ffea-0228-4e2b-850e-92c64391bf61
# ╠═29ebdb59-57c2-4000-b037-2483e8c69d01
# ╠═dc4e4671-a744-4003-adb5-010a3e67628a
# ╠═70d14446-b3f0-4bf3-8b2b-e5836bd46633
# ╠═282ac7f1-1647-447f-af4c-371251524a37
# ╠═66121b10-930f-49a7-a222-2d1e5ea81343
# ╠═0fe93070-e434-43d8-ad95-521ef3d42402
# ╠═0aa995b2-1c7b-40b8-8913-b17852a290f9
# ╠═357bc67d-e581-41a0-a378-905072e55289
# ╠═d480d611-d53e-4f15-b2e9-a22ad993b64e
# ╠═9f4a59b6-a422-457a-a5d9-88beb65f7a78
# ╠═0181d139-7a38-4ba1-ac6b-0c3d76ab6a47
# ╠═1bfd0909-05b6-4782-a34d-1ea1ba649115
# ╠═52ec6011-ecf6-4b26-9bf5-0523f1ad035f
# ╠═d9a14d6d-5edf-47e3-8bf3-92b35854db2b
# ╠═02866544-788a-40cd-a914-285c80bbbf42
# ╠═5466ec0b-d28a-4be7-b686-cb5ad7fea92b
