### A Pluto.jl notebook ###
# v0.19.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try
            Base.loaded_modules[Base.PkgId(
                Base.UUID("6e696c72-6542-2067-7265-42206c756150"),
                "AbstractPlutoDingetjes",
            )].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 360367c3-a17d-4c04-b8eb-e2f0366baa86
import Pkg;

# ╔═╡ a69bf247-acbc-4b4d-889d-3f52ff046ef4
Pkg.activate("/store/users/ehinger/unfoldjl_dev/")

# ╔═╡ cf3dcd54-cfa0-11eb-2e99-756f9df71a4d
using Revise, Unfold, Plots, BSplines

# ╔═╡ a9e6971f-8a50-4564-ace5-1fabe65b1522
import PlutoUI

# ╔═╡ 8a62a281-d91d-4e70-b301-4ebafdd58b38


# ╔═╡ 539da607-f382-4eaa-93ad-31c906aa67ee


# ╔═╡ 17943406-ac61-49bd-bc8d-7d3f62441b8e


# ╔═╡ 89eace50-c0f5-4951-821f-cf010cdb9726
@bind nsplines PlutoUI.Slider(2:10; default = 8, show_value = true)

# ╔═╡ cb62fb85-364d-4142-b28b-c8b0871b4ccc
@bind srate PlutoUI.Slider(1:50; default = 10, show_value = true)

# ╔═╡ 8496a69b-20bc-45c2-92cb-179d0ad7127c
@bind τ1 PlutoUI.Slider(-5:0.1:2; default = -0.5, show_value = true)

# ╔═╡ c8285ed7-bffb-4047-89ff-22a664c19a78
@bind τ2 PlutoUI.Slider(-2:0.1:5; default = 2, show_value = true)

# ╔═╡ 9df0d50d-2a1a-4431-9725-12178c824d62
τ = [τ1, τ2]

# ╔═╡ 94a258ea-febd-4793-89c4-d3755e82c48c
begin
    function splinebasis(τ, sfreq, nsplines, name::String)
        τ = Unfold.round_times(τ, sfreq)
        times = range(τ[1], stop = τ[2], step = 1 ./ sfreq)
        kernel = e -> splinekernel(e, times, nsplines)
        type = "splinebasis"

        shiftOnset = Int64(floor(τ[1] * sfreq))

        return Unfold.BasisFunction(kernel, times, times, type, name, shiftOnset)
    end



    function spl_breakpoints(times, nsplines)
        # calculate the breakpoints, evenly spaced
        return collect(range(minimum(times), stop = maximum(times), length = nsplines))
    end
    function splinekernel(e, times, nsplines)
        breakpoints = spl_breakpoints(times, nsplines)
        basis = BSplineBasis(4, breakpoints)  # 4= cubic
        return Unfold.splFunction(times, basis)
    end
    spl = splinebasis(τ, srate, nsplines, "test")

end

# ╔═╡ 3f612e26-aba8-46ed-af19-67cd1be672fa
plot(spl.times, spl.kernel(0))

# ╔═╡ d375b64c-0e64-4b39-8c3f-5a14d365877d


# ╔═╡ 7086e800-d51a-4fe3-b33f-309acd676ae7


# ╔═╡ 4d6fc512-205e-4e25-b3b5-b7981d9dec75


# ╔═╡ ba9bafc2-9733-43d0-80cd-75b2637ba4a7


# ╔═╡ Cell order:
# ╠═360367c3-a17d-4c04-b8eb-e2f0366baa86
# ╠═a9e6971f-8a50-4564-ace5-1fabe65b1522
# ╠═a69bf247-acbc-4b4d-889d-3f52ff046ef4
# ╠═8a62a281-d91d-4e70-b301-4ebafdd58b38
# ╠═539da607-f382-4eaa-93ad-31c906aa67ee
# ╠═17943406-ac61-49bd-bc8d-7d3f62441b8e
# ╠═94a258ea-febd-4793-89c4-d3755e82c48c
# ╠═3f612e26-aba8-46ed-af19-67cd1be672fa
# ╠═89eace50-c0f5-4951-821f-cf010cdb9726
# ╠═cb62fb85-364d-4142-b28b-c8b0871b4ccc
# ╠═8496a69b-20bc-45c2-92cb-179d0ad7127c
# ╠═c8285ed7-bffb-4047-89ff-22a664c19a78
# ╠═9df0d50d-2a1a-4431-9725-12178c824d62
# ╠═d375b64c-0e64-4b39-8c3f-5a14d365877d
# ╠═cf3dcd54-cfa0-11eb-2e99-756f9df71a4d
# ╠═7086e800-d51a-4fe3-b33f-309acd676ae7
# ╠═4d6fc512-205e-4e25-b3b5-b7981d9dec75
# ╠═ba9bafc2-9733-43d0-80cd-75b2637ba4a7
