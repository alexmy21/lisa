### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ e64fa1f2-0cd0-46e2-8b6d-ea96c78c6105
include("/home/alexmy/JULIA/lisa/src/lisa.jl")

# ╔═╡ bf167655-060c-4544-b81b-cc4bd36215ea
hll = HllSet{5}()

# ╔═╡ e18afb8a-d6c2-497e-80fe-e7433c2088b2
push!(hll, Set(["dima","mash","yan"]))

# ╔═╡ 1c0177a9-f654-40f7-a996-39c5e64e6cfe
count(hll)

# ╔═╡ Cell order:
# ╠═e64fa1f2-0cd0-46e2-8b6d-ea96c78c6105
# ╠═bf167655-060c-4544-b81b-cc4bd36215ea
# ╠═e18afb8a-d6c2-497e-80fe-e7433c2088b2
# ╠═1c0177a9-f654-40f7-a996-39c5e64e6cfe
