### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 0c7e101a-7aed-4ae3-9e4e-f8a5688b873f
using Pkg; Pkg.activate(ENV["JPetit"])

# ╔═╡ 1f214aba-7811-4fed-b238-8c882c99c7c0
begin
	#using PlutoUI
	#using CSV
	using DataFrames
end

# ╔═╡ cacba748-f349-11ed-2e15-5b3f79935a35
#begin
 #   import Pkg
  #  Pkg.activate(mktempdir())
   # Pkg.add([Pkg.PackageSpec(name="DataFrames", version="1.5.0")])
    #using DataFrames
#end

# ╔═╡ 8d84062f-d315-4eaf-bde7-f1f1f023b731
x = DataFrame(id=rand('a':'d', 100), v=rand(100))

# ╔═╡ 60a7a88b-fb8a-42c9-993c-7261f2a50666
combine(groupby(x, :id)) do sdf
    n = size(sdf)[1]
    n < 25 ? DataFrame() : DataFrame(n=n) # drop groups with low number of rows
end

# ╔═╡ Cell order:
# ╠═0c7e101a-7aed-4ae3-9e4e-f8a5688b873f
# ╠═1f214aba-7811-4fed-b238-8c882c99c7c0
# ╠═cacba748-f349-11ed-2e15-5b3f79935a35
# ╠═8d84062f-d315-4eaf-bde7-f1f1f023b731
# ╠═60a7a88b-fb8a-42c9-993c-7261f2a50666
