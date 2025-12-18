### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 58616677-9ed4-4e8c-a196-59855e07da1f
begin
	using Markdown
	using CSV
	using DataFrames
	using Plots
	using Printf
	using HDF5
	using Statistics
	using StatsBase
	using Distributions
	using StatsPlots
	using NearestNeighbors
	using Interpolations 
	using Clustering
	using LinearAlgebra
	using Graphs
	using Peaks
	import Glob
end

# ╔═╡ 2c9f00a8-d81b-11f0-94c3-df0e14c11a77
begin
	cmdir=joinpath(ENV["DATA"], "HD5t/itaca")
	pdir =joinpath(ENV["PROJECTS"], "Petit")
end

# ╔═╡ 29b12a39-dcd0-482a-b8a7-7a20e6886fea
begin
      using Pkg
      Pkg.activate(pdir)
      Pkg.instantiate()
  end

# ╔═╡ f32e3c11-2701-4f55-a6ee-866e81f80b4b
function ingredients(path::String)
    # this is from the Julia source code (evalfile in base/loading.jl)
    # but with the modification that it returns the module instead of the last object
    name = Symbol(basename(path))
    m = Module(name)
    Core.eval(m,
        Expr(:toplevel,
                :(eval(x) = $(Expr(:core, :eval))($name, x)),
                :(include(x) = $(Expr(:top, :include))($name, x)),
                :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
                :(include($path))))
    m
end

# ╔═╡ 23879a19-4886-42e9-a293-1bd90d6074d9
begin
	jn = ingredients(string(pdir,"/src/Petit.jl"))
	
end

# ╔═╡ a6197a11-3fd2-4495-89fc-9dc5443927b4
it = ingredients(string(pdir,"/pluto/kde_blob_analysis_functions.jl"))

# ╔═╡ 73dd0b82-9619-4ca7-b9eb-07872a3e6699
begin
	ionL100 = it.get_data("/Users/jjgomezcadenas/Projects/Petit/ItacaResultsLdrft100/", data_type="ion")
	ionL100.bbm
end

# ╔═╡ ebaffbca-d45f-42fd-9721-ea0bb83bfb8a
begin 
	md"""
	## KDE analysis
	"""
end

# ╔═╡ 954f0c21-9682-4d48-9335-daf36a8cc3a1
md"""
- top: peak1 left position vs peak1 prom for bb0nu and xe137. For bb0nu one expect a peak of high prom near the left edege (cero) corresponding to one of the blobs. For xe137 (single electrons) the position of the blob changes along the trajectory
- bottom: peak2 right position vs peak2 prom, same reasoning. 
- One clear cut is to retain events with peak1 left position below 50 and peak2 right position above 150.
- But also notice that peak2 prom is higher for bb0nu than for xe. Thus peak2 prom > 0.5 is also efficient. 
"""

# ╔═╡ e8b35de7-4274-4c20-b65d-163c91b87252
it.plot_peak_pos_vs_prom(ionL100.bb, ionL100.xe)

# ╔═╡ ad9bb9a7-cd38-4446-a0c0-0e0a3d0312b5
it.plot_kdef_eff_and_fom(ionL100.bb, ionL100.xe; lmax = 50, rmin=150)

# ╔═╡ f6b8664b-59af-495c-b44b-42f1081cfb4d
md"""
- Fixing prom1 left position (lmax) to 50 and prom2 right position (rmin) to 150 the fom for pmin maximizes around 0.4
"""

# ╔═╡ a733c55c-5284-47e7-a5c5-6f18210a3185
begin
	pmin = 0.4
	lmax = 50
	rmin = 150
	ionL100K = it.apply_ked_cut(ionL100; pmin = pmin, lmax = lmax, rmin = rmin)
	md"""
		- #### KDE cuts 
		- pmin= $(pmin) lmax = $(lmax) rmin= $(rmin)
		- KDE eff for bb = $(ionL100K.effbb)
		- KDE eff for Xe = $(ionL100K.effxe)
		"""
end

# ╔═╡ 0a64722e-9f88-449a-9424-fbaa7866442d
md"""
## BLOB analysis
"""

# ╔═╡ 44ee1537-98a2-449e-a389-2e589cb26c8f
let
	bbpeb = it.plot_eb1_vs_eb2(ionL100.bb, title="bb0nu: Eb1 vs Eb2")
	xepeb = it.plot_eb1_vs_eb2(ionL100.xe, title="xe: Eb1 vs Eb2")
	bbpeb2 = it.plot_eb1_vs_eb2(ionL100K.bb, title="bb0nu/KDE: Eb1 vs Eb2")
	xepeb2 = it.plot_eb1_vs_eb2(ionL100K.xe, title="xe/KDE: Eb1 vs Eb2")
	plot(bbpeb, xepeb,bbpeb2,xepeb2, layout=(2,2), size=(1200, 1000), margin=6Plots.mm)
end

# ╔═╡ ee8d81ff-36ea-4671-a4f0-43b325ab2b4d
md"""
### Direct cut on eblob2
"""

# ╔═╡ 262d86d5-1aa8-4346-bdf2-df9cbc152c57
it.plot_eb2_cut_eff_and_fom(ionL100.bb,ionL100.xe)

# ╔═╡ 7d269cf9-7a55-4704-be2f-20f4b9fb22a5
begin
	eb2cut=400.0
	ionL100B = it.apply_blob_cut(ionL100; eb2cut=eb2cut)
	
	md"""
	- #### BLOB cuts 
	- eb2cut = $(eb2cut)
	- BLOB eff for bb = $(ionL100B.effbb)
	- BLOB eff for Xe = $(ionL100B.effxe)
	"""
end

# ╔═╡ 520e44d4-ddde-4545-8a34-fcedf560063f
begin 
ionL100KB = it.apply_kde_and_blob_cut(ionL100; pmin = pmin, 
										   lmax = lmax, 
										   rmin = rmin, 
										   eb2cut=eb2cut)

	md"""
		- #### KDE +  BLOB cuts 
		- eff for bb = $(ionL100KB.effbb)
		- eff for Xe = $(ionL100KB.effxe)
		"""
end

# ╔═╡ 888096e1-1499-40ca-b5f5-edc3b8cf1431
md"""
### KDE vs BLOB
- Combining with KDE does not improve BLOB analysis. 
"""

# ╔═╡ 6acc1fac-6219-45cc-91a8-bc951e341f39
md"""
### ROC Curve
"""

# ╔═╡ 4b59ebe2-b2bf-4117-b8c0-b4afd8763046
jn.Petit.plot_roc(ionL100.bb.Eb2_keV,
                  ionL100.xe.Eb2_keV;
                  cuts = range(0, 800.0, length=500),
                  title = "ion ROC Curve",
                  show_auc = true,
                  show_diagonal = false)

# ╔═╡ 8abbf09e-ac9f-4ea3-8a21-261f2b27de97
ionL100B.xe

# ╔═╡ e2645903-dc16-46e4-9752-a9224b67bdfe
xe_events = ionL100B.xe.event

# ╔═╡ 7f15c6ca-04da-4beb-a6dc-851e0c980a7b
ionL100B.bb

# ╔═╡ 76ec2a91-c83d-4ab2-9f1f-0a3d61b12869
bb_events = ionL100B.bb.event

# ╔═╡ 53218032-59e1-4b47-a8f2-6f69f222a941
function write_vectors(v, filename)
	open(filename, "w") do f
	    for x in v
	        println(f, x)
	    end
	end
end

# ╔═╡ 4bc5a6a9-edf3-4435-b550-1c298c472b70
write_vectors(xe_events, "xe_events.txt")

# ╔═╡ 6498e73d-bbf8-44eb-83f3-72895314605c
write_vectors(bb_events[1:87], "bb_events.txt")

# ╔═╡ d0913fb9-2b91-4ec1-b58d-bb7debfa6d5f
md"""
## Post analysis
"""

# ╔═╡ 867401c3-b0c8-4e75-95d4-f3bc9114692e
dfbb0nu = CSV.read("bb0nu/itaca_single_track_analysis_results.csv", DataFrame)

# ╔═╡ d5a7c941-4a0b-4780-ac23-b10e7a55c824
dfxe137 = CSV.read("xe137/itaca_single_track_analysis_results.csv", DataFrame)

# ╔═╡ 5fe72311-5d94-4f51-a889-1f2bd3a7341d
let
	p1 = histogram(dfbb0nu.asymmetry)
	p2 = histogram(dfxe137.asymmetry)
	plot(p1,p2)
end

# ╔═╡ 57e4ae59-c4d1-44e2-92a4-2a80d0d2231d
function plot_eb1_vs_eb2(df; 
                         title="Blob Energies: Eb1 vs Eb2",
                         xlabel="Eb1 (keV)", 
                         ylabel="Eb2 (keV)",
                         markersize=4,
                         alpha=0.6,
                         add_diagonal=true)
    
    p = scatter(df.Eb1_keV, df.Eb2_keV,
                xlabel=xlabel,
                ylabel=ylabel,
                title=title,
                markersize=markersize,
                alpha=alpha,
                legend=false,
                aspect_ratio=:equal)
    
    if add_diagonal
        # Add diagonal line (Eb1 = Eb2)
        lims = (min(minimum(df.Eb1_keV), minimum(df.Eb2_keV)),
                max(maximum(df.Eb1_keV), maximum(df.Eb2_keV)))
        plot!(p, [lims[1], lims[2]], [lims[1], lims[2]], 
              linestyle=:dash, color=:gray, label="Eb1=Eb2")
    end
    
    return p
end

# ╔═╡ 8fc7aa5f-c4b3-4382-92ac-43a3ea04f726
let
sb12bb =plot_eb1_vs_eb2(dfbb0nu; 
                title="Blob Energies: Eb1 vs Eb2",
                xlabel="Eb2 (keV)", 
                ylabel="Eb1 (keV)",
                markersize=2)
	sb12xe =plot_eb1_vs_eb2(dfxe137; 
                title="Blob Energies: Eb1 vs Eb2",
                xlabel="Eb2 (keV)", 
                ylabel="Eb1 (keV)",
                markersize=2)
	plot(sb12bb, sb12xe)
end

# ╔═╡ Cell order:
# ╠═2c9f00a8-d81b-11f0-94c3-df0e14c11a77
# ╠═29b12a39-dcd0-482a-b8a7-7a20e6886fea
# ╠═58616677-9ed4-4e8c-a196-59855e07da1f
# ╠═f32e3c11-2701-4f55-a6ee-866e81f80b4b
# ╠═23879a19-4886-42e9-a293-1bd90d6074d9
# ╠═a6197a11-3fd2-4495-89fc-9dc5443927b4
# ╠═73dd0b82-9619-4ca7-b9eb-07872a3e6699
# ╠═ebaffbca-d45f-42fd-9721-ea0bb83bfb8a
# ╠═954f0c21-9682-4d48-9335-daf36a8cc3a1
# ╠═e8b35de7-4274-4c20-b65d-163c91b87252
# ╠═ad9bb9a7-cd38-4446-a0c0-0e0a3d0312b5
# ╠═f6b8664b-59af-495c-b44b-42f1081cfb4d
# ╠═a733c55c-5284-47e7-a5c5-6f18210a3185
# ╠═0a64722e-9f88-449a-9424-fbaa7866442d
# ╠═44ee1537-98a2-449e-a389-2e589cb26c8f
# ╠═ee8d81ff-36ea-4671-a4f0-43b325ab2b4d
# ╠═262d86d5-1aa8-4346-bdf2-df9cbc152c57
# ╠═7d269cf9-7a55-4704-be2f-20f4b9fb22a5
# ╠═520e44d4-ddde-4545-8a34-fcedf560063f
# ╠═888096e1-1499-40ca-b5f5-edc3b8cf1431
# ╠═6acc1fac-6219-45cc-91a8-bc951e341f39
# ╠═4b59ebe2-b2bf-4117-b8c0-b4afd8763046
# ╠═8abbf09e-ac9f-4ea3-8a21-261f2b27de97
# ╠═e2645903-dc16-46e4-9752-a9224b67bdfe
# ╠═7f15c6ca-04da-4beb-a6dc-851e0c980a7b
# ╠═76ec2a91-c83d-4ab2-9f1f-0a3d61b12869
# ╠═53218032-59e1-4b47-a8f2-6f69f222a941
# ╠═4bc5a6a9-edf3-4435-b550-1c298c472b70
# ╠═6498e73d-bbf8-44eb-83f3-72895314605c
# ╠═d0913fb9-2b91-4ec1-b58d-bb7debfa6d5f
# ╠═867401c3-b0c8-4e75-95d4-f3bc9114692e
# ╠═d5a7c941-4a0b-4780-ac23-b10e7a55c824
# ╠═8fc7aa5f-c4b3-4382-92ac-43a3ea04f726
# ╠═5fe72311-5d94-4f51-a889-1f2bd3a7341d
# ╠═57e4ae59-c4d1-44e2-92a4-2a80d0d2231d
