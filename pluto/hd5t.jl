### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 349825ff-7ffe-4fa1-ba26-a772041f0323
begin
	using Revise
	using PlutoUI
	using CSV
	using DataFrames
	using Plots
	using Printf
	#using Interpolations
	#using QuadGK
	using Markdown
	using InteractiveUtils
	#using LsqFit
	using Statistics
	#using Chain
	using StatsBase
	using Distributions
	#using Unitful 
	using StatsPlots
	#using DataFramesMeta
	using HDF5
	import Glob
end

# ╔═╡ a0032301-8d33-4fa7-9401-1a0c47841f3f
using Petit

# ╔═╡ e987802c-a598-4bdd-8bff-3e6d5cf50df2
using histos

# ╔═╡ 04b446d6-f34f-11ed-2565-0b15d65b6781
PlutoUI.TableOfContents(title="HD5t analysis", indent=true)


# ╔═╡ 871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
begin
	cmdir=joinpath(ENV["DATA"], "HD5t")
	pdir =joinpath(ENV["PROJECTS"], "Petit")
end

# ╔═╡ 947c237c-9852-40e9-a83f-c23666db90aa
begin
      using Pkg
      Pkg.activate(pdir)
      #Pkg.develop(path=pdir)
  end

# ╔═╡ c44d77dd-22d5-486b-aa64-51da3acc4ec9
names(Petit)

# ╔═╡ 6c59aeae-7990-4b43-8378-0de210a3291a
md"""
# Read the Data Frame and fix the data 
"""

# ╔═╡ 5caa7c21-82e5-420b-b4e7-a0e33543b74b
md"Read the DF? $(@bind readdf CheckBox(default=true))"

# ╔═╡ d0f9e26b-2a3c-42a4-a826-7c8694f5d470
md"""
# Functions
"""

# ╔═╡ 42f9d64e-ca4e-4fa8-ab7f-de5f78726183
function euclidean_distance(x1, y1, z1, x2, y2, z2)
    return sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
end

# ╔═╡ 466b600b-034d-4159-b59c-8c05ebaaef36


# ╔═╡ 43a5b188-b8bc-4c15-83fa-fc64d537df20



# ╔═╡ 971ee5d1-ceaf-4c9b-b951-4080d7823e98

function plot_true_hits2(ghdf::GroupedDataFrame, index::Int; nbins=100)
    df = ghdf[index]

    x = Float64.(df.x)
    y = Float64.(df.y)
    z = Float64.(df.z)
    e = Float64.(df.energy)

    # Compute padded limits (1.3x range)
    xmid, xrange = mean((minimum(x), maximum(x))), maximum(x) - minimum(x)
    ymid, yrange = mean((minimum(y), maximum(y))), maximum(y) - minimum(y)
    zmid, zrange = mean((minimum(z), maximum(z))), maximum(z) - minimum(z)

    xlim = (xmid - 0.65 * xrange, xmid + 0.65 * xrange)
    ylim = (ymid - 0.65 * yrange, ymid + 0.65 * yrange)
    zlim = (zmid - 0.65 * zrange, zmid + 0.65 * zrange)

    cmap = cgrad(:viridis, alpha=1.0)

    # === XY ===
    h_xy = fit(Histogram, (x, y), nbins=nbins)
    wxy = h_xy.weights
    wxy_masked = map(v -> v == 0.0 ? NaN : v, wxy)
    xcenters_xy = diff(h_xy.edges[1]) ./ 2 .+ h_xy.edges[1][1:end-1]
    ycenters_xy = diff(h_xy.edges[2]) ./ 2 .+ h_xy.edges[2][1:end-1]
    p1 = heatmap(xcenters_xy, ycenters_xy, wxy_masked';
        xlabel="x", ylabel="y", title="XY Heatmap",
        xlims=xlim, ylims=ylim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === XZ ===
    h_xz = fit(Histogram, (x, z), nbins=nbins)
    wxz_masked = map(v -> v == 0.0 ? NaN : v, h_xz.weights)
    xcenters_xz = diff(h_xz.edges[1]) ./ 2 .+ h_xz.edges[1][1:end-1]
    zcenters_xz = diff(h_xz.edges[2]) ./ 2 .+ h_xz.edges[2][1:end-1]
    p2 = heatmap(xcenters_xz, zcenters_xz, wxz_masked';
        xlabel="x", ylabel="z", title="XZ Heatmap",
        xlims=xlim, ylims=zlim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === YZ ===
    h_yz = fit(Histogram, (y, z), nbins=nbins)
    wyz_masked = map(v -> v == 0.0 ? NaN : v, h_yz.weights)
    ycenters_yz = diff(h_yz.edges[1]) ./ 2 .+ h_yz.edges[1][1:end-1]
    zcenters_yz = diff(h_yz.edges[2]) ./ 2 .+ h_yz.edges[2][1:end-1]
    p3 = heatmap(ycenters_yz, zcenters_yz, wyz_masked';
        xlabel="y", ylabel="z", title="YZ Heatmap",
        xlims=ylim, ylims=zlim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === 3D scatter ===
    p4 = scatter(x, y, z, marker_z=e, ms=2,
        xlabel="x", ylabel="y", zlabel="z", title="3D Scatter",
        xlims=xlim, ylims=ylim, zlims=zlim,
        colorbar_title="Energy", legend=false, cgrad=cmap)

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
end

# ╔═╡ 8d93ac33-3f2c-4dc6-80ac-6160548cc1bf
function plot_true_hits(ghdf::GroupedDataFrame, index::Int; nbins=100)
    df = ghdf[index]

    x = Float64.(df.x)
    y = Float64.(df.y)
    z = Float64.(df.z)
    e = Float64.(df.energy)

    ### XY ###
    h_xy = fit(Histogram, (x, y), nbins=nbins)
    xcenters_xy = diff(h_xy.edges[1]) ./ 2 .+ h_xy.edges[1][1:end-1]
    ycenters_xy = diff(h_xy.edges[2]) ./ 2 .+ h_xy.edges[2][1:end-1]
    p1 = heatmap(xcenters_xy, ycenters_xy, h_xy.weights',
        xlabel="x", ylabel="y", title="XY Heatmap")

    ### XZ ###
    h_xz = fit(Histogram, (x, z), nbins=nbins)
    xcenters_xz = diff(h_xz.edges[1]) ./ 2 .+ h_xz.edges[1][1:end-1]
    zcenters_xz = diff(h_xz.edges[2]) ./ 2 .+ h_xz.edges[2][1:end-1]
    p2 = heatmap(xcenters_xz, zcenters_xz, h_xz.weights',
        xlabel="x", ylabel="z", title="XZ Heatmap")

    ### YZ ###
    h_yz = fit(Histogram, (y, z), nbins=nbins)
    ycenters_yz = diff(h_yz.edges[1]) ./ 2 .+ h_yz.edges[1][1:end-1]
    zcenters_yz = diff(h_yz.edges[2]) ./ 2 .+ h_yz.edges[2][1:end-1]
    p3 = heatmap(ycenters_yz, zcenters_yz, h_yz.weights',
        xlabel="y", ylabel="z", title="YZ Heatmap")

    ### 3D scatter ###
    p4 = scatter(x, y, z, marker_z=e, ms=2,
        xlabel="x", ylabel="y", zlabel="z", title="3D Scatter", colorbar_title="Energy", legend=false)

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(900, 700))
end

# ╔═╡ 46d3b837-c9df-4e04-98ae-611054c9ad6d
function get_group_contents(filename::String; path::String = "/")
    h5open(filename, "r") do file
        if !haskey(file, path)
            error("Path $path not found in HDF5 file.")
        end
        group = file[path]
        return Dict(name => typeof(group[name]) for name in keys(group))
    end
end

# ╔═╡ a4d68f10-b2f3-4b48-b3d0-2baac6f93ca1
function get_subgroups(filename::String, path::String = "/")
    h5open(filename, "r") do file
        if !haskey(file, path)
            error("Path $path not found in HDF5 file.")
        end
        group = file[path]
        return [name for name in keys(group) if group[name] isa HDF5.Group]
    end
end

# ╔═╡ bbc7767d-2ab4-40ff-ad2a-e1f81b6ad368
function inspect_mc(filename::String)
    h5open(filename, "r") do fid
        dset = fid["MC"]
        println("Type: ", typeof(dset))
        data = read(dset)
        println("Read type: ", typeof(data))
        println("First element: ", data[1])
        return data
    end
end

# ╔═╡ 3a0af47d-c297-4324-bfe8-a267b7cedb6d
#function get_dfs(filename::String)

#	fid     = h5open(filename, "r")
	#dset =fid["MC"]["vertices"]
#	dset =fid["MC"]
#	dda = read(dset)
#	ddic = Dict()
 #   for (i, key) in enumerate(keys(dda[1]))
  #      ddic[key] = [ddi[i] for ddi in dda]
   # end
	#close(fid)
    #DataFrame(ddic)
#end

# ╔═╡ 59c673be-1ead-488c-84bd-21bf14980313
function get_dataset_dfs(filename::String)
    h5open(filename, "r") do fid
        group = fid["MC"]
        dfs = Dict{String, DataFrame}()

        for name in keys(group)
            data = read(group[name])
            # Try to make a DataFrame if it's an array of structs or tuples
            try
                dfs[name] = DataFrame(data)
            catch
                dfs[name] = DataFrame((value=data,))
            end
        end

        return dfs  # Dict of DataFrames
    end
end

# ╔═╡ 70914efe-c290-4f43-8c2a-f521281c977e
function getdirs(bdir::AbstractString)
	fdrs = Glob.glob("*", bdir)
	[split(f,"/")[end] for f in fdrs]
end

# ╔═╡ ae89f5dc-e958-496a-91ac-0bd977355563
let
	let
	readdir(cmdir)
	dirs = getdirs(cmdir)
	md""" Select set  : $(@bind sdata Select(dirs))"""
end
end

# ╔═╡ 4af9a3ef-e883-4bc3-a2f1-212102e4951b
begin
	xfile    = joinpath(cmdir, sdata)
end

# ╔═╡ 37a23226-4e3c-4e21-a4ee-a9aef75b2093
if readdf
	get_subgroups(xfile)
end

# ╔═╡ 8a55b4a3-5cbf-48c3-b150-2bd4ad73f440
begin
	dfs = get_dataset_dfs(xfile)
	hitsdf = dfs["hits"]
	partdf =dfs["particles"]
	md"""
	- Size of particle DF = $(size(partdf))
	- Size of hits DF = $(size(hitsdf))
	"""
end

# ╔═╡ 07b4e4c1-0469-41ca-9890-bb4f990b4645
partdf

# ╔═╡ fccea553-c2b6-4133-bb45-847bf40276ee
gamdf  = filter(:primary => ==(1), partdf)

# ╔═╡ 6df1da0a-074f-485f-8a66-89645c7ab92e
names(hitsdf)

# ╔═╡ 0fa67f2a-50a5-4c82-b42b-07f45f14e914
hitsdf

# ╔═╡ 59226b94-f112-4840-be05-8c72493423ba
histogram(hitsdf.event_id, label="events", bins=100, legend=:topleft, color=:gray)

# ╔═╡ 99b43220-0ce1-41e2-a1f5-447387c51bdf
begin
    ghdf = groupby(hitsdf, :event_id)
    counts = Vector{Int}(undef, length(ghdf))  # preallocate vector

    for (i, subdf) in enumerate(ghdf)
        #println("Event ID: ", subdf.event_id[1])
        #println("Number of entries: ", nrow(subdf))
        counts[i] = nrow(subdf)
    end
	histogram(counts, label="hits per event", bins=100, legend=:topleft, color=:gray)
      # return the result
end

# ╔═╡ 1963d9dc-0cb9-4b8e-98ee-41491c6c784b
plot_true_hits2(ghdf, 7, nbins=100)

# ╔═╡ 457854d7-41c5-4b2b-bee2-b7ce4e3bfffb
begin
	energy_per_event = combine(ghdf, :energy => sum => :total_energy)
	histogram(energy_per_event.total_energy, label="evergy/event", bins=100, legend=:topleft, color=:gray)
end

# ╔═╡ 8c929794-64a4-48eb-bd85-ccfad6133b1e
vdf = voxelize_hits(ghdf, 1.0)

# ╔═╡ 2396b88b-80cc-4a8b-9bb6-d6d5098b1ac5
begin
	gvdf = groupby(vdf, :event_id)
    counts2 = Vector{Int}(undef, length(gvdf))  # preallocate vector

    for (i, subdf) in enumerate(gvdf)
        #println("Event ID: ", subdf.event_id[1])
        #println("Number of entries: ", nrow(subdf))
        counts2[i] = nrow(subdf)
    end
	histogram(counts2, label="hits per event voxelized", bins=100, legend=:topleft, color=:gray)
      # return the result
end

# ╔═╡ d26db229-60a3-48e3-82ec-869c0b05572b
plot_true_hits2(gvdf, 1, nbins=100)

# ╔═╡ a4b5d213-43eb-4810-9059-23c37572aa58
histogram_voxel_distances(gvdf; bins=100, maxevts=100, max_distance=100)

# ╔═╡ 730f2af9-03d0-45aa-a71f-42ba7641a5cf
histogram_closest_distance(gvdf; bins=100, maxevts=100, max_distance=10.0)

# ╔═╡ ceb6c086-fa7a-4265-92db-353b6380f184
trks = build_vgraph(gvdf, 1; max_distance=2.5, energy_threshold=0.001)

# ╔═╡ e2058dc6-d3ad-479f-b4f1-16627a05b213
typeof(trks[1])

# ╔═╡ d98ea3d5-22bd-480d-9218-80cc7220bc07
trks[1].components

# ╔═╡ 6fd2a640-9232-46a6-8464-8e95bfc77641
histogram(hitsdf.x, label="X", bins=50, legend=:topleft, color=:gray)

# ╔═╡ 8ae015df-6212-4113-82a5-7752aa7b0fda
histogram(hitsdf.y, label="Y", bins=50, legend=:topleft, color=:gray)

# ╔═╡ 51916a3a-6c4f-4988-874a-b8e24b4cbe3f
histogram(hitsdf.z, label="Z", bins=50, legend=:topleft, color=:gray)

# ╔═╡ a1b2c3d4-e5f6-7890-1234-567890abcdef
function histogram_voxel_energy(gvdf::GroupedDataFrame; bins=50, title="Voxel Energy Distribution")
    all_energies = Float64[]
    
    for group in gvdf
        append!(all_energies, group.energy)
    end
    
    return histogram(all_energies, bins=bins, xlabel="Energy", ylabel="Count", 
                    title=title, legend=false, color=:blue, alpha=0.7)
end

# ╔═╡ f15e5308-0dbf-4ce2-88ab-b16119642c6e
histogram_voxel_energy(gvdf::GroupedDataFrame; bins=50, title="Voxel Energy Distribution")

# ╔═╡ Cell order:
# ╠═04b446d6-f34f-11ed-2565-0b15d65b6781
# ╠═871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
# ╠═947c237c-9852-40e9-a83f-c23666db90aa
# ╠═349825ff-7ffe-4fa1-ba26-a772041f0323
# ╠═a0032301-8d33-4fa7-9401-1a0c47841f3f
# ╠═c44d77dd-22d5-486b-aa64-51da3acc4ec9
# ╠═e987802c-a598-4bdd-8bff-3e6d5cf50df2
# ╟─6c59aeae-7990-4b43-8378-0de210a3291a
# ╠═ae89f5dc-e958-496a-91ac-0bd977355563
# ╠═4af9a3ef-e883-4bc3-a2f1-212102e4951b
# ╠═5caa7c21-82e5-420b-b4e7-a0e33543b74b
# ╠═37a23226-4e3c-4e21-a4ee-a9aef75b2093
# ╠═8a55b4a3-5cbf-48c3-b150-2bd4ad73f440
# ╠═07b4e4c1-0469-41ca-9890-bb4f990b4645
# ╠═fccea553-c2b6-4133-bb45-847bf40276ee
# ╠═6df1da0a-074f-485f-8a66-89645c7ab92e
# ╠═0fa67f2a-50a5-4c82-b42b-07f45f14e914
# ╠═59226b94-f112-4840-be05-8c72493423ba
# ╠═99b43220-0ce1-41e2-a1f5-447387c51bdf
# ╠═6fd2a640-9232-46a6-8464-8e95bfc77641
# ╠═8ae015df-6212-4113-82a5-7752aa7b0fda
# ╠═51916a3a-6c4f-4988-874a-b8e24b4cbe3f
# ╠═1963d9dc-0cb9-4b8e-98ee-41491c6c784b
# ╠═457854d7-41c5-4b2b-bee2-b7ce4e3bfffb
# ╠═8c929794-64a4-48eb-bd85-ccfad6133b1e
# ╠═2396b88b-80cc-4a8b-9bb6-d6d5098b1ac5
# ╠═d26db229-60a3-48e3-82ec-869c0b05572b
# ╠═f15e5308-0dbf-4ce2-88ab-b16119642c6e
# ╠═a4b5d213-43eb-4810-9059-23c37572aa58
# ╠═730f2af9-03d0-45aa-a71f-42ba7641a5cf
# ╠═ceb6c086-fa7a-4265-92db-353b6380f184
# ╠═e2058dc6-d3ad-479f-b4f1-16627a05b213
# ╠═d98ea3d5-22bd-480d-9218-80cc7220bc07
# ╠═d0f9e26b-2a3c-42a4-a826-7c8694f5d470
# ╠═42f9d64e-ca4e-4fa8-ab7f-de5f78726183
# ╠═466b600b-034d-4159-b59c-8c05ebaaef36
# ╠═43a5b188-b8bc-4c15-83fa-fc64d537df20
# ╠═971ee5d1-ceaf-4c9b-b951-4080d7823e98
# ╠═8d93ac33-3f2c-4dc6-80ac-6160548cc1bf
# ╠═46d3b837-c9df-4e04-98ae-611054c9ad6d
# ╠═a4d68f10-b2f3-4b48-b3d0-2baac6f93ca1
# ╠═bbc7767d-2ab4-40ff-ad2a-e1f81b6ad368
# ╠═3a0af47d-c297-4324-bfe8-a267b7cedb6d
# ╠═59c673be-1ead-488c-84bd-21bf14980313
# ╠═70914efe-c290-4f43-8c2a-f521281c977e
# ╠═a1b2c3d4-e5f6-7890-1234-567890abcdef
