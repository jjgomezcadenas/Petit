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

# ╔═╡ a0032301-8d33-4fa7-9401-1a0c47841f3f
begin
	using Revise
    using Petit
end

# ╔═╡ 349825ff-7ffe-4fa1-ba26-a772041f0323
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Plots
	using Printf
	using HDF5
	using Markdown
	using InteractiveUtils
	using Statistics
	using StatsBase
	using Distributions
	using StatsPlots
	#using DataFramesMeta
	
	import Glob
	#using Interpolations
	#using QuadGK
	#using LsqFit
	#using Chain
	#using Unitful 
end

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
      Pkg.instantiate()
  end

# ╔═╡ c44d77dd-22d5-486b-aa64-51da3acc4ec9
names(Petit)

# ╔═╡ 6c59aeae-7990-4b43-8378-0de210a3291a
md"""
# Read the Data 
"""

# ╔═╡ 5caa7c21-82e5-420b-b4e7-a0e33543b74b
md"Read the DF? $(@bind readdf CheckBox(default=true))"

# ╔═╡ 07b4e4c1-0469-41ca-9890-bb4f990b4645
#partdf

# ╔═╡ fccea553-c2b6-4133-bb45-847bf40276ee
#gamdf  = filter(:primary => ==(1), partdf)

# ╔═╡ 144d60a3-d70f-442b-a252-76178fecdbf7
md"""
# Hits data frame
"""

# ╔═╡ 6624fb61-e8ac-41e5-a602-c8765e21cede
md"""
## Hits global statistics
"""

# ╔═╡ 474b25f8-bd95-4389-867a-bb753dc77d45
md"""
## Event selection
"""

# ╔═╡ 6c58a746-da45-4e39-a178-91e060b2f34b
md"""
- Select event number
"""

# ╔═╡ aa71c155-9bed-4ad6-963e-575100094a9c
@bind nevent NumberField(0:1000000, default=1)


# ╔═╡ d1c1a8cb-a70e-480b-8f65-2cfa1a7022c3
md"""
### True hits
"""

# ╔═╡ 272c88a5-c6c0-41e1-b782-f8e23007eefc
md"""
#### Distances between voxels
"""

# ╔═╡ 7c20e495-0c39-44a5-9c40-8cd9c5d8e3de
md"""
#### Energy of voxels
"""

# ╔═╡ 1e6adccf-b640-4fd7-92cc-02cf61f27a24
md"""
#### Plot event 
"""

# ╔═╡ 37d7c197-0a10-4145-ab85-b5a22eae273e
md"""
## Voxelize hits
"""

# ╔═╡ ea1afc2f-9467-4b8b-a3e3-237925ac9fa1
md"""
- Select voxel size
"""

# ╔═╡ 5248add6-0d0c-446d-8490-0ccc009fd6e9
@bind vsize NumberField(0.0:0.1:10.0, default=5.0)

# ╔═╡ de61fd85-48cf-4c49-8abc-f929bdc4c3fc
md"""
### Select event
"""

# ╔═╡ c6bdf2b5-0b45-4e6a-af2a-b8f1c1ed4b78
md"""
#### Distances between voxels
"""

# ╔═╡ ee47c998-1649-4df2-bdc2-33bf53818e62
md"""
#### Energy of voxels
"""

# ╔═╡ 21f70e07-ef56-4b84-87df-dd4e0d468ed1
md"""
#### Plot event 
"""

# ╔═╡ a60ea1b2-1c4c-4565-ac34-d2580dd016e3
md"""
### Build Tracks
"""

# ╔═╡ 38cb58a7-a752-4d9e-bda9-afa14a9bcf5c
md"""
## Loop analysis
"""

# ╔═╡ 149507ef-73e4-43a9-8623-8b3d65bef952
md"Run analysis? $(@bind run_analysis CheckBox(default=false))"

# ╔═╡ f1b3c4ef-1b0c-4243-be00-1b0a15b6f3a2
if run_analysis
	md"events to run $(@bind events_to_run NumberField(0:1000000, default=100))" 
end

# ╔═╡ b19304c2-600a-42d3-a430-c205f75d8f77
if run_analysis
	md"max distance (mm) $(@bind max_distance_mm NumberField(0.0:1000.0, default=10.0))" 
end

# ╔═╡ 7d65660f-6b52-48e2-afea-0edd4de33faf
if run_analysis
	md"energy threshold (keV) $(@bind energy_threshold_kev NumberField(0.0:1000.0, default=10.0))" 
end

# ╔═╡ 0c725b8e-93f1-45d2-bc17-abbedcdaeab1
md"""
    struct AnalysisResults
        single_track_energies::Vector{Float64}
        two_track_primary::Vector{Float64}
        two_track_secondary::Vector{Float64}
        three_track_primary::Vector{Float64}
        three_track_secondary::Vector{Float64}
        n_events_processed::Int
        n_single_track::Int
        n_two_track::Int
        n_three_plus_track::Int
        n_failed::Int
end
"""

# ╔═╡ d0f9e26b-2a3c-42a4-a826-7c8694f5d470
md"""
# Functions
"""

# ╔═╡ 33efdcc5-cc1d-48f4-81af-b5f65bb348f7
function histogram_stats(hitsdf)
	h_id, p_id = step_hist(hitsdf.event_id;
         nbins = 50,
         xlim   = (0.0, 10000.0),
         xlabel = "event number",
         ylabel = "Frequency",
         title=" Events")
	
	nhits = hits_per_all_events(hitsdf)
	h_nhits, p_nhits = step_hist(nhits;
         nbins = 50,
         xlim   = (0.0, 500.0),
         xlabel = "number of hits per track",
         ylabel = "Frequency",
         title=" Hits per track")

	h_x, p_x = step_hist(hitsdf.x;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "X (mm)",
         ylabel = "Frequency",
         title=" X")

	h_y, p_y = step_hist(hitsdf.y;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "Y (mm)",
         ylabel = "Frequency",
         title=" Y")

	h_z, p_z = step_hist(hitsdf.z;
         nbins = 50,
         xlim   = (0.0, 4000.0),
         xlabel = "Z (mm)",
         ylabel = "Frequency",
         title=" Z ")

	histogram(hitsdf.x, label="X", bins=50, legend=:topleft, color=:gray)
	plot(p_nhits, p_x, p_y, p_z)
end

# ╔═╡ 6fc13d23-cb44-40a6-a52d-32b81855d6a0
function histogram_voxel_energy(df; emx=200.0)
	h_vxe, p_vxe = step_hist(1e+3*df.energy;
         nbins = 50,
         xlim   = (0.0, emx),
         xlabel = "Voxel energy (keV)",
         ylabel = "Frequency",
         title=" Voxel energy ")
	plot(p_vxe)
end

# ╔═╡ 2d4732eb-47b0-498b-b0ff-c9dcf3951ecd
function histogram_distances(df; dmx=100.0, dcmx=20.0)
	vd = voxel_distances(df)
	vcd = voxel_closest_distance(df)
	_, p_vd = step_hist(vd;
         nbins = 50,
         xlim   = (0.0, dmx),
         xlabel = "voxel distance (mm)",
         ylabel = "Frequency",
         title=" Voxel Distance ")
	_, p_vcd = step_hist(vcd;
         nbins = 50,
         xlim   = (0.0, dcmx),
         xlabel = "voxel closest distance (mm)",
         ylabel = "Frequency",
         title=" Voxel Closest Distance ")
	
	plot(p_vd, p_vcd)
end

# ╔═╡ 769f772d-2695-468f-93a0-7655b5d08f10
function histogram_energies_trks(results)
	
	_, p_t1 = step_hist(results.single_track_energies;
         nbins = 50,
         xlim   = (0.0, 2600.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E one track ")
	_, p_t2 = step_hist(results.two_track_primary;
         nbins = 50,
         xlim   = (0.0, 2600.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E two tracks ")
	_, p_t2s = step_hist(results.two_track_secondary;
         nbins = 50,
         xlim   = (0.0, 250.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E secondaries 2t ")
	_, p_t3s = step_hist(results.three_track_secondary;
         nbins = 50,
         xlim   = (0.0, 250.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E secondaries 3 track ")
	
	plot(p_t1, p_t2, p_t2s, p_t3s)
end

# ╔═╡ 4bb62775-b172-4537-9087-c37cf9ea97f4
function plot_true_hits2(ghdf::GroupedDataFrame, index::Int; nbins=100)
    df = get_event(ghdf, index)
	plot_hits_evt(df; nbins)
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

# ╔═╡ 6df1da0a-074f-485f-8a66-89645c7ab92e
names(hitsdf)

# ╔═╡ 0fa67f2a-50a5-4c82-b42b-07f45f14e914
hitsdf

# ╔═╡ e7771756-7cc6-45e7-8be5-98a3ec857cee
histogram_stats(hitsdf)

# ╔═╡ 6c148761-f9e2-4c7d-8d50-ee78bc0a8baf
evtdf = get_event(hitsdf, nevent)

# ╔═╡ 53884bc3-20d5-41af-9d0a-170820447f6f
histogram_distances(evtdf, dmx=100.0, dcmx=5.0)

# ╔═╡ 2416d5e3-8563-45c8-aafd-2f6528c50b85
histogram_voxel_energy(evtdf)

# ╔═╡ 1963d9dc-0cb9-4b8e-98ee-41491c6c784b
plot_hits_evt(hitsdf, nevent; nbins=100)

# ╔═╡ 8c929794-64a4-48eb-bd85-ccfad6133b1e
vdf = voxelize_hits(hitsdf, vsize)

# ╔═╡ d4892a5a-4488-4244-9513-8ece86d59661
begin
	vxdf = get_event(vdf, nevent)
	md"""
	- Voxelize event number $(nevent). Number of voxels = $(hits_per_event(vdf, nevent))
	"""
end

# ╔═╡ b978173f-d0e6-44fa-9447-15cd9fa7f3d9
histogram_distances(vxdf, dmx=100.0, dcmx=20.0)

# ╔═╡ bae05d6e-614b-4c37-ab64-0132d5f0e9dd
histogram_voxel_energy(vxdf, emx=200.0)

# ╔═╡ 9a89ffdc-3889-48c0-a952-ad2c36094b4c
plot_hits_evt(vdf, nevent; nbins=100) 

# ╔═╡ ceb6c086-fa7a-4265-92db-353b6380f184
begin
	trks = build_tracks(vdf, nevent; max_distance=5.0, energy_threshold=0.0001)
	md"""
	- number of traks found = $(length(trks))
	"""
end

# ╔═╡ d3b2a6f8-7c06-4512-8f28-096be0179c0a
if length(trks)>0
	histogram_voxel_energy(trks[1].voxels)
end

# ╔═╡ 9956f72c-ce38-4dc3-8217-65ddb2e1f723
if length(trks)>0
	plot_hits(trks[1].voxels, nbins=100)
end

# ╔═╡ df08eb9c-064d-4bcf-bde2-bcea5ef0d430
begin
      if length(trks) > 0
          ekev = 1e+3 * round(sum(trks[1].voxels.energy), digits=3)
          md"""
          Energy in main track = $(ekev) keV
          
          """
      else
          md"""
          No tracks found
          """
      end
  end

# ╔═╡ d98ea3d5-22bd-480d-9218-80cc7220bc07
if length(trks)>1
	trks[2].voxels
end

# ╔═╡ cc9ed738-cba3-49e2-a62f-c38da8aea5d8
if length(trks)>1
	histogram_voxel_energy(trks[2].voxels)
end

# ╔═╡ 59c79d07-52a7-482e-b4f8-f29a28401b8a
if length(trks)>2
	trks[3].voxels
end

# ╔═╡ 8fffd455-cf05-4a4c-b99e-7f81a2926373
if length(trks)>2
	histogram_voxel_energy(trks[3].voxels)
end

# ╔═╡ 4d0c3fe4-5e76-441d-87c7-3a840791257c
if run_analysis
	results = analysis_loop(hitsdf; 
				 events_to_run=1:events_to_run, 
				 voxel_size_mm=vsize,
				 max_distance_mm=max_distance_mm, 
				 energy_threshold_kev=energy_threshold_kev)

  

  # Access statistics
	md"""
	- Processed $(results.n_events_processed) events
	- Single track events: $(results.n_single_track)
	- Two track events: $(results.n_two_track)
	- Three+ track events: $(results.n_three_plus_track)
	- Failed events: $(results.n_failed)
	- Average energy in peak = $(mean(results.single_track_energies)) keV
	"""
  
end

# ╔═╡ faa0ee0d-447e-43f8-8ebd-d3e6c3c1df24
if run_analysis
	histogram_energies_trks(results)
end

# ╔═╡ Cell order:
# ╠═04b446d6-f34f-11ed-2565-0b15d65b6781
# ╠═871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
# ╠═947c237c-9852-40e9-a83f-c23666db90aa
# ╠═349825ff-7ffe-4fa1-ba26-a772041f0323
# ╠═a0032301-8d33-4fa7-9401-1a0c47841f3f
# ╠═c44d77dd-22d5-486b-aa64-51da3acc4ec9
# ╠═6c59aeae-7990-4b43-8378-0de210a3291a
# ╠═ae89f5dc-e958-496a-91ac-0bd977355563
# ╠═4af9a3ef-e883-4bc3-a2f1-212102e4951b
# ╠═5caa7c21-82e5-420b-b4e7-a0e33543b74b
# ╠═37a23226-4e3c-4e21-a4ee-a9aef75b2093
# ╠═8a55b4a3-5cbf-48c3-b150-2bd4ad73f440
# ╠═07b4e4c1-0469-41ca-9890-bb4f990b4645
# ╠═fccea553-c2b6-4133-bb45-847bf40276ee
# ╠═6df1da0a-074f-485f-8a66-89645c7ab92e
# ╠═144d60a3-d70f-442b-a252-76178fecdbf7
# ╠═0fa67f2a-50a5-4c82-b42b-07f45f14e914
# ╠═6624fb61-e8ac-41e5-a602-c8765e21cede
# ╠═e7771756-7cc6-45e7-8be5-98a3ec857cee
# ╠═474b25f8-bd95-4389-867a-bb753dc77d45
# ╠═6c58a746-da45-4e39-a178-91e060b2f34b
# ╠═aa71c155-9bed-4ad6-963e-575100094a9c
# ╠═6c148761-f9e2-4c7d-8d50-ee78bc0a8baf
# ╠═d1c1a8cb-a70e-480b-8f65-2cfa1a7022c3
# ╠═272c88a5-c6c0-41e1-b782-f8e23007eefc
# ╠═53884bc3-20d5-41af-9d0a-170820447f6f
# ╠═7c20e495-0c39-44a5-9c40-8cd9c5d8e3de
# ╠═2416d5e3-8563-45c8-aafd-2f6528c50b85
# ╠═1e6adccf-b640-4fd7-92cc-02cf61f27a24
# ╠═1963d9dc-0cb9-4b8e-98ee-41491c6c784b
# ╠═37d7c197-0a10-4145-ab85-b5a22eae273e
# ╠═ea1afc2f-9467-4b8b-a3e3-237925ac9fa1
# ╠═5248add6-0d0c-446d-8490-0ccc009fd6e9
# ╠═8c929794-64a4-48eb-bd85-ccfad6133b1e
# ╠═de61fd85-48cf-4c49-8abc-f929bdc4c3fc
# ╠═d4892a5a-4488-4244-9513-8ece86d59661
# ╠═c6bdf2b5-0b45-4e6a-af2a-b8f1c1ed4b78
# ╠═b978173f-d0e6-44fa-9447-15cd9fa7f3d9
# ╠═ee47c998-1649-4df2-bdc2-33bf53818e62
# ╠═bae05d6e-614b-4c37-ab64-0132d5f0e9dd
# ╠═21f70e07-ef56-4b84-87df-dd4e0d468ed1
# ╠═9a89ffdc-3889-48c0-a952-ad2c36094b4c
# ╠═a60ea1b2-1c4c-4565-ac34-d2580dd016e3
# ╠═ceb6c086-fa7a-4265-92db-353b6380f184
# ╠═d3b2a6f8-7c06-4512-8f28-096be0179c0a
# ╠═9956f72c-ce38-4dc3-8217-65ddb2e1f723
# ╠═df08eb9c-064d-4bcf-bde2-bcea5ef0d430
# ╠═d98ea3d5-22bd-480d-9218-80cc7220bc07
# ╠═cc9ed738-cba3-49e2-a62f-c38da8aea5d8
# ╠═59c79d07-52a7-482e-b4f8-f29a28401b8a
# ╠═8fffd455-cf05-4a4c-b99e-7f81a2926373
# ╠═38cb58a7-a752-4d9e-bda9-afa14a9bcf5c
# ╠═149507ef-73e4-43a9-8623-8b3d65bef952
# ╠═f1b3c4ef-1b0c-4243-be00-1b0a15b6f3a2
# ╠═b19304c2-600a-42d3-a430-c205f75d8f77
# ╠═7d65660f-6b52-48e2-afea-0edd4de33faf
# ╠═0c725b8e-93f1-45d2-bc17-abbedcdaeab1
# ╠═4d0c3fe4-5e76-441d-87c7-3a840791257c
# ╠═faa0ee0d-447e-43f8-8ebd-d3e6c3c1df24
# ╠═d0f9e26b-2a3c-42a4-a826-7c8694f5d470
# ╠═33efdcc5-cc1d-48f4-81af-b5f65bb348f7
# ╠═6fc13d23-cb44-40a6-a52d-32b81855d6a0
# ╠═2d4732eb-47b0-498b-b0ff-c9dcf3951ecd
# ╠═769f772d-2695-468f-93a0-7655b5d08f10
# ╠═4bb62775-b172-4537-9087-c37cf9ea97f4
# ╠═46d3b837-c9df-4e04-98ae-611054c9ad6d
# ╠═a4d68f10-b2f3-4b48-b3d0-2baac6f93ca1
# ╠═bbc7767d-2ab4-40ff-ad2a-e1f81b6ad368
# ╠═59c673be-1ead-488c-84bd-21bf14980313
# ╠═70914efe-c290-4f43-8c2a-f521281c977e
