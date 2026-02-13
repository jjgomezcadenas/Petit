### A Pluto.jl notebook ###
# v0.20.21

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

# ╔═╡ 7504d7aa-a780-4956-99a5-08a7f9a462b2
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

# ╔═╡ c9fc0547-0e73-4629-9909-e59c3d75169d
begin
	jn = ingredients(string(pdir,"/src/Petit.jl"))
end

# ╔═╡ 6c59aeae-7990-4b43-8378-0de210a3291a
md"""
# Read the Data 
"""

# ╔═╡ ae89f5dc-e958-496a-91ac-0bd977355563
let
	let
	readdir(cmdir)
	dirs = jn.Petit.getdirs(cmdir)
	md""" Select set  : $(@bind sdata Select(dirs))"""
end
end

# ╔═╡ 4af9a3ef-e883-4bc3-a2f1-212102e4951b
begin
	#xfile    = joinpath(cmdir, sdata)
end

# ╔═╡ 9ea7167f-493b-4715-930e-3718a3dc6216
xfile = "/Users/jjgomezcadenas/Data/HD5t/precdr/copperbkg/bi214/bi214_copper_endcaps_1.next.h5"

# ╔═╡ 5caa7c21-82e5-420b-b4e7-a0e33543b74b
md"Read the DF? $(@bind readdf CheckBox(default=true))"

# ╔═╡ 8a55b4a3-5cbf-48c3-b150-2bd4ad73f440
begin
	dfs = jn.Petit.get_dataset_dfs(xfile)
	hitsdf = dfs["hits"]
	partdf =dfs["particles"]
	md"""
	- Size of particle DF = $(size(partdf))
	- Size of hits DF = $(size(hitsdf))
	"""
end

# ╔═╡ 07b4e4c1-0469-41ca-9890-bb4f990b4645
partdf

# ╔═╡ f60557f4-113e-44ab-ab53-56e967f8fda8
md"""
- number of events in data frame = $(jn.Petit.number_of_events(partdf))
- energy of primary = $(round(1e+3*jn.Petit.energy_primary(partdf,1), digits=1)) keV
"""

# ╔═╡ f84ad689-77ff-4641-9f07-1f04ef718ccc
md"""
### Find events with alpha particles
"""

# ╔═╡ e3d1ed4b-a388-4965-83e9-f21e18c318ec
if occursin("bi214_surface", xfile) || occursin("tl208_surface", xfile) 
	alphadf = jn.Petit.find_events_with_alphas(partdf)
end

# ╔═╡ e3ec7177-e78a-42b7-9f0a-d904af5395a1


# ╔═╡ 49075c9c-7297-4226-ac34-88aeab6834e5
if occursin("bi214_surface", xfile) || occursin("tl208_surface", xfile) 
	md"""
	#### Detection of alpha particles
	- The plots above show where the apha particles were produced.  
	- For the events kept in the DF, the alpha particles end up in the copper or steel and are not detected
	"""
end

# ╔═╡ 144d60a3-d70f-442b-a252-76178fecdbf7
md"""
# Hits data frame
"""

# ╔═╡ 0fa67f2a-50a5-4c82-b42b-07f45f14e914
hitsdf

# ╔═╡ 4fb34131-0a18-453f-93c7-c31f34b959ce
md"""
#### Initial number of events in DF
- These are events contained in the gas which deposit between 2300 and 2700 keV
- number of events in DF =$(jn.Petit.number_of_events(hitsdf))
"""

# ╔═╡ 6624fb61-e8ac-41e5-a602-c8765e21cede
md"""
## Hits global statistics
"""

# ╔═╡ 19e48672-9ad7-405c-9706-60d4699bda2c
begin
h_y, p_y = jn.Petit.step_hist(hitsdf.y;
         nbins = 50,
         xlim   = (-2000.0, -1000.0),
         xlabel = "Y (mm)",
         ylabel = "Frequency",
         title=" Y")

	h_z, p_z = jn.Petit.step_hist(hitsdf.z;
         nbins = 50,
         xlim   = (-200.0, 200.0),
         xlabel = "Z (mm)",
         ylabel = "Frequency",
         title=" Z ")
	plot(p_y, p_z)
end

# ╔═╡ 5bcffb43-a8f7-448f-86c2-f3261d489bbf
md""" 
- Set the values of the fiducial cuts
"""

# ╔═╡ 35d5b13a-9769-410d-99bc-edd8e8cf15fb
@bind xyc NumberField(1500.0:10.0:2000.0, default=1800.0)

# ╔═╡ 490d05a8-a555-48d0-9e86-952e69ecaf8e
@bind zc NumberField(50.0:1.0:150.0, default=100.0)

# ╔═╡ ff6ef332-9b58-4dda-a997-489544d5636e
begin
	#xyc = 1800.0
	#zc = 100.0
md"""
### Fiducial cuts
- abs(x) < $(xyc)  
- abs(y) < $(xyc)
- abs(z) <$(zc)
"""
end

# ╔═╡ 292c8f4c-5739-4b5a-8125-00a9a0ac626f
begin
	hitsdf2 = jn.Petit.filter_fiducial_events(hitsdf, xyc, zc)
end

# ╔═╡ 1b54761b-e05c-4384-8b43-c3e2ac32eb34
jn.Petit.number_of_events(hitsdf2)

# ╔═╡ cea5dce3-4fdc-445d-abdd-117334946199
md"""
#### After fiducial cut
- number of events in DF =$(jn.Petit.number_of_events(hitsdf2))
"""

# ╔═╡ e1845fd0-9c43-48f7-a5f5-772f9ea50660
hitsdf2

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


# ╔═╡ 264be000-9352-4a3d-9238-384da3c06053
unique(hitsdf2.event_id)

# ╔═╡ 6c148761-f9e2-4c7d-8d50-ee78bc0a8baf
evtdf = jn.Petit.get_event(hitsdf2, nevent)

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

# ╔═╡ 1963d9dc-0cb9-4b8e-98ee-41491c6c784b
jn.Petit.plot_hits_evt(hitsdf, nevent; nbins=100)

# ╔═╡ 37d7c197-0a10-4145-ab85-b5a22eae273e
md"""
## Voxelize hits
"""

# ╔═╡ ea1afc2f-9467-4b8b-a3e3-237925ac9fa1
md"""
- Select voxel size
"""

# ╔═╡ 5248add6-0d0c-446d-8490-0ccc009fd6e9
@bind vsize NumberField(0.0:0.1:10.0, default=3.0)

# ╔═╡ 8c929794-64a4-48eb-bd85-ccfad6133b1e
vdf = jn.Petit.voxelize_hits(hitsdf2, vsize)

# ╔═╡ de61fd85-48cf-4c49-8abc-f929bdc4c3fc
md"""
### Select event
"""

# ╔═╡ d4892a5a-4488-4244-9513-8ece86d59661
begin
	vxdf = jn.Petit.get_event(vdf, nevent)
	md"""
	- Voxelize event number $(nevent). Number of voxels = $(jn.Petit.hits_per_event(vdf, nevent))
	"""
end

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

# ╔═╡ 9a89ffdc-3889-48c0-a952-ad2c36094b4c
jn.Petit.plot_hits_evt(vdf, nevent; nbins=100) 

# ╔═╡ a60ea1b2-1c4c-4565-ac34-d2580dd016e3
md"""
### Build Tracks
"""

# ╔═╡ 995555f1-dd99-4903-a3b7-9e48b62d675d
trk_energy_kev(trks, trkid) = 1e+3 * round(sum(trks[trkid].voxels.energy), digits=3)

# ╔═╡ 37110c55-a96e-4c4e-9ca1-9fb464865af0
md"max distance (mm) $(@bind max_distance_mm NumberField(0.0:1000.0, default=10.0))" 

# ╔═╡ d115f04b-1cf2-47b4-81fd-d15aa8d6cd97
md"energy threshold (keV) $(@bind energy_threshold_kev NumberField(0.0:1000.0, default=1.0))" 

# ╔═╡ ceb6c086-fa7a-4265-92db-353b6380f184
begin
	trks = jn.Petit.make_tracks(vdf, nevent;
						max_dist=energy_threshold_kev,
						energy_thr=energy_threshold_kev)

	md"""
	- number of traks found = $(length(trks))
	"""
end

# ╔═╡ d3b2a6f8-7c06-4512-8f28-096be0179c0a
if length(trks)>0
	md"""
	- Track number 1
	"""
end

# ╔═╡ 9956f72c-ce38-4dc3-8217-65ddb2e1f723
if length(trks)>0
	jn.Petit.plot_hits(trks[1].voxels, nbins=100)
end

# ╔═╡ df08eb9c-064d-4bcf-bde2-bcea5ef0d430
begin
      if length(trks) > 0
          ekev = 1e+3 * round(sum(trks[1].voxels.energy), digits=3)
          md"""
          Energy in track number 1 = $(trk_energy_kev(trks, 1)) keVV
          
          """
      else
          md"""
          No tracks found
          """
      end
  end

# ╔═╡ 3b33f43e-f1b0-4740-9946-72bfcd578dac
if length(trks)>0
	md"""
	- Track number 2
	"""
end

# ╔═╡ aba1fd92-ccc2-405a-afb5-80ad1c73baa5
if length(trks)>1
	jn.Petit.plot_hits(trks[2].voxels, nbins=100)
end

# ╔═╡ d7167ba5-5c81-431e-8bec-4f242ebb3b87

if length(trks) > 1
	  #ekev = 1e+3 * round(sum(trks[2].voxels.energy), digits=3)
	  md"""
	  Energy in track number 2 = $(trk_energy_kev(trks, 2)) keV
	  
	  """
else
	  md"""
	  No tracks found
	  """
end

# ╔═╡ be81c6f7-4736-4931-a780-673d0e62d2e8
if length(trks)>0
	md"""
	- Track number 3
	"""
end

# ╔═╡ d98ea3d5-22bd-480d-9218-80cc7220bc07
if length(trks)>2
	jn.Petit.plot_hits(trks[3].voxels, nbins=100)
else
          md"""
          - No tracks found
          """
end

# ╔═╡ d84ed988-a09e-4c80-86d8-89eb1c9d77e0
if length(trks) > 2
	  
	  md"""
	  Energy in track number 2 = $(trk_energy_kev(trks, 3)) keV
	  
	  """
else
	  md"""
	  - No tracks found
	  """
end

# ╔═╡ 727fa030-3d24-4864-9df1-ee157da78d20
begin
	result = jn.Petit.walk_track_from_extremes(trks[1])
	start_voxel, end_voxel = result.extremes
	start_voxel, end_voxel = result.extremes
  	path = result.path_indices
  	ordered_voxels = result.path_voxels
  	track_length = result.total_length

	md"""
	### Find track Extremes
	- start voxel: x = $(start_voxel.x), y = $(start_voxel.y), z = $(start_voxel.z)
	- end voxel: x = $(end_voxel.x), y = $(end_voxel.y), z = $(end_voxel.z)
	- track length L =$(track_length)
	"""
end

# ╔═╡ 0ba05450-c601-4c3e-8147-74c35c6ef997
jn.Petit.plot_track_with_extremes(trks[1], result;
                                 markersize_voxels=4.0,
                                 markersize_extremes=10.0,
                                 show_connections=false,
                                 alpha_connections=0.3)

# ╔═╡ 38cb58a7-a752-4d9e-bda9-afa14a9bcf5c
md"""
## Loop analysis
"""

# ╔═╡ 149507ef-73e4-43a9-8623-8b3d65bef952
md"Run analysis? $(@bind run_analysis CheckBox(default=false))"

# ╔═╡ f1b3c4ef-1b0c-4243-be00-1b0a15b6f3a2
if run_analysis
	md"events to run $(@bind events_to_run NumberField(0:1000000, default=1000))" 
end

# ╔═╡ 8af4e9c9-2cdc-417d-b681-2647ba807b86
xfile

# ╔═╡ 4d0c3fe4-5e76-441d-87c7-3a840791257c
if run_analysis
	results = jn.Petit.event_loop(cmdir; input_file= xfile,
					   events_to_run=events_to_run, 
				 	   voxel_size_mm=vsize,
				 	   max_distance_mm=max_distance_mm, 
				       energy_threshold_kev=energy_threshold_kev,
					   xyc=xyc,
                       zc=zc)


  # Access statistics
	md"""
	- Processed $(results.n_events_processed) events
	- Single track events: $(results.n_single_track)
	- Two track events: $(results.n_two_track)
	- Three+ track events: $(results.n_three_plus_track)
	- Failed events: $(results.n_failed)

	"""
  
end

# ╔═╡ 25533234-5ebf-409e-a6fa-4f52e9a4ac00
results

# ╔═╡ 7f94b90d-3955-4cb0-b2f4-fd37b66d3b20
deb = false

# ╔═╡ 84eac452-abc3-4570-9e58-791b75aa5041
if run_analysis && deb
	println("events with a single track")
	for (i, evt) in enumerate(results.single_track.ids)
		println("event =", evt, "energy (keV) =", results.single_track.energies[i])
	end
end

# ╔═╡ 1eaecf8f-e561-4bf2-af51-b010f924d3a9
if run_analysis
eht1, pht1 = jn.Petit.step_hist(results.single_track.energies;
         nbins = 40,
         xlim   = (2300.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E one track ")
	plot(pht1)
end

# ╔═╡ b5cbfe88-ada1-4160-af2d-0a7b9aab5964
md"energy resolution (keV) $(@bind erex NumberField(0.0:30.0, default=12.5))" 

# ╔═╡ b9995760-1dc7-47dd-9ec6-098f8a1d204c
md"roi low (keV) $(@bind roi_low NumberField(2475.0:2500.0, default=2475.0))" 

# ╔═╡ f64b47c1-0aba-4e9a-96a4-3938fb08279b
md"roi up (keV) $(@bind roi_up NumberField(2475.0:2500.0, default=2500.0))" 

# ╔═╡ a855a625-043c-4744-b247-1dd528ea1dae
if run_analysis
	eres = jn.Petit.smear_histogram(eht1, erex)
	ehrx, phrx = jn.Petit.step_hist(eres;
         nbins = 40,
         xlim   = (2300.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E (resolution) one track ")
	plot(phrx)
end

# ╔═╡ 842ecb21-35a7-463e-9691-9ee3bf15e231
jn.Petit.counts_in_range(ehrx, roi_low, roi_up)

# ╔═╡ 85d7771b-fc92-4844-9086-0d8c3b41e210
function signal_eff(ehrx, rlow, rup; step=10.0) 
	countx = []
	norm = sum(ehrx.weights)
	for rx in rlow:step:rup
		push!(countx, jn.Petit.counts_in_range(ehrx, rx, rup))
	end
	countx = countx /norm
end

# ╔═╡ 0528e076-9994-472c-9d44-4d5b8c2af754
if run_analysis
	step= 10.0
	xx = roi_low:step:roi_up
	countx = signal_eff(ehrx, roi_low, roi_up; step=step)
	plot(xx, countx)
end

# ╔═╡ d0f9e26b-2a3c-42a4-a826-7c8694f5d470
md"""
# Functions
"""

# ╔═╡ e18e2427-5818-4d3a-9dd8-e715aad2958f
begin
	function histogram_energy_primary(partdf)
	energies = 1e+3*jn.Petit.energy_primary(partdf)
	he, pe = jn.Petit.step_hist(energies;
                   nbins = 50,
                   xlim = (2300.0, 2700.0),
                   xlabel = "energy (keV)",
                   ylabel = "Frequency",
                   title ="Energy Primary")
	plot(pe)
end
function histogram_energy_deposited(hitsdf; emin=2300.0, emax=2700.0)
	energies = 1e+3*jn.Petit.energy_deposited(hitsdf)
	_, pe = jn.Petit.step_hist(energies;
                   nbins = 50,
				   logy=true, 
                   xlim = (emin, emax),
                   xlabel = "energy (keV)",
                   ylabel = "Frequency",
                   title ="Energy Deposited")
	plot(pe)
end
end

# ╔═╡ f6726067-9bb6-44c1-bcdf-2959798e1ce2

if occursin("bi214_copper", xfile) || occursin("tl208_copper", xfile)  || occursin("0nubb", xfile)
	histogram_energy_primary(partdf)
end

# ╔═╡ bb097911-6b23-4e21-a3d7-65fffe2f3bf7
histogram_energy_deposited(hitsdf, emax=2800.0)

# ╔═╡ 3984ae1a-f844-434d-841d-9a0f0557f85f
histogram_energy_deposited(hitsdf2, emax=2800.0)

# ╔═╡ beb867b5-711f-47ba-b6d1-c3a6e63dd6f1
function histogram_alphas(df)
	
	h_x, p_x = jn.Petit.step_hist(df.x;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "X (mm)",
         ylabel = "Frequency",
         title=" X")

	h_y, p_y = jn.Petit.step_hist(df.y;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "Y (mm)",
         ylabel = "Frequency",
         title=" Y")

	h_z, p_z = jn.Petit.step_hist(df.z;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "Z (mm)",
         ylabel = "Frequency",
         title=" Z ")

	h_e, p_e = jn.Petit.step_hist(1e+3*df.energy;
         nbins = 50,
         xlim   = (0.0, 8000.0),
         xlabel = "energy (V)",
         ylabel = "Frequency",
         title=" E (keV) ")

	plot(p_x, p_y, p_z, p_e)
end

# ╔═╡ 5502a45f-bba7-47b4-acad-6b14749fb4ed
if occursin("bi214_surface", xfile) || occursin("tl208_surface", xfile) 
	histogram_alphas(alphadf)
end

# ╔═╡ 33efdcc5-cc1d-48f4-81af-b5f65bb348f7
function histogram_stats(hitsdf)
	h_id, p_id = jn.Petit.step_hist(hitsdf.event_id;
         nbins = 50,
         xlim   = (0.0, 10000.0),
         xlabel = "event number",
         ylabel = "Frequency",
         title=" Events")
	
	nhits = jn.Petit.hits_per_all_events(hitsdf)
	h_nhits, p_nhits = jn.Petit.step_hist(nhits;
         nbins = 50,
         xlim   = (0.0, 500.0),
         xlabel = "number of hits per track",
         ylabel = "Frequency",
         title=" Hits per track")

	h_x, p_x = jn.Petit.step_hist(hitsdf.x;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "X (mm)",
         ylabel = "Frequency",
         title=" X")

	h_y, p_y = jn.Petit.step_hist(hitsdf.y;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "Y (mm)",
         ylabel = "Frequency",
         title=" Y")

	h_z, p_z = jn.Petit.step_hist(hitsdf.z;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "Z (mm)",
         ylabel = "Frequency",
         title=" Z ")

	plot(p_nhits, p_x, p_y, p_z)
end

# ╔═╡ e7771756-7cc6-45e7-8be5-98a3ec857cee
histogram_stats(hitsdf)

# ╔═╡ c1f44776-3eae-4006-bbf2-7280e60b8367
histogram_stats(hitsdf2)

# ╔═╡ 6fc13d23-cb44-40a6-a52d-32b81855d6a0
function histogram_voxel_energy(df; emx=200.0, title="Voxel energy")
	h_vxe, p_vxe = jn.Petit.step_hist(1e+3*df.energy;
         nbins = 50,
         xlim   = (0.0, emx),
         xlabel = "Voxel energy (keV)",
         ylabel = "Frequency",
         title=title)
	plot(p_vxe)
end

# ╔═╡ 2416d5e3-8563-45c8-aafd-2f6528c50b85
histogram_voxel_energy(evtdf)

# ╔═╡ bae05d6e-614b-4c37-ab64-0132d5f0e9dd
histogram_voxel_energy(vxdf, emx=200.0)

# ╔═╡ aa73dadb-fa63-40c1-98f3-bdf5cdfa835d
if length(trks)>0
	histogram_voxel_energy(trks[1].voxels)
end

# ╔═╡ ceafa987-894e-4ce1-af81-0fd90ba832dc
if length(trks)>1
	histogram_voxel_energy(trks[2].voxels)
end

# ╔═╡ 59c79d07-52a7-482e-b4f8-f29a28401b8a
if length(trks)>2
	histogram_voxel_energy(trks[2].voxels)
end

# ╔═╡ 2d4732eb-47b0-498b-b0ff-c9dcf3951ecd
function histogram_distances(df; dmx=100.0, dcmx=20.0)
	vd = jn.Petit.voxel_distances(df)
	vcd = jn.Petit.voxel_closest_distance(df)
	_, p_vd = jn.Petit.step_hist(vd;
         nbins = 50,
         xlim   = (0.0, dmx),
         xlabel = "voxel distance (mm)",
         ylabel = "Frequency",
         title=" Voxel Distance ")
	_, p_vcd = jn.Petit.step_hist(vcd;
         nbins = 50,
         xlim   = (0.0, dcmx),
         xlabel = "voxel closest distance (mm)",
         ylabel = "Frequency",
         title=" Voxel Closest Distance ")
	
	plot(p_vd, p_vcd)
end

# ╔═╡ 53884bc3-20d5-41af-9d0a-170820447f6f
histogram_distances(evtdf, dmx=100.0, dcmx=5.0)

# ╔═╡ b978173f-d0e6-44fa-9447-15cd9fa7f3d9
histogram_distances(vxdf, dmx=100.0, dcmx=20.0)

# ╔═╡ 769f772d-2695-468f-93a0-7655b5d08f10
function histogram_energies_trks(results)
	
	_, p_t1 = jn.Petit.step_hist(results.single_track.energies;
         nbins = 50,
         xlim   = (0.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E one track ")
	_, p_t2 = jn.Petit.step_hist(results.two_track_primary.energies;
         nbins = 50,
         xlim   = (0.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E two tracks ")
	_, p_t2s = jn.Petit.step_hist(results.two_track_secondary.energies;
         nbins = 50,
         xlim   = (0.0, 250.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E secondaries 2t ")
	_, p_t3s = jn.Petit.step_hist(results.three_track_secondary.energies;
         nbins = 50,
         xlim   = (0.0, 250.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E secondaries 3 track ")
	
	plot(p_t1, p_t2, p_t2s, p_t3s)
end

# ╔═╡ faa0ee0d-447e-43f8-8ebd-d3e6c3c1df24
if run_analysis
	histogram_energies_trks(results)
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

# ╔═╡ f436311b-b47a-486b-aed0-170c6b7dbc94
function inspect_mc2(filename::String)
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
inspect_mc(xfile)

# ╔═╡ 2019c5f3-b700-4dc4-bc95-06f8cb1e966a
h5open(xfile, "r") do fid
      config_array = read(fid["MC/configuration"])
      param_value = config_array[2].param_value  # Assuming structured array
  end

# ╔═╡ ae3fb816-c8cb-4b37-abf7-e3a978045c2d
function nof_events(xfile)
	
	fid = h5open(xfile, "r") 
    config_array = read(fid["MC/configuration"])
    param_value = config_array[2].param_value  # Assuming structured array
	return parse(Int, param_value)
end

# ╔═╡ c52c7194-dff1-4a97-ac62-8abdd4df2555
nof = nof_events(xfile)

# ╔═╡ 71b5e9cc-94d7-4fcf-8498-6e422d4e8416


# ╔═╡ Cell order:
# ╠═947c237c-9852-40e9-a83f-c23666db90aa
# ╠═04b446d6-f34f-11ed-2565-0b15d65b6781
# ╠═871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
# ╠═349825ff-7ffe-4fa1-ba26-a772041f0323
# ╠═7504d7aa-a780-4956-99a5-08a7f9a462b2
# ╠═c9fc0547-0e73-4629-9909-e59c3d75169d
# ╠═6c59aeae-7990-4b43-8378-0de210a3291a
# ╠═ae89f5dc-e958-496a-91ac-0bd977355563
# ╠═4af9a3ef-e883-4bc3-a2f1-212102e4951b
# ╠═9ea7167f-493b-4715-930e-3718a3dc6216
# ╠═5caa7c21-82e5-420b-b4e7-a0e33543b74b
# ╠═c52c7194-dff1-4a97-ac62-8abdd4df2555
# ╠═8a55b4a3-5cbf-48c3-b150-2bd4ad73f440
# ╠═07b4e4c1-0469-41ca-9890-bb4f990b4645
# ╠═f60557f4-113e-44ab-ab53-56e967f8fda8
# ╠═f6726067-9bb6-44c1-bcdf-2959798e1ce2
# ╠═f84ad689-77ff-4641-9f07-1f04ef718ccc
# ╠═e3d1ed4b-a388-4965-83e9-f21e18c318ec
# ╠═e3ec7177-e78a-42b7-9f0a-d904af5395a1
# ╠═5502a45f-bba7-47b4-acad-6b14749fb4ed
# ╠═49075c9c-7297-4226-ac34-88aeab6834e5
# ╠═144d60a3-d70f-442b-a252-76178fecdbf7
# ╠═bb097911-6b23-4e21-a3d7-65fffe2f3bf7
# ╠═0fa67f2a-50a5-4c82-b42b-07f45f14e914
# ╠═4fb34131-0a18-453f-93c7-c31f34b959ce
# ╠═6624fb61-e8ac-41e5-a602-c8765e21cede
# ╠═e7771756-7cc6-45e7-8be5-98a3ec857cee
# ╠═19e48672-9ad7-405c-9706-60d4699bda2c
# ╠═5bcffb43-a8f7-448f-86c2-f3261d489bbf
# ╠═35d5b13a-9769-410d-99bc-edd8e8cf15fb
# ╠═490d05a8-a555-48d0-9e86-952e69ecaf8e
# ╠═ff6ef332-9b58-4dda-a997-489544d5636e
# ╠═292c8f4c-5739-4b5a-8125-00a9a0ac626f
# ╠═1b54761b-e05c-4384-8b43-c3e2ac32eb34
# ╠═c1f44776-3eae-4006-bbf2-7280e60b8367
# ╠═cea5dce3-4fdc-445d-abdd-117334946199
# ╠═e1845fd0-9c43-48f7-a5f5-772f9ea50660
# ╠═3984ae1a-f844-434d-841d-9a0f0557f85f
# ╠═474b25f8-bd95-4389-867a-bb753dc77d45
# ╠═6c58a746-da45-4e39-a178-91e060b2f34b
# ╠═aa71c155-9bed-4ad6-963e-575100094a9c
# ╠═264be000-9352-4a3d-9238-384da3c06053
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
# ╠═995555f1-dd99-4903-a3b7-9e48b62d675d
# ╠═37110c55-a96e-4c4e-9ca1-9fb464865af0
# ╠═d115f04b-1cf2-47b4-81fd-d15aa8d6cd97
# ╠═ceb6c086-fa7a-4265-92db-353b6380f184
# ╠═d3b2a6f8-7c06-4512-8f28-096be0179c0a
# ╠═9956f72c-ce38-4dc3-8217-65ddb2e1f723
# ╠═df08eb9c-064d-4bcf-bde2-bcea5ef0d430
# ╠═aa73dadb-fa63-40c1-98f3-bdf5cdfa835d
# ╠═3b33f43e-f1b0-4740-9946-72bfcd578dac
# ╠═aba1fd92-ccc2-405a-afb5-80ad1c73baa5
# ╠═d7167ba5-5c81-431e-8bec-4f242ebb3b87
# ╠═ceafa987-894e-4ce1-af81-0fd90ba832dc
# ╠═be81c6f7-4736-4931-a780-673d0e62d2e8
# ╠═d98ea3d5-22bd-480d-9218-80cc7220bc07
# ╠═d84ed988-a09e-4c80-86d8-89eb1c9d77e0
# ╠═59c79d07-52a7-482e-b4f8-f29a28401b8a
# ╠═727fa030-3d24-4864-9df1-ee157da78d20
# ╠═0ba05450-c601-4c3e-8147-74c35c6ef997
# ╠═38cb58a7-a752-4d9e-bda9-afa14a9bcf5c
# ╠═149507ef-73e4-43a9-8623-8b3d65bef952
# ╠═f1b3c4ef-1b0c-4243-be00-1b0a15b6f3a2
# ╠═8af4e9c9-2cdc-417d-b681-2647ba807b86
# ╠═4d0c3fe4-5e76-441d-87c7-3a840791257c
# ╠═25533234-5ebf-409e-a6fa-4f52e9a4ac00
# ╠═faa0ee0d-447e-43f8-8ebd-d3e6c3c1df24
# ╠═7f94b90d-3955-4cb0-b2f4-fd37b66d3b20
# ╠═84eac452-abc3-4570-9e58-791b75aa5041
# ╠═1eaecf8f-e561-4bf2-af51-b010f924d3a9
# ╠═b5cbfe88-ada1-4160-af2d-0a7b9aab5964
# ╠═b9995760-1dc7-47dd-9ec6-098f8a1d204c
# ╠═f64b47c1-0aba-4e9a-96a4-3938fb08279b
# ╠═a855a625-043c-4744-b247-1dd528ea1dae
# ╠═842ecb21-35a7-463e-9691-9ee3bf15e231
# ╠═85d7771b-fc92-4844-9086-0d8c3b41e210
# ╠═0528e076-9994-472c-9d44-4d5b8c2af754
# ╠═d0f9e26b-2a3c-42a4-a826-7c8694f5d470
# ╠═e18e2427-5818-4d3a-9dd8-e715aad2958f
# ╠═beb867b5-711f-47ba-b6d1-c3a6e63dd6f1
# ╠═33efdcc5-cc1d-48f4-81af-b5f65bb348f7
# ╠═6fc13d23-cb44-40a6-a52d-32b81855d6a0
# ╠═2d4732eb-47b0-498b-b0ff-c9dcf3951ecd
# ╠═769f772d-2695-468f-93a0-7655b5d08f10
# ╠═4bb62775-b172-4537-9087-c37cf9ea97f4
# ╠═46d3b837-c9df-4e04-98ae-611054c9ad6d
# ╠═a4d68f10-b2f3-4b48-b3d0-2baac6f93ca1
# ╠═bbc7767d-2ab4-40ff-ad2a-e1f81b6ad368
# ╠═f436311b-b47a-486b-aed0-170c6b7dbc94
# ╠═59c673be-1ead-488c-84bd-21bf14980313
# ╠═2019c5f3-b700-4dc4-bc95-06f8cb1e966a
# ╠═ae3fb816-c8cb-4b37-abf7-e3a978045c2d
# ╠═71b5e9cc-94d7-4fcf-8498-6e422d4e8416
