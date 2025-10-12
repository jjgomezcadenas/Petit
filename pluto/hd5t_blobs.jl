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

# ╔═╡ b604fcf3-7355-43f5-9b85-4c9c1ff77cd7
begin
      df_bi214 = CSV.read("bi214_fm.csv", DataFrame)
      df_tl208 = CSV.read("tl208_fm.csv", DataFrame)
      df_e0bb = CSV.read("electron_0nubb_fm.csv", DataFrame)

      # Plot with both points+error bars and lines
      plot(df_bi214.eblob2, df_bi214.eff,
           line=:solid, lw=1, alpha=0.5,
           label=nothing, color=:blue)
      scatter!(df_bi214.eblob2, df_bi214.eff,
               yerror=df_bi214.err,
               label="Bi-214",
               markersize=4,
               color=:blue,
               markerstrokewidth=1)

      plot!(df_tl208.eblob2, df_tl208.eff,
            line=:solid, lw=1, alpha=0.5,
            label=nothing, color=:orange)
      scatter!(df_tl208.eblob2, df_tl208.eff,
               yerror=df_tl208.err,
               label="Tl-208",
               markersize=4,
               color=:orange,
               markerstrokewidth=1)

      plot!(df_e0bb.eblob2, df_e0bb.eff,
            line=:solid, lw=1, alpha=0.5,
            label=nothing, color=:green)
      scatter!(df_e0bb.eblob2, df_e0bb.eff,
               yerror=df_e0bb.err,
               label="single electron",
               markersize=4,
               color=:green,
               markerstrokewidth=1)

      xlabel!("Blob 2 Energy Cut (keV)")
      ylabel!("Efficiency")
      title!("Efficiency vs Energy Cut")
      ylims!(0, 0.1)
      plot!(grid=true, gridalpha=0.3)  # Correct way to add 
end

# ╔═╡ 9ea7167f-493b-4715-930e-3718a3dc6216
xfile = "/Users/jjgomezcadenas/Data/HD5t/precdr/copperbkg/bi214/bi214_copper_endcaps_1.next.h5"

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

# ╔═╡ 3984ae1a-f844-434d-841d-9a0f0557f85f
#histogram_energy_deposited(hitsdf2, emax=2800.0)

# ╔═╡ 474b25f8-bd95-4389-867a-bb753dc77d45
#md"""
## Event selection
#"""

# ╔═╡ 6c58a746-da45-4e39-a178-91e060b2f34b
md"""
- Select event number
"""

# ╔═╡ aa71c155-9bed-4ad6-963e-575100094a9c
@bind nevent NumberField(0:1000000, default=1)


# ╔═╡ 2416d5e3-8563-45c8-aafd-2f6528c50b85
#histogram_voxel_energy(evtdf)

# ╔═╡ ea1afc2f-9467-4b8b-a3e3-237925ac9fa1
md"""
- Select voxel size
"""

# ╔═╡ 5248add6-0d0c-446d-8490-0ccc009fd6e9
@bind vsize NumberField(0.0:0.1:10.0, default=3.0)

# ╔═╡ de61fd85-48cf-4c49-8abc-f929bdc4c3fc
md"""
### Select event
"""

# ╔═╡ d4892a5a-4488-4244-9513-8ece86d59661
#begin
#	vdf = jn.Petit.voxelize_hits(hitsdf2, vsize)
#	vxdf = jn.Petit.get_event(vdf, nevent)
#	md"""
#	- Voxelize event number $(nevent). Number of voxels = #$(jn.Petit.hits_per_event(vdf, nevent))
#	"""
#end

# ╔═╡ c6bdf2b5-0b45-4e6a-af2a-b8f1c1ed4b78
md"""
#### Distances between voxels
"""

# ╔═╡ b978173f-d0e6-44fa-9447-15cd9fa7f3d9
#histogram_distances(vxdf, dmx=100.0, dcmx=20.0)

# ╔═╡ ee47c998-1649-4df2-bdc2-33bf53818e62
md"""
#### Energy of voxels
"""

# ╔═╡ bae05d6e-614b-4c37-ab64-0132d5f0e9dd
#histogram_voxel_energy(vxdf, emx=200.0)

# ╔═╡ 21f70e07-ef56-4b84-87df-dd4e0d468ed1
md"""
#### Plot event 
"""

# ╔═╡ 9a89ffdc-3889-48c0-a952-ad2c36094b4c
#jn.Petit.plot_hits_evt(vdf, nevent; nbins=100) 

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

# ╔═╡ 38cb58a7-a752-4d9e-bda9-afa14a9bcf5c
md"""
## Loop analysis
"""

# ╔═╡ 149507ef-73e4-43a9-8623-8b3d65bef952
md"Run analysis? $(@bind run_analysis CheckBox(default=false))"

# ╔═╡ f1b3c4ef-1b0c-4243-be00-1b0a15b6f3a2
if run_analysis
	md"events to run $(@bind events_to_run NumberField(0:1000000, default=10))" 
end

# ╔═╡ 8af4e9c9-2cdc-417d-b681-2647ba807b86
xfile

# ╔═╡ c0d93b0b-43f9-4d23-b180-ba813ed7414d
"""
I want to write a script with the following functionality:
1. Reads an h5 file specified by directory and file name. Example: 
cmdir = /Users/jjgomezcadenas/Data/HD5t/precdr/copperbkg/bi214/
input_file = bi214_copper_endcaps_1.next.h5

2. finds out how many events are in the file (ntot). Also take as input the maximum number of events the use wants to run (nmax), Define nmx = minimum(ntot, nmax). Take nmax as negative by default. If nmax is negative, then nmx = ntot.

3. Open an h5 file (specified by out_file parameter in the same directory)

4. Takes as parameter the integer nbatch, which is the number of batches in the loop

for ievent in 1:nbatch:nmx
...
end

4. Inside the loop, calls function event_loop_single_track(cmdir; input_file="0nubb.next.h5",
                    events_to_run=nbatch,
                    voxel_size_mm=5,
                    max_distance_mm=10,
                    energy_threshold_kev=10,
                    xyc::Float64=1950.0,
                    zc::Float64=10.0)

Be careful with the last call, since the number of events may be smaller than nbatch.

5. event_loop_single_track returns a Track structure, with the single tracks found in the batch. Access the h5 file and write the tracks. 

6. Close the file.
7. Provide a function to read back the output file and perform some statistics

"""

# ╔═╡ 4d0c3fe4-5e76-441d-87c7-3a840791257c
if run_analysis
	strks = jn.Petit.event_loop_single_track(cmdir; input_file= xfile,
					   events_to_run=events_to_run, 
				 	   voxel_size_mm=vsize,
				 	   max_distance_mm=max_distance_mm, 
				       energy_threshold_kev=energy_threshold_kev,
					   xyc=xyc,
                       zc=zc)

	md"""
	number of tracks = $(length(strks))
	"""
  
end

# ╔═╡ 88aec10d-6110-4408-9037-29007fc6c5d1
strks

# ╔═╡ 9bae1645-1485-411d-a6d3-61906c7c4194
begin
	i = 3 
	xresult = jn.Petit.walk_track_from_extremes(strks[i])
	xstart_voxel, xend_voxel = xresult.extremes
  	xtrack_length = xresult.total_length
	energy_kev = 1e+3 * sum(strks[i].voxels.energy)

	md"""
	### Find track Extremes
	- confidence = $(xresult.confidence)
	- start voxel: x = $(xstart_voxel.x), y = $(xstart_voxel.y), z = $(xstart_voxel.z)
	- end voxel: x = $(xend_voxel.x), y = $(xend_voxel.y), z = $(xend_voxel.z)
	- track length L =$(xtrack_length)
	- energy = $(energy_kev) KeV
	"""
end

# ╔═╡ 03cd5a27-1ace-43b3-8a31-4e34629b082d
#strks[i]

# ╔═╡ 6137bec5-66e7-431f-ab9b-a9cfdf5762c6
#jn.Petit.plot_track_with_extremes(strks[i], xresult;
#                                 markersize_voxels=4.0,
#                                 markersize_extremes=10.0,
#                                 show_connections=false,
#                                 alpha_connections=0.3)

# ╔═╡ bf77b535-485a-4608-aa4b-1c15950f3f9b
begin
	r= 5.0 # mm
	blobs = jn.Petit.energy_in_spheres_around_extremes(strks[i], xresult, r)
	eb1 = blobs.blob1_energy * 1e+3
	eb2 = blobs.blob2_energy * 1e+3
	nb1 = blobs.blob1_voxel_count
	nb2 = blobs.blob2_voxel_count
	md"""
	- blob 1 energy = $(round(eb1, digits=1)) keV
	- blob 2 energy = $(round(eb2, digits=1)) keV
	- blob 1 # of voxels = $(nb1)
	- blob 2 # of voxels = $(nb2)
	"""
end

# ╔═╡ 1c202304-d6ca-4dde-a7ac-f5cb2b961b93
jn.Petit.plot_track_blobs(strks[i], r;
                         markersize_voxels=3.0,
                         show_connections=true,
                         alpha_connections=0.2,
                         alpha_spheres=0.3,
                         sphere_resolution=20)

# ╔═╡ 6c57e536-6815-437a-aba6-64a47a499666
 #jn.Petit.save_track_with_analysis(strks[i], "track_bends_sharp")

# ╔═╡ 98033f21-5720-4d22-94a8-a72e650a860b
function blob_analysis(strks, r)
	eB1=Float64[]
	eB2=Float64[]
	CON = Float64[]
	E = Float64[]
	TL = Float64[]
	for strk in strks
		xresult = jn.Petit.walk_track_from_extremes(strk)
	  	xtrack_length = xresult.total_length
		confidence = xresult.confidence
		energy_kev = 1e+3 * sum(strks[i].voxels.energy)
		
		blobs = jn.Petit.energy_in_spheres_around_extremes(strk, xresult, r)
		eb1 = blobs.blob1_energy * 1e+3
		eb2 = blobs.blob2_energy * 1e+3
		push!(eB1, eb1)
		push!(eB2, eb2)
		push!(CON, confidence)
		push!(E, energy_kev)
		push!(TL, xtrack_length)
	end
	return CON, TL, E, eB1, eB2
		
	
end

# ╔═╡ 403b3c7f-d7c3-4762-8605-5d24fb087415
xcon, xtl, xe, xeb1, xeb2 = blob_analysis(strks, 15.0)

# ╔═╡ b59257e9-9a5e-467f-bbc3-f04cdf0e9a75
begin
	h_c, p_c = jn.Petit.step_hist(xcon *1.0;
	         nbins = 20,
	        xlim   = (0.0, 1.0),
	         xlabel = "confidence",
	        ylabel = "Frequency",
	         title=" confidence ")
	h_tl, p_tl = jn.Petit.step_hist(xtl;
	         nbins = 20,
	        #xlim   = (0.0, 1.0),
	         xlabel = " Track Length",
	        ylabel = "Frequency",
	         title=" Track Length ")
	h_ekl, p_ek = jn.Petit.step_hist(xe;
	         nbins = 20,
	        #xlim   = (0.0, 1.0),
	         xlabel = " E (keV)",
	        ylabel = "Frequency",
	         title=" E (keV)")
	plot(p_c, p_tl, p_ek)
end

# ╔═╡ 32f432a9-68c1-4a65-9b98-8f0174bf08ae
begin
	h_b1, p_b1 = jn.Petit.step_hist(xeb1;
	         nbins = 20,
	        #xlim   = (0.0, 1.0),
	         xlabel = "Eb1",
	        ylabel = "Frequency",
	         title=" Eb1 ")
	h_b2, p_b2 = jn.Petit.step_hist(xeb2;
	         nbins = 20,
	        #xlim   = (0.0, 1.0),
	         xlabel = "Eb2",
	        ylabel = "Frequency",
	         title=" Eb2 ")
	
	plot(p_b1, p_b2)
end

# ╔═╡ 93ab0564-f931-4121-992b-156de526c0b1
begin
	
	scatter(xeb1, xeb2)
end

# ╔═╡ b9995760-1dc7-47dd-9ec6-098f8a1d204c
md"roi low (keV) $(@bind roi_low NumberField(2475.0:2500.0, default=2475.0))" 

# ╔═╡ f64b47c1-0aba-4e9a-96a4-3938fb08279b
md"roi up (keV) $(@bind roi_up NumberField(2475.0:2500.0, default=2500.0))" 

# ╔═╡ d0f9e26b-2a3c-42a4-a826-7c8694f5d470
md"""
# Functions
"""

# ╔═╡ 85d7771b-fc92-4844-9086-0d8c3b41e210
function signal_eff(ehrx, rlow, rup; step=10.0) 
	countx = []
	norm = sum(ehrx.weights)
	for rx in rlow:step:rup
		push!(countx, jn.Petit.counts_in_range(ehrx, rx, rup))
	end
	countx = countx /norm
end

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

# ╔═╡ 769f772d-2695-468f-93a0-7655b5d08f10
function histogram_energies_trks(results)
	
	_, p_t1 = jn.Petit.step_hist(results.single_track.energy;
         nbins = 50,
         xlim   = (0.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E one track ")
	_, p_t2 = jn.Petit.step_hist(results.two_track_primary.energy;
         nbins = 50,
         xlim   = (0.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E two tracks ")
	_, p_t2s = jn.Petit.step_hist(results.two_track_secondary.energy;
         nbins = 50,
         xlim   = (0.0, 250.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E secondaries 2t ")
	_, p_t3s = jn.Petit.step_hist(results.three_track_secondary.energy;
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


# ╔═╡ 70914efe-c290-4f43-8c2a-f521281c977e


# ╔═╡ Cell order:
# ╠═947c237c-9852-40e9-a83f-c23666db90aa
# ╠═04b446d6-f34f-11ed-2565-0b15d65b6781
# ╠═871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
# ╠═349825ff-7ffe-4fa1-ba26-a772041f0323
# ╠═7504d7aa-a780-4956-99a5-08a7f9a462b2
# ╠═c9fc0547-0e73-4629-9909-e59c3d75169d
# ╠═b604fcf3-7355-43f5-9b85-4c9c1ff77cd7
# ╠═9ea7167f-493b-4715-930e-3718a3dc6216
# ╠═5bcffb43-a8f7-448f-86c2-f3261d489bbf
# ╠═35d5b13a-9769-410d-99bc-edd8e8cf15fb
# ╠═490d05a8-a555-48d0-9e86-952e69ecaf8e
# ╠═ff6ef332-9b58-4dda-a997-489544d5636e
# ╠═3984ae1a-f844-434d-841d-9a0f0557f85f
# ╠═474b25f8-bd95-4389-867a-bb753dc77d45
# ╠═6c58a746-da45-4e39-a178-91e060b2f34b
# ╠═aa71c155-9bed-4ad6-963e-575100094a9c
# ╠═2416d5e3-8563-45c8-aafd-2f6528c50b85
# ╠═ea1afc2f-9467-4b8b-a3e3-237925ac9fa1
# ╠═5248add6-0d0c-446d-8490-0ccc009fd6e9
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
# ╠═38cb58a7-a752-4d9e-bda9-afa14a9bcf5c
# ╠═149507ef-73e4-43a9-8623-8b3d65bef952
# ╠═f1b3c4ef-1b0c-4243-be00-1b0a15b6f3a2
# ╠═8af4e9c9-2cdc-417d-b681-2647ba807b86
# ╠═c0d93b0b-43f9-4d23-b180-ba813ed7414d
# ╠═4d0c3fe4-5e76-441d-87c7-3a840791257c
# ╠═88aec10d-6110-4408-9037-29007fc6c5d1
# ╠═9bae1645-1485-411d-a6d3-61906c7c4194
# ╠═03cd5a27-1ace-43b3-8a31-4e34629b082d
# ╠═6137bec5-66e7-431f-ab9b-a9cfdf5762c6
# ╠═bf77b535-485a-4608-aa4b-1c15950f3f9b
# ╠═1c202304-d6ca-4dde-a7ac-f5cb2b961b93
# ╠═6c57e536-6815-437a-aba6-64a47a499666
# ╠═98033f21-5720-4d22-94a8-a72e650a860b
# ╠═403b3c7f-d7c3-4762-8605-5d24fb087415
# ╠═b59257e9-9a5e-467f-bbc3-f04cdf0e9a75
# ╠═32f432a9-68c1-4a65-9b98-8f0174bf08ae
# ╠═93ab0564-f931-4121-992b-156de526c0b1
# ╠═b9995760-1dc7-47dd-9ec6-098f8a1d204c
# ╠═f64b47c1-0aba-4e9a-96a4-3938fb08279b
# ╠═d0f9e26b-2a3c-42a4-a826-7c8694f5d470
# ╠═85d7771b-fc92-4844-9086-0d8c3b41e210
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
# ╠═59c673be-1ead-488c-84bd-21bf14980313
# ╠═70914efe-c290-4f43-8c2a-f521281c977e
