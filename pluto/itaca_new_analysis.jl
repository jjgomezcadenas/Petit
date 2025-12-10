### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ f8216615-54b7-4f4a-b20c-65681052a16e
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

# ╔═╡ 5a4d088a-8517-4dd2-b0d4-86b86393b9bc
include("plot_analysis.jl")

# ╔═╡ 63889de2-cfc7-11f0-9979-e5f52b76e205
begin
	cmdir=joinpath(ENV["DATA"], "HD5t/itaca")
	pdir =joinpath(ENV["PROJECTS"], "Petit")
end

# ╔═╡ 6099909b-2c73-4b13-8429-979ac06dcef3
begin
      using Pkg
      Pkg.activate(pdir)
      Pkg.instantiate()
  end

# ╔═╡ 374c62af-72c5-4ca5-98e3-bc17ae89561b
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

# ╔═╡ 30c8283e-ef8e-42a4-9a41-422a139c0770
begin
	jn = ingredients(string(pdir,"/src/Petit.jl"))
	it = ingredients(string(pdir,"/src/itaca_aux.jl"))
	#include(string(pdir,"/src/itaca_aux.jl")) 
end

# ╔═╡ 9ba4863f-d64d-41d7-8e72-1c7bf1c3569b
#Check the diff against main, and remove all AI generated slop introduced in this branch.

#This includes:
#- Extra comments that a human wouldn't add or is inconsistent with the rest of the file
#- Extra defensive checks or try/catch blocks that are abnormal for that area of the codebase (especially if called by trusted / validated codepaths)
#- Casts to any to get around type issues
#- Any other style that is inconsistent with the file

#Report at the end with only a 1-3 sentence summary of what you changed

# ╔═╡ 49ee4ee5-f4e4-4d06-9658-e8961a261702
md"""
# Functions
"""

# ╔═╡ d0c1c7c8-754b-4966-bb70-db5b8387f058
function get_sigma(particle_type, ldrft; 
				   dt = 3.5, dl = 0.9, 
				   tK = 297.0, edrift = 500.0, Pbar=15.0)
	
	if particle_type == "ion"
		σt =  jn.Petit.sigma_t_ion_mm(tK, ldrft, edrift)
		σl =  0.0
	else
		σt =  jn.Petit.sigma_t_mm(ldrft, Pbar; dtmm=dt)
		σl =  jn.Petit.sigma_l_mm(ldrft, Pbar; dlmm=dl)
	end
	σt, σl
end

# ╔═╡ 2afa3fc9-13a5-4538-b33f-0f6324f48876
function get_energy_threshold(particle_type; 
							  energy_threshold_ions =  10.0,
							  energy_threshold_keV =  10.0)
	f = 1e+5/2.5 # ions per MeV
	fkeV = f*1e-3 # ions per keV

	if particle_type == "ion"
		energy_threshold_keV = energy_threshold_ions/fkeV
	end
	energy_threshold_keV
end

# ╔═╡ 4f64c34f-bd1d-45df-8a86-a305ea89eb4c
function get_voxel_size_and_distance(ldrft, σt)
	if ldrft > 50.0
		voxel_scale = 1.5 
		voxel_distance_scale = 1.5
	else
		voxel_scale = 3.0 
		voxel_distance_scale = 2.0
	end
	
	voxel_size = σt * voxel_scale
	mcvox_size = 0.5
	max_distance = voxel_size * voxel_distance_scale
	(voxel_size, mcvox_size, max_distance)
end

# ╔═╡ 89d39d4f-2f3e-450f-a733-341c5dd22c2d
function length_and_energy_of_tracks(tracks)
	LT = [length(track.voxels.energy) for track in tracks]
	E = [sum(track.voxels.energy)*1e+3 for track in tracks]
	LT, E
end

# ╔═╡ 9ed059a1-b269-40f9-a0dd-520579b6a9b8
md"""
# Analysis
"""

# ╔═╡ 941446d9-c2f8-45d7-b998-3ebba4a9812e
md"""
## Process a single event
- Process a single event: diffuse, voxelize, make tracks, project onto path, compute KDE.
- Computes KDE energy profiles for both RECO voxels and MC path.

"""

# ╔═╡ dc29ae3b-6725-413c-be12-228cb548f42d
begin
	particle_type = "electron"
	ldrft = 20.0 # in cm
	
	σt,  σl=  get_sigma(particle_type, ldrft; 
				        dt = 3.5, dl = 0.9, 
				        tK = 297.0, edrift = 500.0, Pbar=15.0)

	voxel_size, mcvox_size, max_distance = get_voxel_size_and_distance(ldrft, σt)

	energy_threshold_keV  = get_energy_threshold(particle_type; 
							                     energy_threshold_ions =  10.0,
							                     energy_threshold_keV =  10.0)

	
	nbins = 100
	nsigma = 3.0 
	kde_bandwidth =2*voxel_size
    mc_kde_bandwidth=2*voxel_size
    n_kde_eval = 200
	
	dfpars = jn.Petit.DiffusionParams(ldrft, σt, σl, voxel_size, max_distance, 
									  energy_threshold_keV, nbins, nsigma)


	

end

# ╔═╡ 2390103c-c761-456b-842e-eeb31d9a6c0c
it.diffusion_params_md(dfpars)

# ╔═╡ aa356ade-9c1a-4e73-8739-212379c543f9
begin
	#bbdf = load_data("bb0nu/bb0nu_15bar_p1.h5")	
	#input_file="bb0nu/bb0nu_15bar_p1.h5"
	input_file="xe137/electrons_2400_2500_15bar_100mum.next.h5"
  	hitsdf =jn.Petit.load_data(input_file, cmdir)
end


# ╔═╡ 42d4cb93-1d37-4f9d-8435-28e70de00485
begin
	nevent=5
	md"- Examine event $(nevent)"
end

# ╔═╡ 939cbc40-985c-4963-967f-3eb1f1807c1d
md""" ### Get event hits
"""

# ╔═╡ b748e2a3-0504-4b2f-9d80-365081dc21f6
begin
	event_df = jn.Petit.get_event(hitsdf, nevent)
	jn.Petit.plot_event(event_df)
end

# ╔═╡ 64e87873-f760-4d8f-8837-427231af7951
md"""### Compute MC path 
- MC path (voxelized primary particle trajectory with arc-length)
- First/last rows give MC truth extremes
"""
   

# ╔═╡ 5746835d-4a6a-4680-b97e-5131c0649179
begin
	mc_path = jn.Petit.compute_mc_path(event_df, mcvox_size)
	jn.Petit.plot_event(event_df; mc_path=mc_path)
end

# ╔═╡ 7bd5586b-f465-43b9-bdbd-d221eb6545f5
md"### Transform hits (removes time, particle_id, etc.)"

# ╔═╡ 4aa91af9-cebd-471b-9ecf-786f1ec11690
event_mc = jn.Petit.transform_hits_df(event_df)

# ╔═╡ f93100e5-37ae-4b43-9a8b-c5a5fa00e71f
md"### Diffuse event"

# ╔═╡ c92bce64-39b7-4282-8291-9805112a554f
begin
	diffused_df = jn.Petit.diffuse_xyz_image_mc(event_mc;
                                             sigma_t_mm=σt,
                                             sigma_l_mm=σl,
                                             nbins=nbins,
                                             nsigma=nsigma)
	histogram(diffused_df.electrons)
end

# ╔═╡ 2150e65f-b949-4874-b308-12499c020273
jn.Petit.plot_hits(diffused_df, energy_column=:electrons)


# ╔═╡ 3c3272a8-0368-413e-a4f6-7e233bbcc699
md"### Voxelize event"

# ╔═╡ a831ea92-c29b-4e21-96e0-704ee82f7795
begin
	voxels = jn.Petit.voxelize_event(diffused_df, voxel_size)
	jn.Petit.plot_event(voxels)
end

# ╔═╡ 8d7a02f4-12d2-4501-80fc-52f4ec930703
    md"### Make MC tracks"
    

# ╔═╡ 725637cf-507c-4c39-8bf4-6d55233b1af6
md"#### Voxelize MC "

# ╔═╡ 9e5c887b-bd64-4b12-8231-531499323288
begin
	mcvoxels = jn.Petit.voxelize_event(event_mc, mcvox_size )
	jn.Petit.plot_event(mcvoxels)
end

# ╔═╡ 20922702-cf5f-4bf5-8be9-1a3c96f11e30
begin
mctracks = jn.Petit.make_tracks(mcvoxels;
                               max_distance_mm=2.5 * mcvox_size,
                               energy_threshold_kev=1.0,
                               diffusion=dfpars)
	
md"""
- number of MC tracks found = $(length(mctracks))
"""
end

# ╔═╡ 4afafdfd-b783-4fcd-9942-dcb08dd77549
if length(mctracks) > 1
	LT, E = length_and_energy_of_tracks(mctracks)
	println("# of voxels: ", join([l for l in LT], ", "),)
	println("Energies: ", join([round(e, digits=2) for e in E], ", "), " keV")
end

# ╔═╡ 33d6f97a-ccdf-4b65-8133-3704f5fea584
 md"### Make RECO tracks"

# ╔═╡ 83891e72-df80-4b9d-9989-7768604f2f5a
begin
tracks = jn.Petit.make_tracks(voxels;
                               max_distance_mm=max_distance,
                               energy_threshold_kev=energy_threshold_keV,
                               diffusion=dfpars)
md"""
- number of tracks found = $(length(tracks))
"""
end

# ╔═╡ e28d9f0e-5a7e-46e6-acee-02c8135a88ce
begin
	track = tracks[1]
	walk_result = jn.Petit.walk_track_from_extremes(track)
	it.walk_result_md(walk_result)
end

# ╔═╡ 41024877-b6e8-49f7-8301-5490f276ca85
begin
	path = jn.Petit.get_raw_path(track, walk_result.path_indices)
	it.path_md(path::DataFrame)
end

# ╔═╡ dbdfc123-6bc6-458b-93f8-603d589646d1
begin
	extreme_dists = jn.Petit.compute_extreme_distances(path, mc_path)
	it.extreme_distances_md(extreme_dists)
end

# ╔═╡ 6460e445-f8c0-48c3-bef9-1cefc042ecad
jn.Petit.plot_track_with_paths(track, path, mc_path;
                              show_distances=false)

# ╔═╡ 02ba72f8-6437-45b4-a6cc-b2554a2317c9
md"""
### KDE
"""

# ╔═╡ e4f9a433-fbea-4068-ae2b-66556365971f
begin
	 reco_kde = jn.Petit.get_reco_kde(track, path; bandwidth=kde_bandwidth, 
									  n_eval=n_kde_eval)
	 mc_kde = jn.Petit.get_mc_kde(mc_path; bandwidth=mc_kde_bandwidth, 
								  n_eval=n_kde_eval)
end

# ╔═╡ e55bc8f0-c003-4f2a-b35c-9705fc53fe02


# ╔═╡ 1a94a0d8-8027-4067-8697-11f88b4ecf14
let
    p1 = plot(reco_kde.kde_s, reco_kde.kde_f,
              xlabel="Arc length s (mm)",
              ylabel=" KDE Energy density f(s)",
              title="Longitudinal Energy Density (Event $(nevent))",
              label="RECO KDE (h=$(round(kde_bandwidth, digits=1))mm)",
              linewidth=2, 
              color=:blue)
	p2 = plot!(p1, mc_kde.kde_s, mc_kde.kde_f,
              xlabel="Arc length s (mm)",
              ylabel=" KDE Energy density f(s)",
              title="Longitudinal Energy Density (Event $(nevent))",
              label="MC KDE (h=$(round(kde_bandwidth, digits=1))mm)",
              linewidth=2, 
              color=:red)
	p3 = histogram(reco_kde.s, weights=reco_kde.E .* 1e3,
                   bins=range(0, reco_kde.track_length, length=n_kde_eval+1),
                   xlabel="Arc length s (mm)",
                   ylabel="Energy (keV)",
                   title="Energy vs S (Event $(nevent)",
                   label="RECO",
                   fillalpha=0.7,
                   color=:gray)
	p4 = histogram!(p3, mc_kde.s, weights=mc_kde.E .* 1e3,
                   bins=range(0, mc_kde.track_length, length=n_kde_eval+1),
                   xlabel="Arc length s (mm)",
                   ylabel="Energy (keV)",
                   title="Energy vs S (MC)",
                   label="MC",
                   fillalpha=0.7,
                   color=:green)
	plot(p2,p4, layout=(1, 2), size=(1200, 1000))
	
end

# ╔═╡ d5bda625-3ab9-4eb8-83d7-9192888abe08
md"### Find peaks " 

# ╔═╡ 9c237e6d-e798-4b61-b7d8-21f2fd8f5a6e
begin
	peaks = jn.Petit.find_peaks(reco_kde.kde_f, reco_kde.kde_s; prom_scale=0.2)
	jn.Petit.plot_kde_peaks(reco_kde.kde_s, reco_kde.kde_f, peaks; title="RECO KDE")
end

# ╔═╡ 638ae541-9ba1-4dd8-83b3-a474a13d245e
peaks

# ╔═╡ 30743fc8-632c-436d-9cd9-06a254e59015
minf = minimum(reco_kde.kde_f)

# ╔═╡ e27b3ee8-0a47-4baa-aa8d-9098ac802c5d
 pk = it.kde_peaks(peaks, reco_kde.kde_f)

# ╔═╡ 791ae0ee-1b30-4349-938d-3539a196b549
md"""
### Blob analysis
"""

# ╔═╡ 03e5c952-1e75-4a2d-bee5-1f262d0ae38c
Rb = 30.0

# ╔═╡ e97c03c9-aaeb-497b-bd5d-8ccd63f889ff
blobs = jn.Petit.find_blob_energies(track,
                                    path;
                                    radius=Rb) 

# ╔═╡ 632673fc-2cf5-4c2a-8d4b-f2a634f311f4
 jn.Petit.plot_reco_track_with_voxels_and_spheres(track, path, blobs, Rb)

# ╔═╡ 9ba38a7d-539c-4599-9a98-fd9070408b35
md"""
## End of main analysis
"""

# ╔═╡ 7966be51-97ca-4475-a8ed-592bd0fae39a
md"""
## Offline analysis
"""

# ╔═╡ 8c9d1192-f38f-43a9-ad69-88881be55a95
function select_kde(bbdf, xedf; pmin=0.5, lmax = 50, rmin=150)
	

	# --- Define a helper predicate for signal-like events ---
	function is_signal_like(ev; p_min, L_max, R_min)
	    p1 = ev.peak1_prom
	    p2 = ev.peak2_prom
	    L1 = ev.peak1_left
	    R2 = ev.peak2_right

    	left_ok  = (L1 ≤ lmax)
    	right_ok = (p2 ≥ p_min) && (R2 ≥ rmin)

    	return left_ok && right_ok
	end

	# --- Apply to BB0ν (signal) ---
	bbsel = filter(ev -> is_signal_like(ev; p_min=pmin,
										 L_max=lmax,
										 R_min=rmin),
					                     bbdf)
	xesel = filter(ev -> is_signal_like(ev; p_min=pmin,
										 L_max=lmax,
										 R_min=rmin),
				   						 xedf)
	return bbsel, xesel
end

# ╔═╡ 17be29e6-9b6d-4934-b5fe-0265eaafd2b2
function kde_eff(bbdf, xedf)
	EFBB = []
	EFXE = []
	range = 0.2:0.1:1.0
	for p in 0.2:0.1:1.0
		bbsel, xesel =select_kde(bbdf, xedf; pmin=p, lmax = 50, rmin=150)
		push!(EFBB, size(bbsel)[1]/size(bbdf)[1])
		push!(EFXE, size(xesel)[1]/size(xedf)[1])
	end
	collect(range), EFBB, EFXE
	
end

# ╔═╡ 8c03d148-60d0-4a74-a5ab-c70f86001b82
	function select_blob(bb, xe; eb2cut=400.0)
		bbsel = filter(row -> row.Eb2_keV > eb2cut, bb )
		xesel = filter(row -> row.Eb2_keV > eb2cut, xe )
		return bbsel, xesel
	end

# ╔═╡ ca4db9f1-8c14-4ade-ac59-3ed939f699cd
begin
	xe_ion_ldrft_100 =  load_analysis_results("/Users/jjgomezcadenas/Projects/Petit/ItacaResultsMT/xe_p15_l100_ion_analysis.csv")
	xe_ion_ldrft_100_meta =  load_analysis_results("/Users/jjgomezcadenas/Projects/Petit/ItacaResultsMT/xe_p15_l100_ion_metadata.csv")
end

# ╔═╡ 802bbbbe-f4ff-49ec-93eb-174696b7516b
begin 
	 plot_all_analysis(xe_ion_ldrft_100)
end

# ╔═╡ 86d124a2-be19-43f0-b1de-5a914c08cbf0
begin
	bb_ion_ldrft_100 =  load_analysis_results("/Users/jjgomezcadenas/Projects/Petit/ItacaResultsMT/bb_p15_l100_ion_analysis.csv")
	bb_ion_ldrft_100_meta =  load_analysis_results("/Users/jjgomezcadenas/Projects/Petit/ItacaResultsMT/bb_p15_l100_ion_metadata.csv")
end

# ╔═╡ 1886f26c-0b24-48ff-adac-36849bc3bd01
plot_all_analysis(bb_ion_ldrft_100)

# ╔═╡ 4665cf3d-f0c7-4e89-aca7-cd3f18e574d6
scatter(bb_ion_ldrft_100.peak1_left, bb_ion_ldrft_100.peak2_right;
            xlabel="peak1 left",
            ylabel="peak2 right",
            title="bb0nu peak1 vs peak2 pos",
            legend=false,
            markersize=2,
            alpha=0.6)

# ╔═╡ dc2e7af8-78aa-46b3-b116-a5e22661ec74
scatter(xe_ion_ldrft_100.peak1_left, xe_ion_ldrft_100.peak2_right;
            xlabel="peak1 left",
            ylabel="peak2 right",
            title="xe137 peak1 vs peak2 pos",
            legend=false,
            markersize=2,
            alpha=0.6)

# ╔═╡ 8325e8eb-9fe3-4d14-9997-ef1739b843f7
let
	sbb = scatter(bb_ion_ldrft_100.peak1_left, bb_ion_ldrft_100.peak1_prom;
            xlabel="peak1 left",
            ylabel="peak1 prom",
            title="bb0nu peak1 vs prom1",
            legend=false,
            markersize=2,
            alpha=0.6,  ylims=(0, 1.5))
	sxe = scatter(xe_ion_ldrft_100.peak1_left, xe_ion_ldrft_100.peak1_prom;
            xlabel="peak1 left",
            ylabel="peak1 prom",
            title="xe137 peak1 vs prom1",
            legend=false,
            markersize=2,
            alpha=0.6, ylims=(0, 1.5))
	plot(sbb, sxe)
end

# ╔═╡ 570754ec-4321-4278-95b4-105d548e9743
let
	sbb = scatter(bb_ion_ldrft_100.peak2_right, bb_ion_ldrft_100.peak2_prom;
            xlabel="peak2 right",
            ylabel="peak2 prom",
            title="bb0nu peak2 vs prom2",
            legend=false,
            markersize=2,
            alpha=0.6, ylims=(0, 1.5))
	sxe = scatter(xe_ion_ldrft_100.peak2_right, xe_ion_ldrft_100.peak2_prom;
            xlabel="peak2 right",
            ylabel="peak2 prom",
            title="xe137 peak2 vs prom2",
            legend=false,
            markersize=2,
            alpha=0.6, ylims=(0, 1.5))
	plot(sbb, sxe)
end

# ╔═╡ de04dfd0-d1f8-4b4b-a28c-e0fa97419190
let
	cuts, effbb, effxe = kde_eff(bb_ion_ldrft_100, xe_ion_ldrft_100)
	p = plot(cuts, effbb,
             label="bb0nu",
             color=:blue,
             linewidth=2,
             xlabel="prom",
             ylabel="Efficiency",
             title="Efficiency KDE",
             legend=:right,
             grid=true)

    p2 = plot!(p, cuts, effxe,
          label="Xe137",
          color=:red,
          linewidth=2)
	p3 = plot(cuts, effbb./sqrt.(effxe),
             label="fomd",
             color=:blue,
             linewidth=2,
             xlabel="prom",
             ylabel="Efficiency",
             title="fom KDE",
             legend=:right,
             grid=true)
	plot(p2, p3)
end

# ╔═╡ 7fd09206-2652-4f0d-adb8-85c4b078a2c9
begin
bb_ion_ldrft_100_pp, xe_ion_ldrft_100_pp =select_kde(bb_ion_ldrft_100, 
													   xe_ion_ldrft_100; 
													   pmin=0.5, lmax = 50, rmin=150)
	eff_kde_bb = size(bb_ion_ldrft_100_pp)[1]/size(bb_ion_ldrft_100)[1]
	eff_kde_xe = size(xe_ion_ldrft_100_pp)[1]/size(xe_ion_ldrft_100)[1]
	md"""
	- KDE eff for bb = $(eff_kde_bb)
	- KDE eff for Xe = $(eff_kde_xe)
	"""
end

# ╔═╡ 4a34f718-432f-4c74-97d3-22a866c0ceef
size(bb_ion_ldrft_100_pp)[1]/size(bb_ion_ldrft_100)[1]

# ╔═╡ ce7d355b-5866-40e4-910b-155a4f9587a8
begin
	
	sbb1 = scatter(bb_ion_ldrft_100_pp.peak1_prom, bb_ion_ldrft_100_pp.peak2_prom;
	            xlabel="peak1 prom",
	            ylabel="peak2 prom",
	            title="bb0nu peak1 vs peak2 proms",
	            legend=false,
	            markersize=2,
	            alpha=0.6, ylims=(0, 1.5))
	sxe1 = scatter(xe_ion_ldrft_100_pp.peak1_prom, xe_ion_ldrft_100_pp.peak2_prom;
	            xlabel="peak1 prom",
	            ylabel="peak2 prom",
	            title="xe137 peak1 vs peak2 proms",
	            legend=false,
	            markersize=2,
	            alpha=0.6, ylims=(0, 1.5))
	plot(sbb1, sxe1)
end

# ╔═╡ ecd2c6e3-816f-447d-a50e-8c7416bb1c8b
let
	bbpeb = plot_eb1_vs_eb2(bb_ion_ldrft_100)
	xepeb = plot_eb1_vs_eb2(xe_ion_ldrft_100)
	bbpeb2 = plot_eb1_vs_eb2(bb_ion_ldrft_100_pp)
	xepeb2 = plot_eb1_vs_eb2(xe_ion_ldrft_100_pp)
	plot(bbpeb, xepeb,bbpeb2,xepeb2)
end

# ╔═╡ cbc2014e-4803-49c2-a46a-810bfdabae49
let
	eb2_sig = bb_ion_ldrft_100.Eb2_keV
	eb2_bkg = xe_ion_ldrft_100.Eb2_keV
	p1 = jn.Petit.plot_efficiency_vs_cut(eb2_sig,
                                eb2_bkg;
                                cuts = range(200, 600, length=40),
                                xlabel = "Cut Value",
                                title = "Efficiency vs Cut",
                                signal_label = "bb0nu",
                                background_label = "Xe137")
	p2 = jn.Petit.plot_fom_vs_cut(eb2_sig,
                         eb2_bkg;
                         cuts = range(200, 600, length=40),
                         xlabel = "Cut Value",
                         title= "FOM vs Cut",
                         mark_optimal = true)
	plot(p1, p2)

end

# ╔═╡ e59e67e9-3476-44a9-b233-e88bfed26a00
begin
	bb_ion_ldrft_100_bl, xe_ion_ldrft_100_bl =select_blob(bb_ion_ldrft_100, 
													  xe_ion_ldrft_100; 
													  eb2cut=400.0)
	eff_bl_bb = size(bb_ion_ldrft_100_bl)[1]/size(bb_ion_ldrft_100)[1]
	eff_bl_xe = size(xe_ion_ldrft_100_bl)[1]/size(xe_ion_ldrft_100)[1]
	md"""
	- BLOB eff for bb = $(eff_bl_bb)
	- BLOB eff for Xe = $(eff_bl_xe)
	"""
end

# ╔═╡ 3927b2da-867e-49e4-b755-ecda0f0904f3
let
	eb2_sig = bb_ion_ldrft_100_pp.Eb2_keV
	eb2_bkg = xe_ion_ldrft_100_pp.Eb2_keV
	p1 = jn.Petit.plot_efficiency_vs_cut(eb2_sig,
                                eb2_bkg;
                                cuts = range(200, 600, length=40),
                                xlabel = "Cut Value",
                                title = "Efficiency vs Cut",
                                signal_label = "bb0nu",
                                background_label = "Xe137")
	p2 = jn.Petit.plot_fom_vs_cut(eb2_sig,
                         eb2_bkg;
                         cuts = range(200, 600, length=40),
                         xlabel = "Cut Value",
                         title= "FOM vs Cut",
                         mark_optimal = true)
	plot(p1, p2)

end

# ╔═╡ 1800b4a7-2a65-4da0-a6d3-b04196dbc334
let
	
	bb_ion_ldrft_100_bl, xe_ion_ldrft_100_bl =select_blob(bb_ion_ldrft_100_pp, 
													  xe_ion_ldrft_100_pp; 
													  eb2cut=350.0)
	eff_bl_bb = size(bb_ion_ldrft_100_bl)[1]/size(bb_ion_ldrft_100)[1]
	eff_bl_xe = size(xe_ion_ldrft_100_bl)[1]/size(xe_ion_ldrft_100)[1]
	eff_tot_bb = eff_kde_bb * eff_bl_bb
	eff_tot_xe =  eff_kde_xe * eff_bl_xe
	md"""
	- KDE eff for bb = $(eff_kde_bb)
	- KDE eff for Xe = $(eff_kde_xe)
	- BLOB eff for bb = $(eff_bl_bb)
	- BLOB eff for Xe = $(eff_bl_xe)
	- TOT eff for bb = $(eff_tot_bb)
	- TOT eff for Xe = $(eff_tot_xe)
	"""
end

# ╔═╡ 88b09dfd-dfb6-4e91-a29e-ed9e46e00ff0
jn.Petit.roc_curve(bb_ion_ldrft_100.Eb2_keV,
                  xe_ion_ldrft_100.Eb2_keV;
                   cuts = range(0, 600.0, length=50))

# ╔═╡ 03b93cba-bd2a-42a6-ae12-204e9f1c6772
prion = jn.Petit.plot_roc(bb_ion_ldrft_100.Eb2_keV,
                  xe_ion_ldrft_100.Eb2_keV;
                  cuts = range(0, 800.0, length=500),
                  title = "ROC Curve",
                  show_auc = true,
                  show_diagonal = false)

# ╔═╡ 86353163-59db-44e2-b6dd-6d0c378d57ef
md"""
### Electrons
"""

# ╔═╡ 23d2b5c1-2821-430b-acf2-e4f78145c5f3
begin
	bb_elec_ldrft_100 =  load_analysis_results("/Users/jjgomezcadenas/Projects/Petit/ItacaResultsMT/bb_p15_l100_elec_analysis.csv")
	bb_elec_ldrft_100_meta =  load_analysis_results("/Users/jjgomezcadenas/Projects/Petit/ItacaResultsMT/bb_p15_l100_elec_metadata.csv")
end

# ╔═╡ 6d0b6cd0-9e46-4fbd-be82-624ab083b89c
plot_all_analysis(bb_elec_ldrft_100)

# ╔═╡ a6a05b47-f857-45ae-95d2-e0198a0ac53d
begin
	xe_elec_ldrft_100 =  load_analysis_results("/Users/jjgomezcadenas/Projects/Petit/ItacaResultsMT/xe_p15_l100_elec_analysis.csv")
	xe_elec_ldrft_100_meta =  load_analysis_results("/Users/jjgomezcadenas/Projects/Petit/ItacaResultsMT/xe_p15_l100_elec_metadata.csv")
end

# ╔═╡ c4c0bb36-d255-4a35-88df-c69118babf54
plot_all_analysis(xe_elec_ldrft_100)

# ╔═╡ 8dcc809e-abc5-4e7d-a2d4-2806ad4428a7
let
	sbb = scatter(bb_elec_ldrft_100.peak1_left, bb_ion_ldrft_100.peak1_prom;
            xlabel="peak1 left",
            ylabel="peak1 prom",
            title="bb0nu peak1 vs prom1",
            legend=false,
            markersize=2,
            alpha=0.6,  ylims=(0, 1.5))
	sxe = scatter(xe_elec_ldrft_100.peak1_left, xe_ion_ldrft_100.peak1_prom;
            xlabel="peak1 left",
            ylabel="peak1 prom",
            title="xe137 peak1 vs prom1",
            legend=false,
            markersize=2,
            alpha=0.6, ylims=(0, 1.5))
	plot(sbb, sxe)
end

# ╔═╡ 19c207b5-59e0-43e7-b385-985d1510b448
let
	sbb = scatter(bb_elec_ldrft_100.peak2_right, bb_ion_ldrft_100.peak2_prom;
            xlabel="peak2 right",
            ylabel="peak2 prom",
            title="bb0nu peak2 vs prom2",
            legend=false,
            markersize=2,
            alpha=0.6, ylims=(0, 1.5))
	sxe = scatter(xe_elec_ldrft_100.peak2_right, xe_ion_ldrft_100.peak2_prom;
            xlabel="peak2 right",
            ylabel="peak2 prom",
            title="xe137 peak2 vs prom2",
            legend=false,
            markersize=2,
            alpha=0.6, ylims=(0, 1.5))
	plot(sbb, sxe)
end

# ╔═╡ 935ac741-dd8f-4566-a9ca-8c85e83fe284
let
	eb2_sig = bb_elec_ldrft_100.Eb2_keV
	eb2_bkg = xe_elec_ldrft_100.Eb2_keV
	p1 = jn.Petit.plot_efficiency_vs_cut(eb2_sig,
                                eb2_bkg;
                                cuts = range(200, 600, length=40),
                                xlabel = "Cut Value",
                                title = "Efficiency vs Cut",
                                signal_label = "bb0nu",
                                background_label = "Xe137")
	p2 = jn.Petit.plot_fom_vs_cut(eb2_sig,
                         eb2_bkg;
                         cuts = range(200, 600, length=40),
                         xlabel = "Cut Value",
                         title= "FOM vs Cut",
                         mark_optimal = true)
	plot(p1, p2)

end

# ╔═╡ ccc8aa47-0795-4f9b-8719-f4de220302ee
begin
	bb_elec_ldrft_100_bl, xe_elec_ldrft_100_bl =select_blob(bb_elec_ldrft_100,   
												            xe_elec_ldrft_100; 
													          eb2cut=400.0)
	eff_blb_bb = size(bb_elec_ldrft_100_bl)[1]/size(bb_elec_ldrft_100)[1]
	eff_blb_xe = size(xe_elec_ldrft_100_bl)[1]/size(xe_elec_ldrft_100)[1]
	md"""
	- BLOB eff for bb = $(eff_blb_bb)
	- BLOB eff for Xe = $(eff_blb_xe)
	"""
end

# ╔═╡ d6ece5e0-287f-4e3a-b7f3-1138a60eea3e
begin
prele = jn.Petit.plot_roc(bb_elec_ldrft_100.Eb2_keV,
                  xe_elec_ldrft_100.Eb2_keV;
                  cuts = range(0, 800.0, length=500),
                  title = "ROC Curve",
                  show_auc = true,
                  show_diagonal = false)
plot(prion, prele)
end

# ╔═╡ 7a621f23-28bb-4710-879f-43af85b2fd30


# ╔═╡ Cell order:
# ╠═63889de2-cfc7-11f0-9979-e5f52b76e205
# ╠═6099909b-2c73-4b13-8429-979ac06dcef3
# ╠═f8216615-54b7-4f4a-b20c-65681052a16e
# ╠═374c62af-72c5-4ca5-98e3-bc17ae89561b
# ╠═30c8283e-ef8e-42a4-9a41-422a139c0770
# ╠═9ba4863f-d64d-41d7-8e72-1c7bf1c3569b
# ╠═49ee4ee5-f4e4-4d06-9658-e8961a261702
# ╠═d0c1c7c8-754b-4966-bb70-db5b8387f058
# ╠═2afa3fc9-13a5-4538-b33f-0f6324f48876
# ╠═4f64c34f-bd1d-45df-8a86-a305ea89eb4c
# ╠═89d39d4f-2f3e-450f-a733-341c5dd22c2d
# ╠═9ed059a1-b269-40f9-a0dd-520579b6a9b8
# ╠═941446d9-c2f8-45d7-b998-3ebba4a9812e
# ╠═dc29ae3b-6725-413c-be12-228cb548f42d
# ╠═2390103c-c761-456b-842e-eeb31d9a6c0c
# ╠═aa356ade-9c1a-4e73-8739-212379c543f9
# ╠═42d4cb93-1d37-4f9d-8435-28e70de00485
# ╠═939cbc40-985c-4963-967f-3eb1f1807c1d
# ╠═b748e2a3-0504-4b2f-9d80-365081dc21f6
# ╠═64e87873-f760-4d8f-8837-427231af7951
# ╠═5746835d-4a6a-4680-b97e-5131c0649179
# ╠═7bd5586b-f465-43b9-bdbd-d221eb6545f5
# ╠═4aa91af9-cebd-471b-9ecf-786f1ec11690
# ╠═f93100e5-37ae-4b43-9a8b-c5a5fa00e71f
# ╠═c92bce64-39b7-4282-8291-9805112a554f
# ╠═2150e65f-b949-4874-b308-12499c020273
# ╠═3c3272a8-0368-413e-a4f6-7e233bbcc699
# ╠═a831ea92-c29b-4e21-96e0-704ee82f7795
# ╠═8d7a02f4-12d2-4501-80fc-52f4ec930703
# ╠═725637cf-507c-4c39-8bf4-6d55233b1af6
# ╠═9e5c887b-bd64-4b12-8231-531499323288
# ╠═20922702-cf5f-4bf5-8be9-1a3c96f11e30
# ╠═4afafdfd-b783-4fcd-9942-dcb08dd77549
# ╠═33d6f97a-ccdf-4b65-8133-3704f5fea584
# ╠═83891e72-df80-4b9d-9989-7768604f2f5a
# ╠═e28d9f0e-5a7e-46e6-acee-02c8135a88ce
# ╠═41024877-b6e8-49f7-8301-5490f276ca85
# ╠═dbdfc123-6bc6-458b-93f8-603d589646d1
# ╠═6460e445-f8c0-48c3-bef9-1cefc042ecad
# ╠═02ba72f8-6437-45b4-a6cc-b2554a2317c9
# ╠═e4f9a433-fbea-4068-ae2b-66556365971f
# ╠═e55bc8f0-c003-4f2a-b35c-9705fc53fe02
# ╠═1a94a0d8-8027-4067-8697-11f88b4ecf14
# ╠═d5bda625-3ab9-4eb8-83d7-9192888abe08
# ╠═9c237e6d-e798-4b61-b7d8-21f2fd8f5a6e
# ╠═638ae541-9ba1-4dd8-83b3-a474a13d245e
# ╠═30743fc8-632c-436d-9cd9-06a254e59015
# ╠═e27b3ee8-0a47-4baa-aa8d-9098ac802c5d
# ╠═791ae0ee-1b30-4349-938d-3539a196b549
# ╠═03e5c952-1e75-4a2d-bee5-1f262d0ae38c
# ╠═e97c03c9-aaeb-497b-bd5d-8ccd63f889ff
# ╠═632673fc-2cf5-4c2a-8d4b-f2a634f311f4
# ╠═9ba38a7d-539c-4599-9a98-fd9070408b35
# ╠═7966be51-97ca-4475-a8ed-592bd0fae39a
# ╠═5a4d088a-8517-4dd2-b0d4-86b86393b9bc
# ╠═8c9d1192-f38f-43a9-ad69-88881be55a95
# ╠═17be29e6-9b6d-4934-b5fe-0265eaafd2b2
# ╠═8c03d148-60d0-4a74-a5ab-c70f86001b82
# ╠═ca4db9f1-8c14-4ade-ac59-3ed939f699cd
# ╠═802bbbbe-f4ff-49ec-93eb-174696b7516b
# ╠═86d124a2-be19-43f0-b1de-5a914c08cbf0
# ╠═1886f26c-0b24-48ff-adac-36849bc3bd01
# ╠═4665cf3d-f0c7-4e89-aca7-cd3f18e574d6
# ╠═dc2e7af8-78aa-46b3-b116-a5e22661ec74
# ╠═8325e8eb-9fe3-4d14-9997-ef1739b843f7
# ╠═570754ec-4321-4278-95b4-105d548e9743
# ╠═4a34f718-432f-4c74-97d3-22a866c0ceef
# ╠═de04dfd0-d1f8-4b4b-a28c-e0fa97419190
# ╠═7fd09206-2652-4f0d-adb8-85c4b078a2c9
# ╠═ce7d355b-5866-40e4-910b-155a4f9587a8
# ╠═ecd2c6e3-816f-447d-a50e-8c7416bb1c8b
# ╠═cbc2014e-4803-49c2-a46a-810bfdabae49
# ╠═e59e67e9-3476-44a9-b233-e88bfed26a00
# ╠═3927b2da-867e-49e4-b755-ecda0f0904f3
# ╠═1800b4a7-2a65-4da0-a6d3-b04196dbc334
# ╠═88b09dfd-dfb6-4e91-a29e-ed9e46e00ff0
# ╠═03b93cba-bd2a-42a6-ae12-204e9f1c6772
# ╠═86353163-59db-44e2-b6dd-6d0c378d57ef
# ╠═23d2b5c1-2821-430b-acf2-e4f78145c5f3
# ╠═6d0b6cd0-9e46-4fbd-be82-624ab083b89c
# ╠═a6a05b47-f857-45ae-95d2-e0198a0ac53d
# ╠═c4c0bb36-d255-4a35-88df-c69118babf54
# ╠═8dcc809e-abc5-4e7d-a2d4-2806ad4428a7
# ╠═19c207b5-59e0-43e7-b385-985d1510b448
# ╠═935ac741-dd8f-4566-a9ca-8c85e83fe284
# ╠═ccc8aa47-0795-4f9b-8719-f4de220302ee
# ╠═d6ece5e0-287f-4e3a-b7f3-1138a60eea3e
# ╠═7a621f23-28bb-4710-879f-43af85b2fd30
