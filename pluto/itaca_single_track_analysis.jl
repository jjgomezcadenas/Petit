### A Pluto.jl notebook ###
# v0.20.21

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
	#it = ingredients(string(pdir,"/src/itaca_aux.jl"))
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

# ╔═╡ aa356ade-9c1a-4e73-8739-212379c543f9
function get_data(; datatype="bb0nu")
	nbb = [3,6,9,14,15,20,21,22,26,28]
	nxe = [12,22,86,190,192,292,334,341,348,377]
	if datatype=="bb0nu"
		input_file="bb0nu/selected_events.h5"
		is_double_beta=true
		
		hitsdf =jn.Petit.load_data(input_file, cmdir)
		return hitsdf, is_double_beta, nbb
	elseif datatype=="xe137"
		input_file="xe137/selected_events.h5"
		is_double_beta=false
  		hitsdf =jn.Petit.load_data(input_file, cmdir)
		return hitsdf, is_double_beta, nxe
	else
		println("not yet implemented")
	end
end


# ╔═╡ e549464d-6e19-4d15-a64a-0f62419f80b3
function print_diagnostics(label, track, res, R_cover, coords)
    u, v, path = res[1], res[2], res[3]
    println("\n--- $label ---")
    println("endpoints = ($u, $v)")

    # A) path efficiency
    η, L_path, D_end = jn.Petit.diagnose_path_efficiency(coords, path)
    println("path_efficiency η = L_path/D_end = $(round(η, digits=3))  (L_path=$(round(L_path,digits=2)) mm, D_end=$(round(D_end,digits=2)) mm)")

    # B) endpoint degrees (interior-ness)
    deg_u, deg_v = jn.Petit.diagnose_endpoint_degrees(track.graph, u, v)
    println("deg(u)=$deg_u, deg(v)=$deg_v")

    # C) skeleton coverage (how much energy lies near skeleton)
    f, Ein, Etot = jn.Petit.diagnose_skeleton_coverage(track, coords, path; R_cover=R_cover, energy_col=:energy)
    println("coverage f = Ein/Etot = $(round(f,digits=3))  (Ein=$(round(Ein,digits=3)), Etot=$(round(Etot,digits=3)))")
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
	particle_type = "ion"
	ldrft = 100.0 # in cm
	
	σt,  σl=  jn.Petit.get_sigma(particle_type, ldrft; 
				        dt = 3.5, dl = 0.9, 
				        tK = 297.0, edrift = 500.0, Pbar=15.0)

	voxel_size, mcvox_size, max_distance = jn.Petit.get_voxel_size_and_distance(ldrft, σt)

	energy_threshold_keV  = jn.Petit.get_energy_threshold(particle_type; 
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

# ╔═╡ 57895dc3-45d5-497b-8034-00950ce45fde
function diagonse(track; use_energy_weighting=true,
	                     use_edge_energy_weighting=true,
				         use_mst_fallback=false)
	
	coords = jn.Petit.extract_coords(track)
	# (a) Combined method
	res_comb = jn.Petit.find_extremes_combined(track, coords;
	                             use_energy_weighting=use_energy_weighting,
	                             use_edge_energy_weighting=use_edge_energy_weighting,
	                             use_mst_fallback=use_mst_fallback=false)
	# energy weighted
	res_edge = jn.Petit.find_extremes_edge_energy_weighted_opt(track, coords)
	
	# mst
	res_mst  = jn.Petit.find_extremes_mst_diameter(track, coords)
	
	print_diagnostics("COMBINED", track, res_comb, max_distance, coords)
	#println(res_comb)
	print_diagnostics("EDGE-WEIGHTED", track, res_edge, max_distance, coords)
	#println(res_mst)
	print_diagnostics("MST",  track, res_mst, max_distance, coords)

	# Stability test (best fallback trigger) ------------------------------
	# Measure endpoint displacement between two solutions (allowing endpoint swap).
	Δ_comb_vs_edge = jn.Petit.diagnose_endpoint_stability(res_comb, res_edge, coords)
	Δ_comb_vs_mst  = jn.Petit.diagnose_endpoint_stability(res_comb, res_mst,  coords)

	println("\n--- Stability ---")
	
	println("Δ(combined, edge) = $(round(Δ_comb_vs_edge, digits=2)) mm")
	println("Δ(combined, mst)  = $(round(Δ_comb_vs_mst,  digits=2)) mm")	
end

# ╔═╡ 2390103c-c761-456b-842e-eeb31d9a6c0c
jn.Petit.diffusion_params_md(dfpars)

# ╔═╡ 0cdaf271-c87e-4cb4-aadf-428e4af466bc
function get_walk_mst(res_mst,track)
	extreme1, extreme2, path, confidence = res_mst
	# Get voxel data for extremes
    start_voxel = track.voxels[extreme1, :]
    end_voxel = track.voxels[extreme2, :]

    # Get voxels along the path
    path_voxels = track.voxels[path, :]

	# Calculate total path length
    total_length = 0.0
    for i in 1:length(path)-1
        v1, v2 = path[i], path[i+1]
        total_length += jn.Petit.euclidean_distance(
            track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
            track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
        )
    end
	
	return (extremes = (start_voxel, end_voxel),
            path_indices = path,
            path_voxels = path_voxels,
            total_length = total_length,
            confidence = confidence)
	

end

# ╔═╡ 847bb24d-403e-438a-8b55-9df47f118d65
hitsdf, is_double_beta, nx =get_data(datatype="bb0nu")

# ╔═╡ 42d4cb93-1d37-4f9d-8435-28e70de00485
begin
	method="KDT"
	energy_weighting=true
	edge_energy_weighting=true
	kdn = 10
	
	
	nevent=nx[9]
	md"- Examine event $(nevent)"
end

# ╔═╡ 76fea07a-2868-43f2-a7ca-6e89ddeab85a
md"""
# Monte Carlo
"""

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
	mc_path = jn.Petit.compute_mc_path(event_df, mcvox_size, 
									   is_double_beta=is_double_beta)
	jn.Petit.plot_event(event_df; mc_path=mc_path)
end

# ╔═╡ 113296b3-baeb-4507-a4cf-9c71af51169d
#mc_path

# ╔═╡ 857eff1f-bcb8-46c8-b9f0-7d3fedf01098
md"""
# MC-RECO
"""

# ╔═╡ 725637cf-507c-4c39-8bf4-6d55233b1af6
md"### Voxelize MC "

# ╔═╡ 8d7a02f4-12d2-4501-80fc-52f4ec930703
md"### Reconstruct MC tracks with standard method (KDT)"
    

# ╔═╡ a605b506-e7d8-45dc-bd93-3860b78da9fa
md"""
### Reconstruct with MST
"""

# ╔═╡ e3164fe2-51d8-41f2-b998-48e76ced6e8f
md"""
# Reconstructed Tracks
"""

# ╔═╡ f93100e5-37ae-4b43-9a8b-c5a5fa00e71f
md"### Diffuse event"

# ╔═╡ c92bce64-39b7-4282-8291-9805112a554f
begin
	event_mc = jn.Petit.transform_hits_df(event_df)
	diffused_df = jn.Petit.diffuse_xyz_image_mc(event_mc;
                                             sigma_t_mm=σt,
                                             sigma_l_mm=σl,
                                             nbins=nbins,
                                             nsigma=nsigma)
	jn.Petit.plot_hits(diffused_df, energy_column=:electrons)
end

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
                               diffusion=dfpars, method=method)
	
md"""
- number of MC tracks found with KDT = $(length(mctracks))
"""
end

# ╔═╡ 4afafdfd-b783-4fcd-9942-dcb08dd77549
let
	stats = jn.Petit.track_stats(mctracks)
	for (i, n) in enumerate(stats.n_voxels)
            println("Track $i: $n voxels, $(stats.energies_keV[i]) keV")
    end
end

# ╔═╡ f0af7b62-f5fb-4dc5-a2b2-7ef546a963c7
begin
	mctrack = mctracks[1]
	mcwalk = jn.Petit.walk_track_from_extremes(mctrack, 
									use_energy_weighting=energy_weighting,
									use_edge_energy_weighting=edge_energy_weighting)
	
	println("Track 1: $(ne(mctrack.graph)) edges, $(nv(mctrack.graph)) vertices")

end

# ╔═╡ 6cbc7f44-9869-4573-9f0d-b2209b2cb307
jn.Petit.walk_result_md(mcwalk)

# ╔═╡ 86e27850-77a2-4754-89a4-a0d280aaa61f
begin
	xmcpath = jn.Petit.get_raw_path(mctrack, mcwalk.path_indices)
	xextreme_dists = jn.Petit.compute_extreme_distances(xmcpath, mc_path)
	
	jn.Petit.extreme_distances_md(xextreme_dists)
end

# ╔═╡ 79bb7e9a-a46e-4aa6-aa5f-ef5fcd0b4075
jn.Petit.plot_track_with_paths(mctrack, xmcpath, mc_path;
                              show_distances=false)

# ╔═╡ a4832767-d1f8-4b6a-8476-0051fe5c40a5
begin
	RbMC = 7.0
	blobs_mc = jn.Petit.find_blob_energies(mctrack,
                                    xmcpath;
                                    radius=RbMC) 
	jn.Petit.blobs_md(blobs_mc)
end


# ╔═╡ 7472ee45-ba69-4e4e-9fea-d545a743ff92
jn.Petit.plot_reco_track_with_voxels_and_spheres(mctrack, xmcpath, blobs_mc, RbMC)

# ╔═╡ 58a98611-46de-4815-906b-275acd8e5564
begin
	mccoords = jn.Petit.extract_coords(mctrack)
	mstmc  = jn.Petit.find_extremes_mst_diameter(mctrack, mccoords)
	mcwalk_mst = get_walk_mst(mstmc, mctrack)
	xmcpath_mst = jn.Petit.get_raw_path(mctrack, mcwalk_mst.path_indices)
	jn.Petit.path_md(xmcpath_mst)
end

# ╔═╡ 26f73590-3ef2-4703-b04d-af7fea6898bb
begin
	xextreme_dists_mst = jn.Petit.compute_extreme_distances(xmcpath_mst, mc_path)
	jn.Petit.extreme_distances_md(xextreme_dists_mst)
end

# ╔═╡ 046202f1-0a1a-4dfd-b2ce-dea30028c246
jn.Petit.plot_track_with_paths(mctrack, xmcpath_mst, mc_path;
                              show_distances=false)

# ╔═╡ a2399ccc-36c6-410f-a733-a8a9ae99ac48
begin
	blobs_mc_mst = jn.Petit.find_blob_energies(mctrack,
                                    xmcpath_mst;
                                    radius=RbMC) 
	jn.Petit.blobs_md(blobs_mc_mst)
end

# ╔═╡ 99123acc-b0e2-4634-8c22-29cf0e4f76e3
jn.Petit.plot_reco_track_with_voxels_and_spheres(mctrack, xmcpath_mst, blobs_mc_mst, RbMC)

# ╔═╡ 12c6c5dc-9174-4d9e-b7c2-a79d7d51452e
diagonse(mctrack)

# ╔═╡ 3c3272a8-0368-413e-a4f6-7e233bbcc699
md"### Voxelize event"

# ╔═╡ a831ea92-c29b-4e21-96e0-704ee82f7795
begin
	voxels = jn.Petit.voxelize_event(diffused_df, voxel_size)
	jn.Petit.plot_event(voxels)
end

# ╔═╡ 33d6f97a-ccdf-4b65-8133-3704f5fea584
 md"### Make RECO tracks"

# ╔═╡ 83891e72-df80-4b9d-9989-7768604f2f5a
begin
tracks = jn.Petit.make_tracks(voxels;
                               max_distance_mm= max_distance,
                               energy_threshold_kev=energy_threshold_keV,
                               diffusion=dfpars,
							   method=method)
md"""
- number of tracks found = $(length(tracks))
- default method = $(method)
"""
end

# ╔═╡ e28d9f0e-5a7e-46e6-acee-02c8135a88ce
begin
	track = tracks[1]
	walk_result = jn.Petit.walk_track_from_extremes(track, 
									use_energy_weighting=energy_weighting,
									use_edge_energy_weighting=edge_energy_weighting)
	jn.Petit.walk_result_md(walk_result)
end

# ╔═╡ 41024877-b6e8-49f7-8301-5490f276ca85
begin
	path = jn.Petit.get_raw_path(track, walk_result.path_indices)
	jn.Petit.path_md(path)
end

# ╔═╡ dbdfc123-6bc6-458b-93f8-603d589646d1
begin
	extreme_dists = jn.Petit.compute_extreme_distances(path, mc_path)
	jn.Petit.extreme_distances_md(extreme_dists)
end

# ╔═╡ 6460e445-f8c0-48c3-bef9-1cefc042ecad
jn.Petit.plot_track_with_paths(track, path, mc_path;
                              show_distances=false)

# ╔═╡ 84cbaa35-78cb-410f-a76e-b6ddd1633927
md"""
### MST
"""

# ╔═╡ 07b4f340-5776-492f-8aaa-c6d3a062ec29
begin
	coords = jn.Petit.extract_coords(track)
	mst   = jn.Petit.find_extremes_mst_diameter(track, coords)
	walk_mst = get_walk_mst(mst, track)
	path_mst = jn.Petit.get_raw_path(track, walk_mst.path_indices)
	jn.Petit.path_md(path_mst)
end

# ╔═╡ 5a15d002-bfa7-43df-9b26-afb7eaaed107
begin
	extreme_dists_mst = jn.Petit.compute_extreme_distances(path_mst, mc_path)
	jn.Petit.extreme_distances_md(extreme_dists_mst)
end

# ╔═╡ 35c63b2b-a3b1-4460-926e-ae4d9a63d785
jn.Petit.plot_track_with_paths(track, path_mst, mc_path;
                              show_distances=false)

# ╔═╡ e92133ec-7197-43b5-a8dd-7cdb2f87b9b7
diagonse(track)

# ╔═╡ 9192e2e3-8380-40ba-b661-f7492ccb2034
md"## Event classification"

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
	peaks = jn.Petit.find_peaks(reco_kde.kde_f, reco_kde.kde_s; prom_scale=0.2)
	jn.Petit.plot_kde_peaks(reco_kde.kde_s, reco_kde.kde_f, peaks; title="RECO KDE")
end

# ╔═╡ e27b3ee8-0a47-4baa-aa8d-9098ac802c5d
begin
	pk = jn.Petit.kde_peaks(peaks, reco_kde.kde_f)
	jn.Petit.kde_peaks_md(pk)
end

# ╔═╡ 791ae0ee-1b30-4349-938d-3539a196b549
md"""
### Blob analysis
"""

# ╔═╡ 03e5c952-1e75-4a2d-bee5-1f262d0ae38c
Rb = 10.0

# ╔═╡ e97c03c9-aaeb-497b-bd5d-8ccd63f889ff
begin
	blobs = jn.Petit.find_blob_energies(track,
                                    path;
                                    radius=Rb) 
	jn.Petit.blobs_md(blobs)
end

# ╔═╡ f0d8b3c4-db58-4210-b05f-71a8d591ec32
jn.Petit.blobs_print(blobs, method="DEAULT")

# ╔═╡ 632673fc-2cf5-4c2a-8d4b-f2a634f311f4
 jn.Petit.plot_reco_track_with_voxels_and_spheres(track, path, blobs, Rb)

# ╔═╡ 22e993d6-ae60-40bc-a15b-cba840a1e257
md"""
### MST
"""

# ╔═╡ d0c3bec7-67d0-4ca8-bfc7-df3eb040177f
begin
	blobs_mst = jn.Petit.find_blob_energies(track,
                                    path_mst;
                                    radius=Rb) 
	jn.Petit.blobs_md(blobs_mst)
end

# ╔═╡ 7f070a6e-f1be-4c9a-a870-3de17ddbf7da
jn.Petit.blobs_print(blobs_mst, method="MST")

# ╔═╡ 16dcb4a2-c690-4d1b-bfed-f2784f847001
jn.Petit.plot_reco_track_with_voxels_and_spheres(track, path_mst, blobs_mst, Rb)

# ╔═╡ 4505b84f-8556-4ce1-82ea-cb2875ffd87e
md"## Functions"

# ╔═╡ 0db0d46c-240f-4c16-8e34-2eedd80f3557
function plot_kde(reco_kde, mc_kde)
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

# ╔═╡ Cell order:
# ╠═63889de2-cfc7-11f0-9979-e5f52b76e205
# ╠═6099909b-2c73-4b13-8429-979ac06dcef3
# ╠═f8216615-54b7-4f4a-b20c-65681052a16e
# ╠═374c62af-72c5-4ca5-98e3-bc17ae89561b
# ╠═30c8283e-ef8e-42a4-9a41-422a139c0770
# ╠═9ba4863f-d64d-41d7-8e72-1c7bf1c3569b
# ╠═49ee4ee5-f4e4-4d06-9658-e8961a261702
# ╠═aa356ade-9c1a-4e73-8739-212379c543f9
# ╠═e549464d-6e19-4d15-a64a-0f62419f80b3
# ╠═57895dc3-45d5-497b-8034-00950ce45fde
# ╠═9ed059a1-b269-40f9-a0dd-520579b6a9b8
# ╠═941446d9-c2f8-45d7-b998-3ebba4a9812e
# ╠═dc29ae3b-6725-413c-be12-228cb548f42d
# ╠═2390103c-c761-456b-842e-eeb31d9a6c0c
# ╠═0cdaf271-c87e-4cb4-aadf-428e4af466bc
# ╠═847bb24d-403e-438a-8b55-9df47f118d65
# ╠═42d4cb93-1d37-4f9d-8435-28e70de00485
# ╠═f0d8b3c4-db58-4210-b05f-71a8d591ec32
# ╠═7f070a6e-f1be-4c9a-a870-3de17ddbf7da
# ╠═76fea07a-2868-43f2-a7ca-6e89ddeab85a
# ╠═939cbc40-985c-4963-967f-3eb1f1807c1d
# ╠═b748e2a3-0504-4b2f-9d80-365081dc21f6
# ╠═64e87873-f760-4d8f-8837-427231af7951
# ╠═5746835d-4a6a-4680-b97e-5131c0649179
# ╠═113296b3-baeb-4507-a4cf-9c71af51169d
# ╠═857eff1f-bcb8-46c8-b9f0-7d3fedf01098
# ╠═725637cf-507c-4c39-8bf4-6d55233b1af6
# ╠═9e5c887b-bd64-4b12-8231-531499323288
# ╠═8d7a02f4-12d2-4501-80fc-52f4ec930703
# ╠═20922702-cf5f-4bf5-8be9-1a3c96f11e30
# ╠═4afafdfd-b783-4fcd-9942-dcb08dd77549
# ╠═f0af7b62-f5fb-4dc5-a2b2-7ef546a963c7
# ╠═6cbc7f44-9869-4573-9f0d-b2209b2cb307
# ╠═86e27850-77a2-4754-89a4-a0d280aaa61f
# ╠═79bb7e9a-a46e-4aa6-aa5f-ef5fcd0b4075
# ╠═a4832767-d1f8-4b6a-8476-0051fe5c40a5
# ╠═7472ee45-ba69-4e4e-9fea-d545a743ff92
# ╠═a605b506-e7d8-45dc-bd93-3860b78da9fa
# ╠═58a98611-46de-4815-906b-275acd8e5564
# ╠═26f73590-3ef2-4703-b04d-af7fea6898bb
# ╠═046202f1-0a1a-4dfd-b2ce-dea30028c246
# ╠═a2399ccc-36c6-410f-a733-a8a9ae99ac48
# ╠═99123acc-b0e2-4634-8c22-29cf0e4f76e3
# ╠═12c6c5dc-9174-4d9e-b7c2-a79d7d51452e
# ╠═e3164fe2-51d8-41f2-b998-48e76ced6e8f
# ╠═f93100e5-37ae-4b43-9a8b-c5a5fa00e71f
# ╠═c92bce64-39b7-4282-8291-9805112a554f
# ╠═3c3272a8-0368-413e-a4f6-7e233bbcc699
# ╠═a831ea92-c29b-4e21-96e0-704ee82f7795
# ╠═33d6f97a-ccdf-4b65-8133-3704f5fea584
# ╠═83891e72-df80-4b9d-9989-7768604f2f5a
# ╠═e28d9f0e-5a7e-46e6-acee-02c8135a88ce
# ╠═41024877-b6e8-49f7-8301-5490f276ca85
# ╠═dbdfc123-6bc6-458b-93f8-603d589646d1
# ╠═6460e445-f8c0-48c3-bef9-1cefc042ecad
# ╠═84cbaa35-78cb-410f-a76e-b6ddd1633927
# ╠═07b4f340-5776-492f-8aaa-c6d3a062ec29
# ╠═5a15d002-bfa7-43df-9b26-afb7eaaed107
# ╠═35c63b2b-a3b1-4460-926e-ae4d9a63d785
# ╠═e92133ec-7197-43b5-a8dd-7cdb2f87b9b7
# ╠═9192e2e3-8380-40ba-b661-f7492ccb2034
# ╠═02ba72f8-6437-45b4-a6cc-b2554a2317c9
# ╠═e4f9a433-fbea-4068-ae2b-66556365971f
# ╠═e27b3ee8-0a47-4baa-aa8d-9098ac802c5d
# ╠═791ae0ee-1b30-4349-938d-3539a196b549
# ╠═03e5c952-1e75-4a2d-bee5-1f262d0ae38c
# ╠═e97c03c9-aaeb-497b-bd5d-8ccd63f889ff
# ╠═632673fc-2cf5-4c2a-8d4b-f2a634f311f4
# ╠═22e993d6-ae60-40bc-a15b-cba840a1e257
# ╠═d0c3bec7-67d0-4ca8-bfc7-df3eb040177f
# ╠═16dcb4a2-c690-4d1b-bfed-f2784f847001
# ╠═4505b84f-8556-4ce1-82ea-cb2875ffd87e
# ╠═0db0d46c-240f-4c16-8e34-2eedd80f3557
