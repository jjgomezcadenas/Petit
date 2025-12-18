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

# ╔═╡ 2390103c-c761-456b-842e-eeb31d9a6c0c
jn.Petit.diffusion_params_md(dfpars)

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


# ╔═╡ 847bb24d-403e-438a-8b55-9df47f118d65
hitsdf, is_double_beta, nx =get_data(datatype="xe137")

# ╔═╡ 44ee1931-6b7d-49b5-a1b6-d78c41c739b9
md"""
### Events bb0nu
- Event 3: OK with both.
- Event 6: Better with edge\_energy\_weighting=false
- Event 9: Better with edge\_energy\_weighting=false
- Event 14: Better with edge\_energy\_weighting=false
- Event 15: Better with both
- Event 20: OK with both (better with KDT)
- Event 21: OK with both (better with KDT)
- Event 22: OK with both 
- Event 26: OK with both 
- Event 28: OK with both 
"""

# ╔═╡ fb85a17b-010d-4906-bc75-66624f01958f
md"""
### Events bb0nu
- Event 377: Border case better with edge\_energy\_weighting=true.
- Event 348: Reject with both.
- Event 341: Reject with both. better with edge\_energy\_weighting=true.
- Event 334: Better with edge\_energy\_weighting=true.
- Event 334: Better with edge\_energy\_weighting=true.
- Event 192: Fails with both
- Event 190: Better with edge\_energy\_weighting=true.
- Event 86:  Better with edge\_energy\_weighting=true.
- Event 26: OK with both 
"""

# ╔═╡ 42d4cb93-1d37-4f9d-8435-28e70de00485
begin
	method="KDT"
	energy_weighting=true
	edge_energy_weighting=true
	kdn = 10
	
	
	nevent=nx[2]
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

# ╔═╡ 36ef9ef8-dd58-44c7-b404-4fab5a5030c0
event_df.event_id[1]

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

# ╔═╡ 725637cf-507c-4c39-8bf4-6d55233b1af6
md"### Voxelize MC "

# ╔═╡ 8d7a02f4-12d2-4501-80fc-52f4ec930703
    md"### Reconstruct MC tracks"
    

# ╔═╡ f919c0c0-4c57-41d2-8a64-33683f90dae3
md"- kNN"

# ╔═╡ 4466af45-f004-4131-9e92-54f5ad6685f1
md"- KDT"

# ╔═╡ f5da0289-3728-47a2-bf6b-61cf216d8f11
md"- kNN"

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
	
	mctracks_kdt = jn.Petit.make_tracks(mcvoxels;
                               max_distance_mm=2.5 * mcvox_size,
                               energy_threshold_kev=1.0,
                               diffusion=dfpars, method=method,
                    		   k=kdn)
	
md"""
- number of MC tracks found with KDT = $(length(mctracks_kdt))
"""
end

# ╔═╡ 4afafdfd-b783-4fcd-9942-dcb08dd77549
let
	stats = jn.Petit.track_stats(mctracks_kdt)
	for (i, n) in enumerate(stats.n_voxels)
                println("Track $i: $n voxels, $(stats.energies_keV[i]) keV")
    end
end

# ╔═╡ 373e2a60-9da5-415e-bb4c-a4bea00f5b2c
begin
	mctracks_knn = jn.Petit.make_tracks(mcvoxels;
                               max_distance_mm=2.5 * mcvox_size,
                               energy_threshold_kev=1.0,
                               diffusion=dfpars, method="kNN",
                    		   k=kdn)
	
md"""
- number of MC tracks found with kNN = $(length(mctracks_knn))
"""
end

# ╔═╡ 36d6ff14-d9be-4704-8f34-1a041afd0d4b
let
	stats = jn.Petit.track_stats(mctracks_knn)
	for (i, n) in enumerate(stats.n_voxels)
                println("Track $i: $n voxels, $(stats.energies_keV[i]) keV")
    end
end

# ╔═╡ f0af7b62-f5fb-4dc5-a2b2-7ef546a963c7
begin
	mctrack_kdt = mctracks_kdt[1]
	mctrack_knn = mctracks_knn[1]
	mcwalk_kdt = jn.Petit.walk_track_from_extremes(mctrack_kdt, 
									use_energy_weighting=energy_weighting,
									use_edge_energy_weighting=edge_energy_weighting)
	mcwalk_knn = jn.Petit.walk_track_from_extremes(mctrack_knn,
									use_energy_weighting=energy_weighting,
									use_edge_energy_weighting=edge_energy_weighting)
	md"- KDT"
end

# ╔═╡ cae81030-d08b-44d7-b685-26da591e4960
println("KDT: $(ne(mctrack_kdt.graph)) edges, $(nv(mctrack_kdt.graph)) vertices")


# ╔═╡ ef2b5ebc-455c-4713-ad0a-7fd10cd282b9
println("kNN: $(ne(mctrack_knn.graph)) edges, $(nv(mctrack_knn.graph)) vertices")


# ╔═╡ 6a37fd01-798d-4043-9442-db8f6b7a77f4
jn.Petit.walk_result_md(mcwalk_kdt)

# ╔═╡ 6cbc7f44-9869-4573-9f0d-b2209b2cb307
jn.Petit.walk_result_md(mcwalk_knn)

# ╔═╡ 86e27850-77a2-4754-89a4-a0d280aaa61f
begin
	xmcpath_kdt = jn.Petit.get_raw_path(mctrack_kdt, mcwalk_kdt.path_indices)
	xextreme_dists_kdt = jn.Petit.compute_extreme_distances(xmcpath_kdt, mc_path)
	
	md"- KDT"
end

# ╔═╡ 12af48a9-86cd-47c3-9464-3598f845ad7d
jn.Petit.extreme_distances_md(xextreme_dists_kdt)

# ╔═╡ e42c689f-c966-40c4-a5a5-7d09d99773f5
begin
	xmcpath_knn = jn.Petit.get_raw_path(mctrack_knn, mcwalk_knn.path_indices)
	xextreme_dists_knn = jn.Petit.compute_extreme_distances(xmcpath_knn, mc_path)
	
	md"- kNN"
end

# ╔═╡ c5b1ddce-d041-4afc-939e-c93e2134eade
jn.Petit.extreme_distances_md(xextreme_dists_knn)

# ╔═╡ 79bb7e9a-a46e-4aa6-aa5f-ef5fcd0b4075
jn.Petit.plot_track_with_paths(mctrack_kdt, xmcpath_kdt, mc_path;
                              show_distances=false)

# ╔═╡ 94c67d3c-3ca9-47d2-998b-93b55baa28b7
jn.Petit.plot_track_with_paths(mctrack_knn, xmcpath_knn, mc_path;
                              show_distances=false)

# ╔═╡ 678f4e07-65af-4e8c-a4c1-84cb3bf06175
begin
	
  #@time diffusedK_df = jn.Petit.diffuse_xyz_image_kernel(event_mc; sigma_t_mm=σt,
   #                                          sigma_l_mm=σl,
    #                                         nbins=nbins,
     #                                        nsigma=nsigma)

end

# ╔═╡ 1a696ba7-6b03-4aa4-a70b-64cf012ef441
#jn.Petit.plot_hits(diffusedK_df, energy_column=:electrons)

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
							   method=method,
                    		   k=kdn)
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

# ╔═╡ 3ffe4b80-2694-469f-abc5-376cd6268f17
md"""
#### kNN
"""

# ╔═╡ 437695a7-23f9-4ab7-8a7a-8628e24586e0
begin
tracks_knn = jn.Petit.make_tracks(voxels;
                               max_distance_mm= max_distance,
                               energy_threshold_kev=energy_threshold_keV,
                               diffusion=dfpars,
							   method="kNN",
                    		   k=kdn)
md"""
- method = kNN
- number of tracks found = $(length(tracks))

"""
end

# ╔═╡ 92cfa67b-2bd8-403e-99c8-fc640a53fe4b
begin
	track_knn = tracks_knn[1]
	walk_result_knn = jn.Petit.walk_track_from_extremes(track_knn, 
									use_energy_weighting=energy_weighting,
									use_edge_energy_weighting=edge_energy_weighting)
	jn.Petit.walk_result_md(walk_result_knn)
end

# ╔═╡ 33b3c97b-0975-4a19-bd3f-25c496ca35d1
begin
	path_knn = jn.Petit.get_raw_path(track_knn, walk_result_knn.path_indices)
	jn.Petit.path_md(path_knn)
end

# ╔═╡ 02db2168-4dae-42e9-a87c-e4dceda037ff
begin
	extreme_dists_knn = jn.Petit.compute_extreme_distances(path_knn, mc_path)
	jn.Petit.extreme_distances_md(extreme_dists_knn)
end

# ╔═╡ 912a30e9-f265-4897-bccc-f3a6a12ba41f
jn.Petit.plot_track_with_paths(track_knn, path_knn, mc_path;
                              show_distances=false)

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

# ╔═╡ 0f6364a8-a683-4f88-a523-a26df92c34d9
#plot_kde(reco_kde, mc_kde)

# ╔═╡ d5bda625-3ab9-4eb8-83d7-9192888abe08
#md"### Find peaks " 

# ╔═╡ 9c237e6d-e798-4b61-b7d8-21f2fd8f5a6e
begin
	#peaks = jn.Petit.find_peaks(reco_kde.kde_f, reco_kde.kde_s; prom_scale=0.2)
	#jn.Petit.plot_kde_peaks(reco_kde.kde_s, reco_kde.kde_f, peaks; title="RECO KDE")
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

# ╔═╡ 1bf746ef-a41b-4875-a92a-572ac7a770b0
jn.Petit.blobs_md(blobs)

# ╔═╡ 632673fc-2cf5-4c2a-8d4b-f2a634f311f4
 jn.Petit.plot_reco_track_with_voxels_and_spheres(track, path, blobs, Rb)

# ╔═╡ 55519b5a-f26d-47ee-b7bd-f62c6064a9d0
md"""
### kNN
"""

# ╔═╡ 7c9f9889-48f3-48c7-965b-0c3a16607ed1
begin
	 reco_kde_knn = jn.Petit.get_reco_kde(track_knn, path_knn; 
										  bandwidth=kde_bandwidth, 
									      n_eval=n_kde_eval)
	peaks_knn = jn.Petit.find_peaks(reco_kde_knn.kde_f, reco_kde_knn.kde_s;
									prom_scale=0.2)
	jn.Petit.plot_kde_peaks(reco_kde_knn.kde_s, reco_kde_knn.kde_f, peaks_knn; 
							title="RECO KDE: kNN")
end

# ╔═╡ 9de306eb-a4a4-4f26-9966-cc48231dfbe8
begin
	pk_knn = jn.Petit.kde_peaks(peaks_knn, reco_kde_knn.kde_f)
	jn.Petit.kde_peaks_md(pk_knn)
end

# ╔═╡ bb59522c-4858-430e-8e41-a366f57d6386
begin
	blobs_knn = jn.Petit.find_blob_energies(track_knn,
                                    path_knn;
                                    radius=Rb) 
	jn.Petit.blobs_md(blobs_knn)
end

# ╔═╡ 1441f1b7-8da1-49e0-827e-186b76f5c62e
 jn.Petit.plot_reco_track_with_voxels_and_spheres(track_knn, path_knn, blobs_knn, Rb)

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
# ╠═9ed059a1-b269-40f9-a0dd-520579b6a9b8
# ╠═941446d9-c2f8-45d7-b998-3ebba4a9812e
# ╠═dc29ae3b-6725-413c-be12-228cb548f42d
# ╠═2390103c-c761-456b-842e-eeb31d9a6c0c
# ╠═aa356ade-9c1a-4e73-8739-212379c543f9
# ╠═847bb24d-403e-438a-8b55-9df47f118d65
# ╠═44ee1931-6b7d-49b5-a1b6-d78c41c739b9
# ╠═fb85a17b-010d-4906-bc75-66624f01958f
# ╠═42d4cb93-1d37-4f9d-8435-28e70de00485
# ╠═36ef9ef8-dd58-44c7-b404-4fab5a5030c0
# ╠═1bf746ef-a41b-4875-a92a-572ac7a770b0
# ╠═939cbc40-985c-4963-967f-3eb1f1807c1d
# ╠═b748e2a3-0504-4b2f-9d80-365081dc21f6
# ╠═64e87873-f760-4d8f-8837-427231af7951
# ╠═5746835d-4a6a-4680-b97e-5131c0649179
# ╠═113296b3-baeb-4507-a4cf-9c71af51169d
# ╠═725637cf-507c-4c39-8bf4-6d55233b1af6
# ╠═9e5c887b-bd64-4b12-8231-531499323288
# ╠═8d7a02f4-12d2-4501-80fc-52f4ec930703
# ╠═20922702-cf5f-4bf5-8be9-1a3c96f11e30
# ╠═373e2a60-9da5-415e-bb4c-a4bea00f5b2c
# ╠═4afafdfd-b783-4fcd-9942-dcb08dd77549
# ╠═36d6ff14-d9be-4704-8f34-1a041afd0d4b
# ╠═f0af7b62-f5fb-4dc5-a2b2-7ef546a963c7
# ╠═cae81030-d08b-44d7-b685-26da591e4960
# ╠═ef2b5ebc-455c-4713-ad0a-7fd10cd282b9
# ╠═6a37fd01-798d-4043-9442-db8f6b7a77f4
# ╟─f919c0c0-4c57-41d2-8a64-33683f90dae3
# ╠═6cbc7f44-9869-4573-9f0d-b2209b2cb307
# ╠═86e27850-77a2-4754-89a4-a0d280aaa61f
# ╠═12af48a9-86cd-47c3-9464-3598f845ad7d
# ╠═e42c689f-c966-40c4-a5a5-7d09d99773f5
# ╠═c5b1ddce-d041-4afc-939e-c93e2134eade
# ╠═4466af45-f004-4131-9e92-54f5ad6685f1
# ╠═79bb7e9a-a46e-4aa6-aa5f-ef5fcd0b4075
# ╠═f5da0289-3728-47a2-bf6b-61cf216d8f11
# ╠═94c67d3c-3ca9-47d2-998b-93b55baa28b7
# ╠═f93100e5-37ae-4b43-9a8b-c5a5fa00e71f
# ╠═c92bce64-39b7-4282-8291-9805112a554f
# ╠═678f4e07-65af-4e8c-a4c1-84cb3bf06175
# ╠═1a696ba7-6b03-4aa4-a70b-64cf012ef441
# ╠═3c3272a8-0368-413e-a4f6-7e233bbcc699
# ╠═a831ea92-c29b-4e21-96e0-704ee82f7795
# ╠═33d6f97a-ccdf-4b65-8133-3704f5fea584
# ╠═83891e72-df80-4b9d-9989-7768604f2f5a
# ╠═e28d9f0e-5a7e-46e6-acee-02c8135a88ce
# ╠═41024877-b6e8-49f7-8301-5490f276ca85
# ╠═dbdfc123-6bc6-458b-93f8-603d589646d1
# ╠═6460e445-f8c0-48c3-bef9-1cefc042ecad
# ╠═3ffe4b80-2694-469f-abc5-376cd6268f17
# ╠═437695a7-23f9-4ab7-8a7a-8628e24586e0
# ╠═92cfa67b-2bd8-403e-99c8-fc640a53fe4b
# ╠═33b3c97b-0975-4a19-bd3f-25c496ca35d1
# ╠═02db2168-4dae-42e9-a87c-e4dceda037ff
# ╠═912a30e9-f265-4897-bccc-f3a6a12ba41f
# ╠═9192e2e3-8380-40ba-b661-f7492ccb2034
# ╠═02ba72f8-6437-45b4-a6cc-b2554a2317c9
# ╠═e4f9a433-fbea-4068-ae2b-66556365971f
# ╠═0f6364a8-a683-4f88-a523-a26df92c34d9
# ╠═d5bda625-3ab9-4eb8-83d7-9192888abe08
# ╠═9c237e6d-e798-4b61-b7d8-21f2fd8f5a6e
# ╠═e27b3ee8-0a47-4baa-aa8d-9098ac802c5d
# ╠═791ae0ee-1b30-4349-938d-3539a196b549
# ╠═03e5c952-1e75-4a2d-bee5-1f262d0ae38c
# ╠═e97c03c9-aaeb-497b-bd5d-8ccd63f889ff
# ╠═632673fc-2cf5-4c2a-8d4b-f2a634f311f4
# ╠═55519b5a-f26d-47ee-b7bd-f62c6064a9d0
# ╠═7c9f9889-48f3-48c7-965b-0c3a16607ed1
# ╠═9de306eb-a4a4-4f26-9966-cc48231dfbe8
# ╠═bb59522c-4858-430e-8e41-a366f57d6386
# ╠═1441f1b7-8da1-49e0-827e-186b76f5c62e
# ╠═4505b84f-8556-4ce1-82ea-cb2875ffd87e
# ╠═0db0d46c-240f-4c16-8e34-2eedd80f3557
