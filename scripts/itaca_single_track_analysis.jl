#!/usr/bin/env julia

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)

using Printf
using DataFrames
using CSV
using Glob
using ArgParse
using Plots
using Dates

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit
# Helper functions (get_sigma, get_voxel_size_and_distance, etc.)
# are now part of Petit module via itaca_functions.jl

const cmdir = joinpath(ENV["DATA"], "HD5t/itaca")

#=============================================================================
# Helper Functions
=============================================================================#

"""
    save_plot(plt, outdir, nevent, name)

Save plot to file without displaying.
"""
function save_plot(plt, outdir, nevent, name)
    mkpath(outdir)
    filename = joinpath(outdir, "event_$(nevent)_$(name).png")
    savefig(plt, filename)
    println("  Saved: $filename")
end

"""
    save_and_continue(plt, outdir, nevent, name)

Save plot to file, display, and wait for Enter.
"""
function save_and_continue(plt, outdir, nevent, name)
    mkpath(outdir)
    filename = joinpath(outdir, "event_$(nevent)_$(name).png")
    savefig(plt, filename)
    display(plt)
    println("\nSaved: $filename")
    println("Press Enter to continue...")
    readline()
end

# Note: get_sigma, get_voxel_size_and_distance, get_energy_threshold
# are now available via Petit module

"""
    get_walk_mst(res_mst, track)

Convert MST extreme finding result to walk_result format.
"""
function get_walk_mst(res_mst, track)
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
        total_length += Petit.euclidean_distance(
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

"""
    print_event_results(nevent, track_length, blobs, extreme_dists)

Pretty print the results for an event (CSV row contents).
Note: blobs.Eb1 and blobs.Eb2 are already in keV.
"""
function print_event_results(nevent, track_length, blobs, extreme_dists)
    println("┌─────────────────────────────────────────────────────────────────┐")
    println("│  Event $nevent Results                                          ")
    println("├─────────────────────────────────────────────────────────────────┤")
    println("│  Track length:    $(lpad(@sprintf("%.2f", track_length), 8)) mm                          │")
    println("├─────────────────────────────────────────────────────────────────┤")
    println("│  Blob Analysis:                                                 │")
    println("│    Eb1 (high):    $(lpad(@sprintf("%.1f", blobs.Eb1), 8)) keV                        │")
    println("│    Eb2 (low):     $(lpad(@sprintf("%.1f", blobs.Eb2), 8)) keV                        │")
    println("│    Asymmetry:     $(lpad(@sprintf("%.4f", blobs.asymmetry), 8))                            │")
    println("├─────────────────────────────────────────────────────────────────┤")
    println("│  Extreme Distances (RECO vs MC):                                │")
    println("│    d1:            $(lpad(@sprintf("%.2f", extreme_dists.d1), 8)) mm                          │")
    println("│    d2:            $(lpad(@sprintf("%.2f", extreme_dists.d2), 8)) mm                          │")
    println("└─────────────────────────────────────────────────────────────────┘")
end

"""
    plot_kde_comparison(reco_kde, mc_kde, nevent, kde_bandwidth, n_kde_eval)

Plot KDE and histogram comparison between RECO and MC.
"""
function plot_kde_comparison(reco_kde, mc_kde, nevent, kde_bandwidth, n_kde_eval)
    # KDE comparison
    p1 = plot(reco_kde.kde_s, reco_kde.kde_f,
              xlabel="Arc length s (mm)",
              ylabel="KDE Energy density f(s)",
              title="Longitudinal Energy Density (Event $nevent)",
              label="RECO KDE (h=$(round(kde_bandwidth, digits=1))mm)",
              linewidth=2,
              color=:blue)
    plot!(p1, mc_kde.kde_s, mc_kde.kde_f,
          label="MC KDE",
          linewidth=2,
          color=:red)

    # Histogram comparison
    p2 = histogram(reco_kde.s, weights=reco_kde.E .* 1e3,
                   bins=range(0, reco_kde.track_length, length=n_kde_eval+1),
                   xlabel="Arc length s (mm)",
                   ylabel="Energy (keV)",
                   title="Energy vs S (Event $nevent)",
                   label="RECO",
                   fillalpha=0.5,
                   color=:blue)
    histogram!(p2, mc_kde.s, weights=mc_kde.E .* 1e3,
               bins=range(0, mc_kde.track_length, length=n_kde_eval+1),
               label="MC",
               fillalpha=0.5,
               color=:red)

    return plot(p1, p2, layout=(1, 2), size=(1200, 500))
end

"""
    save_metadata(outdir, params, n_processed, n_single_track)

Save metadata CSV file with run parameters.
"""
function save_metadata(outdir::String, params::Dict, n_processed::Int, n_single_track::Int)
    mkpath(outdir)

    metadata_df = DataFrame(
        parameter = String[],
        value = String[]
    )

    # Add all parameters
    for (k, v) in params
        push!(metadata_df, (parameter=k, value=string(v)))
    end

    # Add processing stats
    push!(metadata_df, (parameter="n_events_processed", value=string(n_processed)))
    push!(metadata_df, (parameter="n_events_single_track", value=string(n_single_track)))
    push!(metadata_df, (parameter="timestamp", value=string(now())))

    CSV.write(joinpath(outdir, "metadata.csv"), metadata_df)
    println("\nMetadata saved to: $(joinpath(outdir, "metadata.csv"))")
end

#=============================================================================
# Main Script
=============================================================================#

"""
    main(; event_list, input_file, particle_type, lmin, lmax, ...)

Run ITACA single track analysis on specified event_ids.

# Drift Length Parameters
- `lmin`: Minimum drift length in cm (default: 100.0)
- `lmax`: Maximum drift length in cm (default: 100.0)

Behavior:
- If `lmin == lmax`: Fixed drift length mode. Diffusion parameters are computed
  once and used for all events.
- If `lmin < lmax`: Random drift length mode. For each event, a drift length is
  sampled uniformly from [lmin, lmax] and diffusion parameters are computed
  per-event.
- If `lmin > lmax`: Error - invalid configuration.

# Output
- CSV file with analysis results including per-event diffusion parameters:
  - ldrft_cm, sigma_t_mm, sigma_l_mm, voxel_size_mm, max_distance_mm, kde_bandwidth_mm
  - track_length_mm, Eb1_keV, Eb2_keV, asymmetry, d1_mm, d2_mm
  - n_peaks, peak1_left, peak1_right, peak1_prom, peak2_left, peak2_right, peak2_prom
- metadata.csv with run configuration (lmin, lmax, particle_type, etc.)
- Optional plots saved to plots/ subdirectory
"""
function main(; event_list::Vector{Int}=[1],
               input_file::String="bb0nu/bb0nu_15bar_p1.h5",
               particle_type::String="ion",
               lmin::Float64=100.0,
               lmax::Float64=100.0,
               dt::Float64=3.5, dl::Float64=0.9,
               tK::Float64=297.0, edrift::Float64=500.0, Pbar::Float64=15.0,
               energy_threshold_ions::Float64=10.0,
               energy_threshold_keV::Float64=10.0,
               nbins::Int=100, nsigma::Float64=3.0,
               n_kde_eval::Int=200, Rb::Float64=10.0,
               outdir::String="results",
               interactive::Bool=false,
               nplot::Int=0,
               print_level::String="verbose",
               graph_method::String="KDT",
               knn_k::Int=10,
               use_energy_weighting::Bool=true,
               use_edge_energy_weighting::Bool=true,
               use_mst_fallback::Bool=false)

    println("╔═══════════════════════════════════════════════════════════════╗")
    println("║         ITACA Single Track Analysis                          ║")
    println("╚═══════════════════════════════════════════════════════════════╝")

    # Validate lmin/lmax
    if lmin > lmax
        error("lmin ($lmin) must be less than or equal to lmax ($lmax)")
    end
    use_random_ldrft = (lmin < lmax)

    # Print level flags
    is_quiet = (print_level == "quiet")
    is_verbose = (print_level == "verbose")

    # If fixed drift (lmin == lmax), compute diffusion parameters once
    # Otherwise, they will be computed per event
    if !use_random_ldrft
        ldrft_fixed = lmin  # lmin == lmax
        σt_fixed, σl_fixed = Petit.get_sigma(particle_type, ldrft_fixed;
                                        dt=dt, dl=dl, tK=tK, edrift=edrift, Pbar=Pbar)
        voxel_size_fixed, mcvox_size_fixed, max_distance_fixed = Petit.get_voxel_size_and_distance(ldrft_fixed, σt_fixed)
        eth_fixed = Petit.get_energy_threshold(particle_type;
                                          energy_threshold_ions=energy_threshold_ions,
                                          energy_threshold_keV=energy_threshold_keV)
        kde_bandwidth_fixed = 2.0 * voxel_size_fixed
        mc_kde_bandwidth_fixed = 2.0 * mcvox_size_fixed
        dfpars_fixed = Petit.DiffusionParams(ldrft_fixed, σt_fixed, σl_fixed,
                                              voxel_size_fixed, max_distance_fixed,
                                              eth_fixed, nbins, nsigma)

        if is_verbose
            println("\n### Diffusion Parameters (fixed drift)")
            Petit.diffusion_params_print(dfpars_fixed)
        end
    end

    if is_verbose
        if use_random_ldrft
            println("\nDrift mode: RANDOM [$lmin, $lmax] cm")
        else
            println("Drift mode: FIXED $lmin cm")
        end
        println("Graph method: $graph_method$(graph_method in ["kNN", "kNN_mutual"] ? " (k=$knn_k)" : "")")
        println("Interactive mode: $interactive")
        println("Plot every N events: $(nplot == 0 ? "disabled" : nplot)")
        println("Print level: $print_level")
    end

    # Store parameters for metadata (only fixed params, per-event params go to CSV)
    params = Dict(
        "input_file" => input_file,
        "particle_type" => particle_type,
        "lmin_cm" => lmin,
        "lmax_cm" => lmax,
        "random_drift" => use_random_ldrft,
        "blob_radius_mm" => Rb,
        "nbins" => nbins,
        "nsigma" => nsigma,
        "n_events" => length(event_list),
        "event_list" => join(event_list, ","),
        "graph_method" => graph_method,
        "knn_k" => knn_k,
        "use_energy_weighting" => use_energy_weighting,
        "use_edge_energy_weighting" => use_edge_energy_weighting,
        "use_mst_fallback" => use_mst_fallback
    )

    # Create output directories
    mkpath(outdir)
    plots_dir = joinpath(outdir, "plots")
    mkpath(plots_dir)

    # Load data
    if is_verbose
        println("\n### Loading data from: $input_file")
    end
    hitsdf = Petit.load_data(input_file, cmdir)

    # Determine if double-beta (bb0nu) or single-electron (xe137) based on filename
    is_double_beta = occursin("bb", lowercase(input_file))
    if is_verbose
        println("  Event type: $(is_double_beta ? "double-beta (bb0nu)" : "single-electron (xe137)")")
    end

    # Initialize results DataFrame (includes per-event diffusion params)
    results_df = DataFrame(
        event = Int[],
        ldrft_cm = Float64[],
        sigma_t_mm = Float64[],
        sigma_l_mm = Float64[],
        voxel_size_mm = Float64[],
        max_distance_mm = Float64[],
        track_length_mm = Float64[],
        track_length_mst_mm = Float64[],
        # MC blob analysis (default method)
        Eb1_mc_keV = Float64[],
        Eb2_mc_keV = Float64[],
        asymmetry_mc = Float64[],
        # MC blob analysis (MST method)
        Eb1_mc_mst_keV = Float64[],
        Eb2_mc_mst_keV = Float64[],
        asymmetry_mc_mst = Float64[],
        # RECO blob analysis (default method)
        Eb1_keV = Float64[],
        Eb2_keV = Float64[],
        asymmetry = Float64[],
        # RECO blob analysis (MST method)
        Eb1_mst_keV = Float64[],
        Eb2_mst_keV = Float64[],
        asymmetry_mst = Float64[],
        # MC extreme distances (default method)
        d1_mc_mm = Float64[],
        d2_mc_mm = Float64[],
        # MC extreme distances (MST method)
        d1_mc_mst_mm = Float64[],
        d2_mc_mst_mm = Float64[],
        # RECO extreme distances (default method)
        d1_mm = Float64[],
        d2_mm = Float64[],
        # RECO extreme distances (MST method)
        d1_mst_mm = Float64[],
        d2_mst_mm = Float64[]
    )

    n_processed = 0
    n_single_track = 0

    # Process events
    for nevent in event_list
        n_processed += 1

        # Determine if we should save plots for this event
        save_plots_this_event = (nplot > 0) && (mod(nevent, nplot) == 0)

        # Determine drift length and diffusion parameters for this event
        if use_random_ldrft
            # Sample random drift length and compute diffusion parameters
            event_ldrft = lmin + rand() * (lmax - lmin)
            event_σt, event_σl = Petit.get_sigma(particle_type, event_ldrft;
                                           dt=dt, dl=dl, tK=tK, edrift=edrift, Pbar=Pbar)
            event_voxel_size, event_mcvox_size, event_max_distance = Petit.get_voxel_size_and_distance(event_ldrft, event_σt)
            event_eth = Petit.get_energy_threshold(particle_type;
                                             energy_threshold_ions=energy_threshold_ions,
                                             energy_threshold_keV=energy_threshold_keV)
            event_kde_bandwidth = 2.0 * event_voxel_size
            event_mc_kde_bandwidth = 2.0 * event_mcvox_size
            event_dfpars = Petit.DiffusionParams(event_ldrft, event_σt, event_σl,
                                                  event_voxel_size, event_max_distance,
                                                  event_eth, nbins, nsigma)
        else
            # Use precomputed fixed parameters
            event_ldrft = ldrft_fixed
            event_σt = σt_fixed
            event_σl = σl_fixed
            event_voxel_size = voxel_size_fixed
            event_mcvox_size = mcvox_size_fixed
            event_max_distance = max_distance_fixed
            event_eth = eth_fixed
            event_kde_bandwidth = kde_bandwidth_fixed
            event_mc_kde_bandwidth = mc_kde_bandwidth_fixed
            event_dfpars = dfpars_fixed
        end

        # Print event header based on print_level
        if is_verbose
            print("Event $nevent: l=$(round(event_ldrft, digits=1))cm, σt=$(round(event_σt, digits=2))mm, σl=$(round(event_σl, digits=2))mm, vox=$(round(event_voxel_size, digits=2))mm ... ")
        elseif save_plots_this_event
            # quiet mode: only print on plot events
            print("Event $nevent... ")
        end

        # Get event
        event_df = Petit.get_event(hitsdf, nevent)

        if interactive || save_plots_this_event
            plt = Petit.plot_event(event_df)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "01_mc_hits")
            else
                save_plot(plt, plots_dir, nevent, "01_mc_hits")
            end
        end

        # Compute MC path
        if is_verbose
            println("\n### Computing MC path")
        end
        mc_path = Petit.compute_mc_path(event_df, event_mcvox_size; is_double_beta=is_double_beta)

        if interactive || save_plots_this_event
            plt = Petit.plot_event(event_df; mc_path=mc_path)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "02_mc_path")
            else
                save_plot(plt, plots_dir, nevent, "02_mc_path")
            end
        end

        # Voxelize MC
        mcvoxels = Petit.voxelize_event(event_df, event_mcvox_size)

        if interactive || save_plots_this_event
            plt = Petit.plot_event(mcvoxels)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "03_mc_voxels")
            else
                save_plot(plt, plots_dir, nevent, "03_mc_voxels")
            end
        end

        # Reconstruct MC tracks (KDT by default)
        mctracks = Petit.make_tracks(mcvoxels;
                                     max_distance_mm=2.5 * event_mcvox_size,
                                     energy_threshold_kev=1.0,
                                     diffusion=event_dfpars)
        stats = Petit.track_stats(mctracks)
        if is_verbose
            println("- number of mc tracks found = $(length(mctracks))")
            for (i, n) in enumerate(stats.n_voxels)
                println("Track $i: $n voxels, $(stats.energies_keV[i]) keV")
            end
        end

        # Walk track for MC

        mctrack = mctracks[1]
	    mcwalk  = Petit.walk_track_from_extremes(mctrack, 
									use_energy_weighting=use_energy_weighting,
									use_edge_energy_weighting=use_edge_energy_weighting)

        if is_verbose
            Petit.walk_result_print(mcwalk; track=mctrack, method="MC-KDT")
        end

        ## Reconstructed MC path
        xmcpath = Petit.get_raw_path(mctrack, mcwalk.path_indices)
	    xextreme_dists = Petit.compute_extreme_distances(xmcpath, mc_path)

        if is_verbose
            Petit.extreme_distances_print(xextreme_dists; method="MC-KDT")
        end

        if interactive || save_plots_this_event
            plt = Petit.plot_track_with_paths(mctrack, xmcpath, mc_path;
                              show_distances=false)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "04_mc_reco_path")
            else
                save_plot(plt, plots_dir, nevent, "04_mc_reco_path")
            end
        end

        # Blob analysis MC Default

        if is_verbose
            println("### Blob Analysis MC DEFAULT")
        end

        blobs_mc = Petit.find_blob_energies(mctrack, xmcpath; radius=Rb)
        if is_verbose
            Petit.blobs_print(blobs_mc; method=graph_method)
        end

        if interactive || save_plots_this_event
            plt = Petit.plot_reco_track_with_voxels_and_spheres(mctrack, 
                                                                xmcpath, 
                                                                blobs_mc, 
                                                                Rb)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "05_blobs_MC")
            else
                save_plot(plt, plots_dir, nevent, "05_blobs_MC")
            end
        end

        ####
        # Walk track for MC/MST

        mccoords = Petit.extract_coords(mctrack)
	    mstmc    = Petit.find_extremes_mst_diameter(mctrack, mccoords)
	    mcwalk_mst = get_walk_mst(mstmc, mctrack)

        if is_verbose
            Petit.walk_result_print(mcwalk_mst; track=mctrack, method="MC-MST")
        end

	    xmcpath_mst = Petit.get_raw_path(mctrack, mcwalk_mst.path_indices)
        xextreme_dists_mst = Petit.compute_extreme_distances(xmcpath_mst, mc_path)

        if is_verbose
            Petit.extreme_distances_print(xextreme_dists_mst; method="MC-MST")
        end

        if interactive || save_plots_this_event
            plt = Petit.plot_track_with_paths(mctrack, xmcpath_mst, mc_path;
                              show_distances=false)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "06_mc_path_mst")
            else
                save_plot(plt, plots_dir, nevent, "06_mc_path_mst")
            end
        end

        # Blob analysis MC Default

        if is_verbose
            println("### Blob Analysis MC/MST")
        end

        blobs_mc_mst = Petit.find_blob_energies(mctrack, xmcpath_mst; radius=Rb)
        if is_verbose
            Petit.blobs_print(blobs_mc_mst; method="MC/MST")
        end

        if interactive || save_plots_this_event
            plt = Petit.plot_reco_track_with_voxels_and_spheres(mctrack, 
                                                                xmcpath_mst, 
                                                                blobs_mc_mst, 
                                                                Rb)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "07_blobs_MC_MST")
            else
                save_plot(plt, plots_dir, nevent, "07_blobs_MC_MST")
            end
        end
	
        # RECO

        # Transform hits
        if is_verbose
            println("### Transforming hits")
        end
        event_mc = Petit.transform_hits_df(event_df)

        # Diffuse event
        if is_verbose
            println("### Diffusing event")
        end
        diffused_df = Petit.diffuse_xyz_image_mc(event_mc;
                                                 sigma_t_mm=event_σt,
                                                 sigma_l_mm=event_σl,
                                                 nbins=nbins,
                                                 nsigma=nsigma)

        if interactive || save_plots_this_event
            plt = Petit.plot_hits(diffused_df, energy_column=:electrons)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "08_diffused")
            else
                save_plot(plt, plots_dir, nevent, "08_diffused")
            end
        end


        # Voxelize event
        if is_verbose
            println("### Voxelizing event")
        end
        voxels = Petit.voxelize_event(diffused_df, event_voxel_size)

        if interactive || save_plots_this_event
            plt = Petit.plot_event(voxels)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "09_voxels")
            else
                save_plot(plt, plots_dir, nevent, "09_voxels")
            end
        end

        # Make tracks
        if is_verbose
            println("### Making tracks")
        end
        tracks = Petit.make_tracks(voxels;
                                   max_distance_mm=event_max_distance,
                                   energy_threshold_kev=event_eth,
                                   diffusion=event_dfpars,
                                   method=graph_method,
                                   k=knn_k)

        n_tracks = length(tracks)
        if is_verbose
            println("Number of tracks found in RECO: $n_tracks")
        end

        if n_tracks != 1
            println("$n_tracks tracks - rejected")
            continue
        end

        n_single_track += 1
        track = tracks[1]

        # Walk track
        if is_verbose
            println("### Walking track from extremes")
        end
        walk_result = Petit.walk_track_from_extremes(track;
            use_energy_weighting=use_energy_weighting,
            use_edge_energy_weighting=use_edge_energy_weighting,
            use_mst_fallback=use_mst_fallback)

        if is_verbose
            Petit.walk_result_print(walk_result; track=track, method=graph_method)
        end

        # Get path
        path = Petit.get_raw_path(track, walk_result.path_indices)
        track_length = path.s[end]
        if is_verbose
           Petit.path_print(path; method=graph_method)
        end

        # Compute extreme distances
        if is_verbose
            println("### Computing extreme distances (RECO vs MC)")
        end
        extreme_dists = Petit.compute_extreme_distances(path, mc_path)

        if is_verbose
            Petit.extreme_distances_print(extreme_dists; method=graph_method)
        end

        if interactive || save_plots_this_event
            plt = Petit.plot_track_with_paths(track, path, mc_path; show_distances=true)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "10_track_paths")
            else
                save_plot(plt, plots_dir, nevent, "10_track_paths")
            end
        end

        ##### MST
        coords = Petit.extract_coords(track)
	    mst    = Petit.find_extremes_mst_diameter(track, coords)
	    walk_mst = get_walk_mst(mst, track)

        if is_verbose
            Petit.walk_result_print(walk_mst; track=track, method="MST")
        end

        # Get path
        path_mst = Petit.get_raw_path(track, walk_mst.path_indices)
        track_length_mst = path_mst.s[end]

        if is_verbose
           Petit.path_print(path_mst; method="MST")
        end

        # Compute extreme distances
        if is_verbose
            println("### Computing extreme distances (MST vs MC)")
        end
        extreme_dists_mst = Petit.compute_extreme_distances(path_mst, mc_path)

        if is_verbose
            Petit.extreme_distances_print(extreme_dists_mst; method="MST")
        end

        if interactive || save_plots_this_event
            plt = Petit.plot_track_with_paths(track, path_mst, mc_path;
                              show_distances=false)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "11_track_path_mst")
            else
                save_plot(plt, plots_dir, nevent, "11_track_path_mst")
            end
        end
        ####
        # KDE analysis (commented for now)
        # if is_very_verbose
        #     println("### KDE Analysis")
        # end
        # reco_kde = Petit.get_reco_kde(track, path; bandwidth=event_kde_bandwidth,
        #                               n_eval=n_kde_eval)
        # mc_kde = Petit.get_mc_kde(mc_path; bandwidth=event_mc_kde_bandwidth,
        #                           n_eval=n_kde_eval)

        # if interactive || save_plots_this_event
        #     plt = plot_kde_comparison(reco_kde, mc_kde, nevent, event_kde_bandwidth, n_kde_eval)
        #     if interactive
        #         save_and_continue(plt, plots_dir, nevent, "06_kde_comparison")
        #     else
        #         save_plot(plt, plots_dir, nevent, "06_kde_comparison")
        #     end
        # end

        # # Find peaks
        # if is_very_verbose
        #     println("### Finding peaks")
        # end
        # peaks = Petit.find_peaks(reco_kde.kde_f, reco_kde.kde_s; prom_scale=0.2)
        # if interactive
        #     print_peaks(peaks)
        # end

        # # Extract peak1 (leftmost) and peak2 (rightmost) info
        # pk = kde_peaks(peaks, reco_kde.kde_f)
        # if is_very_verbose
        #     println("    Peak1: left=$(round(pk.peak1_left, digits=1)), right=$(round(pk.peak1_right, digits=1)), prom=$(round(pk.peak1_prom, digits=3))")
        #     println("    Peak2: left=$(round(pk.peak2_left, digits=1)), right=$(round(pk.peak2_right, digits=1)), prom=$(round(pk.peak2_prom, digits=3))")
        # end

        # if interactive || save_plots_this_event
        #     plt = Petit.plot_kde_peaks(reco_kde.kde_s, reco_kde.kde_f, peaks;
        #                                title="RECO KDE (Event $nevent)")
        #     if interactive
        #         save_and_continue(plt, plots_dir, nevent, "07_kde_peaks")
        #     else
        #         save_plot(plt, plots_dir, nevent, "07_kde_peaks")
        #     end
        # end

        # Blob analysis

        if is_verbose
            println("### Blob Analysis DEFAULT")
        end

        blobs = Petit.find_blob_energies(track, path; radius=Rb)
        if is_verbose
            Petit.blobs_print(blobs; method=graph_method)
        end

        if interactive || save_plots_this_event
            plt = Petit.plot_reco_track_with_voxels_and_spheres(track, path, blobs, Rb)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "12_blobs_reco")
            else
                save_plot(plt, plots_dir, nevent, "12_blobs_reco")
            end
        end


        ####
        # Blob analysis MST

        if is_verbose
            println("### Blob Analysis MST")
        end

        blobs_mst = Petit.find_blob_energies(track, path_mst; radius=Rb)
        if is_verbose
            Petit.blobs_print(blobs_mst; method="MST")
        end

        if interactive || save_plots_this_event
            plt = Petit.plot_reco_track_with_voxels_and_spheres(track, 
                                                                path_mst, 
                                                                blobs_mst, 
                                                                Rb)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "13_blobs_reco_mst")
            else
                save_plot(plt, plots_dir, nevent, "13_blobs_reco_mst")
            end
        end
        ###

        # Add row to results DataFrame (including diffusion params)
        push!(results_df, (
            event = nevent,
            ldrft_cm = event_ldrft,
            sigma_t_mm = event_σt,
            sigma_l_mm = event_σl,
            voxel_size_mm = event_voxel_size,
            max_distance_mm = event_max_distance,
            track_length_mm = track_length,
            track_length_mst_mm = track_length_mst,
            Eb1_mc_keV = blobs_mc.Eb1,
            Eb2_mc_keV = blobs_mc.Eb2,
            asymmetry_mc = blobs_mc.asymmetry,
            Eb1_mc_mst_keV = blobs_mc_mst.Eb1,
            Eb2_mc_mst_keV = blobs_mc_mst.Eb2,
            asymmetry_mc_mst = blobs_mc_mst.asymmetry,
            Eb1_keV = blobs.Eb1,
            Eb2_keV = blobs.Eb2,
            asymmetry = blobs.asymmetry,
            Eb1_mst_keV = blobs_mst.Eb1,
            Eb2_mst_keV = blobs_mst.Eb2,
            asymmetry_mst = blobs_mst.asymmetry,
            d1_mc_mm = xextreme_dists.d1,
            d2_mc_mm = xextreme_dists.d2,
            d1_mc_mst_mm = xextreme_dists_mst.d1,
            d2_mc_mst_mm = xextreme_dists_mst.d2,
            d1_mm = extreme_dists.d1,
            d2_mm = extreme_dists.d2,
            d1_mst_mm = extreme_dists_mst.d1,
            d2_mst_mm = extreme_dists_mst.d2
        ))

        # Print completion message based on print_level
        if is_verbose
            print_event_results(nevent, track_length, blobs, extreme_dists)
            println("✓ Event $nevent completed")
        elseif save_plots_this_event
            println("done")
        end
    end

    # Save results CSV
    results_file = joinpath(outdir, "itaca_single_track_analysis_results.csv")
    CSV.write(results_file, results_df)
    println("\n✓ Results saved to: $results_file")

    # Save metadata
    save_metadata(outdir, params, n_processed, n_single_track)

    println("\n╔═══════════════════════════════════════════════════════════════╗")
    println("║  Analysis Complete                                            ║")
    println("╠═══════════════════════════════════════════════════════════════╣")
    println("║  Events processed:    $(lpad(n_processed, 5))                              ║")
    println("║  Single-track events: $(lpad(n_single_track, 5))                              ║")
    println("║  Output directory:    $(outdir)")
    println("╚═══════════════════════════════════════════════════════════════╝")
end


#=============================================================================
# Command Line Interface
=============================================================================#

if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParseSettings(description="ITACA Single Track Analysis")

    @add_arg_table! s begin
        "--event", "-i"
            help = "Single event_id to process"
            arg_type = Int
            default = -1
        "--list", "-l"
            help = "Comma-separated list of event_ids (e.g., '4,10,17,25')"
            arg_type = String
            default = ""
        "--ni"
            help = "Start event index (1-based) from H5 file"
            arg_type = Int
            default = 0
        "--nl"
            help = "End event index from H5 file (0 = to end)"
            arg_type = Int
            default = 0
        "--input", "-f"
            help = "Input file (relative to DATA/HD5t/itaca)"
            arg_type = String
            default = "bb0nu/bb0nu_15bar_p1.h5"
        "--particle", "-p"
            help = "Particle type: 'ion' or 'electron'"
            arg_type = String
            default = "ion"
        "--lmin"
            help = "Minimum drift length (cm). If lmin == lmax, fixed drift. If lmin < lmax, random drift in [lmin, lmax]"
            arg_type = Float64
            default = 100.0
        "--lmax"
            help = "Maximum drift length (cm). If lmin == lmax, fixed drift. If lmin < lmax, random drift in [lmin, lmax]"
            arg_type = Float64
            default = 100.0
        "--dt"
            help = "Transverse diffusion coefficient (mm/√cm)"
            arg_type = Float64
            default = 3.5
        "--dl"
            help = "Longitudinal diffusion coefficient (mm/√cm)"
            arg_type = Float64
            default = 0.9
        "--tK"
            help = "Temperature in Kelvin"
            arg_type = Float64
            default = 297.0
        "--edrift"
            help = "Drift field in V/cm"
            arg_type = Float64
            default = 500.0
        "--Pbar"
            help = "Pressure in bar"
            arg_type = Float64
            default = 15.0
        "--eth-ion"
            help = "Energy threshold for ions (keV)"
            arg_type = Float64
            default = 10.0
        "--eth"
            help = "Energy threshold (keV)"
            arg_type = Float64
            default = 10.0
        "--nbins"
            help = "Number of bins for diffusion"
            arg_type = Int
            default = 100
        "--nsigma"
            help = "Number of sigma for diffusion"
            arg_type = Float64
            default = 3.0
        "--nkde"
            help = "Number of KDE evaluation points"
            arg_type = Int
            default = 200
        "--Rb"
            help = "Blob radius in mm"
            arg_type = Float64
            default = 10.0
        "--outdir", "-o"
            help = "Output directory for results"
            arg_type = String
            default = "results"
        "--interactive"
            help = "Enable interactive mode (display plots and wait for Enter)"
            action = :store_true
        "--nplot"
            help = "Save plots every N events (0 = disabled)"
            arg_type = Int
            default = 0
        "--print_level"
            help = "Print level: 'quiet' (minimal) or 'verbose' (detailed output)"
            arg_type = String
            default = "verbose"
        "--graph-method", "-g"
            help = "Graph construction method: 'KDT' (radius), 'kNN' (k nearest neighbors), 'kNN_mutual' (mutual kNN)"
            arg_type = String
            default = "KDT"
        "--knn-k"
            help = "Number of neighbors for kNN methods (default: 10)"
            arg_type = Int
            default = 10
        "--no-energy-weighting"
            help = "Disable energy-weighted extreme finding"
            action = :store_true
        "--no-edge-energy-weighting"
            help = "Disable edge-energy-weighted Dijkstra method"
            action = :store_true
        "--mst-fallback"
            help = "Enable MST-based fallback method for extreme finding"
            action = :store_true
    end

    args = parse_args(s)

    # Build event list from arguments
    event_list = Int[]
    if args["ni"] > 0 || args["nl"] > 0
        # Load data to get available events
        input_file = args["input"]
        println("Loading H5 file to get event list...")
        hitsdf = Petit.load_data(input_file, cmdir)
        all_events = sort(unique(hitsdf.event_id))
        n_total = length(all_events)

        ni = args["ni"] > 0 ? args["ni"] : 1
        nl = args["nl"] > 0 ? args["nl"] : n_total

        if ni > n_total
            error("Start index ni=$ni exceeds available events ($n_total)")
        end
        nl = min(nl, n_total)
        if nl < ni
            error("End index nl=$nl must be >= start index ni=$ni")
        end
        event_list = all_events[ni:nl]
        println("Processing events at indices $ni:$nl ($(length(event_list)) events out of $n_total)")
    elseif !isempty(args["list"])
        # Parse comma-separated list
        event_list = [parse(Int, strip(s)) for s in split(args["list"], ",")]
    elseif args["event"] >= 0
        # Single event
        event_list = [args["event"]]
    else
        error("Must specify --event (-i), --list (-l), or --ni/--nl range")
    end

    main(; event_list=event_list,
           input_file=args["input"],
           particle_type=args["particle"],
           lmin=args["lmin"],
           lmax=args["lmax"],
           dt=args["dt"],
           dl=args["dl"],
           tK=args["tK"],
           edrift=args["edrift"],
           Pbar=args["Pbar"],
           energy_threshold_ions=args["eth-ion"],
           energy_threshold_keV=args["eth"],
           nbins=args["nbins"],
           nsigma=args["nsigma"],
           n_kde_eval=args["nkde"],
           Rb=args["Rb"],
           outdir=args["outdir"],
           interactive=args["interactive"],
           nplot=args["nplot"],
           print_level=args["print_level"],
           graph_method=args["graph-method"],
           knn_k=args["knn-k"],
           use_energy_weighting=!args["no-energy-weighting"],
           use_edge_energy_weighting=!args["no-edge-energy-weighting"],
           use_mst_fallback=args["mst-fallback"])
end
