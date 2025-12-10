#!/usr/bin/env julia

"""
Multi-threaded ITACA track reconstruction and analysis script.

Combines:
- Analysis pipeline from itaca_single_track_analysis.jl (diffusion, voxelization,
  track finding, blob analysis, KDE peaks)
- Multi-threading and HDF5 persistence from track_reco_mt.jl

Outputs:
1. CSV summary file with analysis metrics for all events
2. HDF5 files with full track data for events passing cuts
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using HDF5
using DataFrames
using CSV
using Statistics
using Graphs
using Dates

# Load Petit module and helpers
include(joinpath(pdir, "src", "Petit.jl"))
include(joinpath(pdir, "src", "itaca_aux.jl"))
import .Petit

#=============================================================================
# Helper Functions (from itaca_single_track_analysis.jl)
=============================================================================#

"""
    get_sigma(particle_type, ldrft; dt, dl, tK, edrift, Pbar)

Compute transverse and longitudinal sigma based on particle type.
"""
function get_sigma(particle_type::String, ldrft::Float64;
                   dt::Float64=3.5, dl::Float64=0.9,
                   tK::Float64=297.0, edrift::Float64=500.0, Pbar::Float64=15.0)
    if particle_type == "ion"
        σt = Petit.sigma_t_ion_mm(tK, ldrft, edrift)
        σl = 0.0
    else
        σt = Petit.sigma_t_mm(ldrft, Pbar; dtmm=dt)
        σl = Petit.sigma_l_mm(ldrft, Pbar; dlmm=dl)
    end
    return σt, σl
end

"""
    get_voxel_size_and_distance(ldrft, σt)

Compute voxel size, MC voxel size, and max distance based on diffusion.
"""
function get_voxel_size_and_distance(ldrft::Float64, σt::Float64)
    if ldrft > 50.0
        voxel_scale = 2.0
        voxel_distance_scale = 1.5
    else
        voxel_scale = 3.0
        voxel_distance_scale = 2.0
    end

    voxel_size = σt * voxel_scale
    mcvox_size = voxel_size
    max_distance = voxel_size * voxel_distance_scale
    return (voxel_size, mcvox_size, max_distance)
end

"""
    get_energy_threshold(particle_type; energy_threshold_ions, energy_threshold_keV)

Get energy threshold based on particle type.
"""
function get_energy_threshold(particle_type::String;
                              energy_threshold_ions::Float64=10.0,
                              energy_threshold_keV::Float64=10.0)
    f = 1e+5 / 2.5  # ions per MeV
    fkeV = f * 1e-3  # ions per keV

    if particle_type == "ion"
        energy_threshold_keV = energy_threshold_ions / fkeV
    end
    return energy_threshold_keV
end

#=============================================================================
# Data Structures
=============================================================================#

"""
    ItacaResult

Structure to hold ITACA analysis results for a single event.
Combines RecoResult fields with blob/peak analysis.
"""
struct ItacaResult
    # Event identification
    event_id::Int
    thread_id::Int

    # Track data
    track::Petit.Tracks
    path::DataFrame
    track_length::Float64
    confidence::Float64

    # MC path
    mc_path::DataFrame

    # KDE data
    reco_s::Vector{Float64}
    kde_s::Vector{Float64}
    reco_kde_f::Vector{Float64}
    mc_kde_s::Vector{Float64}
    mc_kde_f::Vector{Float64}
    kde_bandwidth::Float64
    mc_kde_bandwidth::Float64

    # Extreme distances
    d1::Float64
    d2::Float64

    # Blob analysis
    Eb1::Float64
    Eb2::Float64
    asymmetry::Float64
    blob1_pos::NTuple{3,Float64}
    blob2_pos::NTuple{3,Float64}
    blob_radius::Float64

    # Peak analysis
    n_peaks::Int
    peak1_left::Float64
    peak1_right::Float64
    peak1_prom::Float64
    peak2_left::Float64
    peak2_right::Float64
    peak2_prom::Float64
end

#=============================================================================
# HDF5 Persistence (adapted from track_reco_mt.jl)
=============================================================================#

"""
    save_path_to_hdf5(path, group)

Save a raw path DataFrame to an HDF5 group.
"""
function save_path_to_hdf5(path::DataFrame, group)
    if nrow(path) == 0
        group["path_data"] = zeros(Float64, 0, 0)
        group["path_columns"] = String[]
        return
    end
    path_data = Matrix(path)
    group["path_data"] = path_data
    group["path_columns"] = String.(names(path))
end

"""
    save_mc_path_to_hdf5(mc_path, group)

Save a MC path DataFrame to an HDF5 group.
"""
function save_mc_path_to_hdf5(mc_path::DataFrame, group)
    if nrow(mc_path) == 0
        group["mc_path_data"] = zeros(Float64, 0, 0)
        group["mc_path_columns"] = String[]
        return
    end
    mc_data = Matrix(mc_path)
    group["mc_path_data"] = mc_data
    group["mc_path_columns"] = String.(names(mc_path))
end

"""
    save_itaca_results_to_hdf5(results, output_path, metadata)

Save ITACA analysis results to HDF5 file.
Extended from track_reco_mt.jl to include blob and peak data.
"""
function save_itaca_results_to_hdf5(results::Vector{ItacaResult}, output_path::String, metadata::Dict)
    h5open(output_path, "w") do fid
        # Save metadata as attributes
        for (key, val) in metadata
            attrs(fid)[key] = val
        end
        attrs(fid)["total_tracks_saved"] = length(results)

        if isempty(results)
            return
        end

        for (idx, result) in enumerate(results)
            track_group_name = "batch_1/track_$(idx)"
            g = create_group(fid, track_group_name)

            track = result.track

            # Event metadata
            g["event_id"] = result.event_id
            g["track_length"] = result.track_length
            g["confidence"] = result.confidence

            # Voxels data
            voxels_data = Matrix(track.voxels)
            g["voxels"] = voxels_data
            g["voxel_columns"] = String.(names(track.voxels))

            # Graph edges
            n_edges = ne(track.graph)
            if n_edges > 0
                edge_matrix = zeros(Int, n_edges, 2)
                for (i, edge) in enumerate(edges(track.graph))
                    edge_matrix[i, 1] = src(edge)
                    edge_matrix[i, 2] = dst(edge)
                end
                g["graph_edges"] = edge_matrix
            else
                g["graph_edges"] = zeros(Int, 0, 2)
            end
            g["n_vertices"] = nv(track.graph)

            # Components
            if !isempty(track.components)
                max_comp_len = maximum(length.(track.components))
                comp_matrix = zeros(Int, length(track.components), max_comp_len)
                for (i, comp) in enumerate(track.components)
                    for (j, v) in enumerate(comp)
                        comp_matrix[i, j] = v
                    end
                end
                g["components"] = comp_matrix
            else
                g["components"] = zeros(Int, 0, 0)
            end

            # Raw path
            save_path_to_hdf5(result.path, g)

            # MC path
            save_mc_path_to_hdf5(result.mc_path, g)

            # Projected voxel arc-lengths
            g["reco_s"] = result.reco_s

            # KDE results
            g["kde_s"] = result.kde_s
            g["reco_kde_f"] = result.reco_kde_f
            g["mc_kde_s"] = result.mc_kde_s
            g["mc_kde_f"] = result.mc_kde_f
            g["kde_bandwidth"] = result.kde_bandwidth
            g["mc_kde_bandwidth"] = result.mc_kde_bandwidth

            # Extreme distances
            g["d1"] = result.d1
            g["d2"] = result.d2

            # Blob analysis
            g["Eb1_keV"] = result.Eb1
            g["Eb2_keV"] = result.Eb2
            g["asymmetry"] = result.asymmetry
            g["blob1_pos"] = collect(result.blob1_pos)
            g["blob2_pos"] = collect(result.blob2_pos)
            g["blob_radius"] = result.blob_radius

            # Peak analysis
            g["n_peaks"] = result.n_peaks
            g["peak1_left"] = result.peak1_left
            g["peak1_right"] = result.peak1_right
            g["peak1_prom"] = result.peak1_prom
            g["peak2_left"] = result.peak2_left
            g["peak2_right"] = result.peak2_right
            g["peak2_prom"] = result.peak2_prom
        end
    end
end

#=============================================================================
# Event Processing (from itaca_single_track_analysis.jl)
=============================================================================#

"""
    process_single_event_itaca(hitsdf, event_id, thread_id, params)

Process a single event with full ITACA analysis pipeline.
Returns ItacaResult if successful, nothing otherwise.
"""
function process_single_event_itaca(hitsdf::DataFrame, event_id::Int, thread_id::Int;
                                    σt::Float64, σl::Float64,
                                    voxel_size::Float64, max_distance::Float64,
                                    mcvox_size::Float64,
                                    energy_threshold_keV::Float64,
                                    dfpars::Petit.DiffusionParams,
                                    kde_bandwidth::Float64,
                                    mc_kde_bandwidth::Float64,
                                    n_kde_eval::Int,
                                    nbins::Int, nsigma::Float64,
                                    Rb::Float64)

    # Step 1: Get event hits
    event_df = Petit.get_event(hitsdf, event_id)
    if nrow(event_df) == 0
        return nothing
    end

    # Step 2: Compute MC path
    mc_path = Petit.compute_mc_path(event_df, mcvox_size)
    if nrow(mc_path) == 0
        return nothing
    end

    # Step 3: Transform hits
    event_mc = Petit.transform_hits_df(event_df)

    # Step 4: Diffuse event
    diffused_df = Petit.diffuse_xyz_image_mc(event_mc;
                                              sigma_t_mm=σt,
                                              sigma_l_mm=σl,
                                              nbins=nbins,
                                              nsigma=nsigma)

    # Step 5: Voxelize
    voxels = Petit.voxelize_event(diffused_df, voxel_size)

    # Step 6: Make tracks
    tracks = Petit.make_tracks(voxels;
                               max_distance_mm=max_distance,
                               energy_threshold_kev=energy_threshold_keV,
                               diffusion=dfpars)

    # Step 7: Filter - only single-track events
    if length(tracks) != 1
        return nothing
    end

    track = tracks[1]

    # Step 8: Walk track from extremes
    walk_result = Petit.walk_track_from_extremes(track)
    if isnothing(walk_result.extremes[1]) || isempty(walk_result.path_indices)
        return nothing
    end

    # Step 9: Get raw path
    path = Petit.get_raw_path(track, walk_result.path_indices)
    if nrow(path) < 2
        return nothing
    end

    track_length = path.s[end]

    # Step 10: Compute extreme distances
    extreme_dists = Petit.compute_extreme_distances(path, mc_path)

    # Project voxels onto path for KDE
    reco_s = Petit.project_voxels_to_path(track.voxels, path)
    reco_E = Vector{Float64}(track.voxels.energy)

    # MC path data
    mc_s = Vector{Float64}(mc_path.s)
    mc_E = Vector{Float64}(mc_path.energy)
    mc_track_length = mc_path.s[end]

    # Step 11: RECO KDE
    kde_s = collect(range(0.0, track_length, length=n_kde_eval))
    reco_kde_f, _ = Petit.energy_weighted_kde(reco_s, reco_E, kde_s; bandwidth=kde_bandwidth)

    # Step 12: MC KDE
    mc_kde_s = collect(range(0.0, mc_track_length, length=n_kde_eval))
    mc_kde_f, _ = Petit.energy_weighted_kde(mc_s, mc_E, mc_kde_s; bandwidth=mc_kde_bandwidth)

    # Step 13: Find peaks
    peaks = Petit.find_peaks(reco_kde_f, kde_s; prom_scale=0.2)

    # Step 14: Extract peak1/peak2 info
    pk = kde_peaks(peaks, reco_kde_f)

    # Step 15: Blob analysis
    blobs = Petit.find_blob_energies(track, path; radius=Rb)

    return ItacaResult(
        event_id, thread_id,
        track, path, track_length, walk_result.confidence,
        mc_path,
        reco_s, kde_s, reco_kde_f, mc_kde_s, mc_kde_f,
        kde_bandwidth, mc_kde_bandwidth,
        extreme_dists.d1, extreme_dists.d2,
        blobs.Eb1, blobs.Eb2, blobs.asymmetry,
        (blobs.blob1.x, blobs.blob1.y, blobs.blob1.z),
        (blobs.blob2.x, blobs.blob2.y, blobs.blob2.z),
        Rb,
        length(peaks.indices),
        pk.peak1_left, pk.peak1_right, pk.peak1_prom,
        pk.peak2_left, pk.peak2_right, pk.peak2_prom
    )
end

#=============================================================================
# Multi-threaded Analysis Loop (from track_reco_mt.jl pattern)
=============================================================================#

"""
    analysis_loop_itaca_mt(hitsdf, thread_id; kwargs...)

Multi-threaded analysis loop for ITACA reconstruction.
Returns vector of ItacaResult.
"""
function analysis_loop_itaca_mt(hitsdf::DataFrame, thread_id::Int;
                                events_to_run=nothing,
                                initial_event::Int=1,
                                σt::Float64, σl::Float64,
                                voxel_size::Float64, max_distance::Float64,
                                mcvox_size::Float64,
                                energy_threshold_keV::Float64,
                                dfpars::Petit.DiffusionParams,
                                kde_bandwidth::Float64,
                                mc_kde_bandwidth::Float64,
                                n_kde_eval::Int,
                                nbins::Int, nsigma::Float64,
                                Rb::Float64,
                                print_level::String="quiet")

    event_ids = Petit.get_events_to_process(hitsdf, events_to_run, initial_event)
    n_single_track = 0
    n_events_processed = 0
    n_failed = 0
    results = ItacaResult[]

    is_verbose = (print_level == "verbose")
    is_very_verbose = (print_level == "very_verbose")

    for event_id in event_ids
        n_events_processed += 1

        if is_very_verbose
            println("Thread $thread_id: Processing event $event_id")
        end

        try
            result = process_single_event_itaca(hitsdf, event_id, thread_id;
                                                σt=σt, σl=σl,
                                                voxel_size=voxel_size,
                                                max_distance=max_distance,
                                                mcvox_size=mcvox_size,
                                                energy_threshold_keV=energy_threshold_keV,
                                                dfpars=dfpars,
                                                kde_bandwidth=kde_bandwidth,
                                                mc_kde_bandwidth=mc_kde_bandwidth,
                                                n_kde_eval=n_kde_eval,
                                                nbins=nbins, nsigma=nsigma,
                                                Rb=Rb)
            if !isnothing(result)
                push!(results, result)
                n_single_track += 1
            end
        catch e
            if is_verbose || is_very_verbose
                println("Thread $thread_id: Error processing event $event_id: $e")
            end
            n_failed += 1
        end
    end

    if print_level != "quiet"
        println("Thread $thread_id: Completed! Processed $n_events_processed events, $n_single_track single-track, $n_failed failed.")
    end

    return results
end

#=============================================================================
# Main Entry Point (from track_reco_mt.jl pattern)
=============================================================================#

"""
    event_loop_itaca_mt(cmdir; kwargs...)

Multi-threaded event processing for ITACA analysis.

Produces:
- CSV file with analysis summary
- HDF5 files containing reconstructed tracks with blob/peak data
"""
function event_loop_itaca_mt(cmdir::String;
                             input_file::String,
                             outdir::String="results",
                             outbase::String="itaca",
                             ievt::Int=1,
                             levt::Int=-1,
                             nthreads::Int=1,
                             particle_type::String="ion",
                             ldrft::Float64=100.0,
                             tK::Float64=297.0,
                             edrift::Float64=500.0,
                             Pbar::Float64=15.0,
                             dt::Float64=3.5,
                             dl::Float64=0.9,
                             energy_threshold_ions::Float64=10.0,
                             nbins::Int=100,
                             nsigma::Float64=3.0,
                             n_kde_eval::Int=200,
                             Rb::Float64=10.0,
                             print_level::String="quiet")

    # Compute diffusion parameters (from itaca_single_track_analysis.jl)
    σt, σl = get_sigma(particle_type, ldrft; dt=dt, dl=dl, tK=tK, edrift=edrift, Pbar=Pbar)
    voxel_size, mcvox_size, max_distance = get_voxel_size_and_distance(ldrft, σt)
    eth = get_energy_threshold(particle_type; energy_threshold_ions=energy_threshold_ions)

    kde_bandwidth = 2.0 * voxel_size
    mc_kde_bandwidth = 2.0 * mcvox_size

    dfpars = Petit.DiffusionParams(ldrft, σt, σl, voxel_size, max_distance, eth, nbins, nsigma)

    # Load and validate input (from track_reco_mt.jl)
    hitsdf, ntot, nevents_config, last_event = Petit.load_and_validate_input(cmdir, input_file, ievt, levt)
    nevents_to_process = last_event - ievt + 1
    optimal_nthreads = Petit.get_optimal_threads(nthreads)
    thread_ranges = Petit.split_events_for_threads(nevents_to_process, optimal_nthreads, ievt)

    # Create output directory
    mkpath(outdir)

    # Print configuration
    println("\n" * "="^70)
    println("ITACA TRACK RECONSTRUCTION + ANALYSIS (MT)")
    println("="^70)
    println("Configuration:")
    println("  Input file:         $input_file")
    println("  Input directory:    $cmdir")
    println("  Output directory:   $outdir")
    println("  Output base name:   $outbase")
    println("  Events:             $ievt to $last_event ($nevents_to_process total)")
    println("  Threads:            $optimal_nthreads")
    println("  Particle type:      $particle_type")
    println("  Drift length:       $ldrft cm")
    println("  σ_t:                $(round(σt, digits=3)) mm")
    println("  σ_l:                $(round(σl, digits=3)) mm")
    println("  Voxel size:         $(round(voxel_size, digits=3)) mm")
    println("  Max distance:       $(round(max_distance, digits=3)) mm")
    println("  KDE bandwidth:      $(round(kde_bandwidth, digits=3)) mm")
    println("  MC KDE bandwidth:   $(round(mc_kde_bandwidth, digits=3)) mm")
    println("  Energy threshold:   $(round(eth, digits=3)) keV")
    println("  Blob radius:        $(round(Rb, digits=3)) mm")
    println("  Print level:        $print_level")
    println("="^70)

    # Process in parallel
    println("\nStarting multi-threaded processing...")
    all_results = Vector{Any}(undef, length(thread_ranges))
    all_itaca_results = Vector{Vector{ItacaResult}}(undef, length(thread_ranges))

    Threads.@threads for i in 1:length(thread_ranges)
        start_evt, num_evts = thread_ranges[i]

        results = analysis_loop_itaca_mt(hitsdf, i;
                                         events_to_run=num_evts,
                                         initial_event=start_evt,
                                         σt=σt, σl=σl,
                                         voxel_size=voxel_size,
                                         max_distance=max_distance,
                                         mcvox_size=mcvox_size,
                                         energy_threshold_keV=eth,
                                         dfpars=dfpars,
                                         kde_bandwidth=kde_bandwidth,
                                         mc_kde_bandwidth=mc_kde_bandwidth,
                                         n_kde_eval=n_kde_eval,
                                         nbins=nbins, nsigma=nsigma,
                                         Rb=Rb,
                                         print_level=print_level)

        all_itaca_results[i] = results

        # Save HDF5 for this thread
        output_file = "$(outbase)_th_$i.h5"
        output_path = joinpath(outdir, output_file)

        metadata = Dict(
            "input_file" => input_file,
            "thread_id" => i,
            "first_event_processed" => start_evt,
            "last_event_processed" => start_evt + num_evts - 1,
            "events_processed" => num_evts,
            "particle_type" => particle_type,
            "ldrift_cm" => ldrft,
            "sigma_t_mm" => σt,
            "sigma_l_mm" => σl,
            "voxel_size_mm" => voxel_size,
            "max_distance_mm" => max_distance,
            "kde_bandwidth_mm" => kde_bandwidth,
            "mc_kde_bandwidth_mm" => mc_kde_bandwidth,
            "n_kde_eval" => n_kde_eval,
            "energy_threshold_kev" => eth,
            "blob_radius_mm" => Rb,
            "nbins" => nbins,
            "nsigma" => nsigma
        )

        save_itaca_results_to_hdf5(results, output_path, metadata)

        all_results[i] = Dict(
            "thread_id" => i,
            "output_file" => output_path,
            "events_processed" => num_evts,
            "tracks_saved" => length(results),
            "first_event" => start_evt,
            "last_event" => start_evt + num_evts - 1
        )

        println("Thread $i: Saved $(length(results)) tracks to $output_file")
    end

    # Merge CSV results from all threads
    println("\nMerging results...")

    csv_df = DataFrame(
        event = Int[],
        thread_id = Int[],
        track_length_mm = Float64[],
        confidence = Float64[],
        Eb1_keV = Float64[],
        Eb2_keV = Float64[],
        asymmetry = Float64[],
        d1_mm = Float64[],
        d2_mm = Float64[],
        n_peaks = Int[],
        peak1_left = Float64[],
        peak1_right = Float64[],
        peak1_prom = Float64[],
        peak2_left = Float64[],
        peak2_right = Float64[],
        peak2_prom = Float64[]
    )

    for thread_results in all_itaca_results
        for r in thread_results
            push!(csv_df, (
                event = r.event_id,
                thread_id = r.thread_id,
                track_length_mm = r.track_length,
                confidence = r.confidence,
                Eb1_keV = r.Eb1,
                Eb2_keV = r.Eb2,
                asymmetry = r.asymmetry,
                d1_mm = r.d1,
                d2_mm = r.d2,
                n_peaks = r.n_peaks,
                peak1_left = r.peak1_left,
                peak1_right = r.peak1_right,
                peak1_prom = r.peak1_prom,
                peak2_left = r.peak2_left,
                peak2_right = r.peak2_right,
                peak2_prom = r.peak2_prom
            ))
        end
    end

    # Sort by event number
    sort!(csv_df, :event)

    # Save CSV
    csv_path = joinpath(outdir, "$(outbase)_analysis.csv")
    CSV.write(csv_path, csv_df)
    println("CSV saved to: $csv_path")

    # Save metadata
    total_tracks = sum(r["tracks_saved"] for r in all_results)

    metadata_df = DataFrame(
        parameter = String[],
        value = String[]
    )

    push!(metadata_df, (parameter="input_file", value=input_file))
    push!(metadata_df, (parameter="particle_type", value=particle_type))
    push!(metadata_df, (parameter="ldrft_cm", value=string(ldrft)))
    push!(metadata_df, (parameter="sigma_t_mm", value=string(round(σt, digits=4))))
    push!(metadata_df, (parameter="sigma_l_mm", value=string(round(σl, digits=4))))
    push!(metadata_df, (parameter="voxel_size_mm", value=string(round(voxel_size, digits=4))))
    push!(metadata_df, (parameter="max_distance_mm", value=string(round(max_distance, digits=4))))
    push!(metadata_df, (parameter="energy_threshold_keV", value=string(round(eth, digits=4))))
    push!(metadata_df, (parameter="kde_bandwidth_mm", value=string(round(kde_bandwidth, digits=4))))
    push!(metadata_df, (parameter="blob_radius_mm", value=string(Rb)))
    push!(metadata_df, (parameter="n_events_processed", value=string(nevents_to_process)))
    push!(metadata_df, (parameter="n_single_track", value=string(total_tracks)))
    push!(metadata_df, (parameter="nthreads", value=string(optimal_nthreads)))
    push!(metadata_df, (parameter="timestamp", value=string(now())))

    metadata_path = joinpath(outdir, "$(outbase)_metadata.csv")
    CSV.write(metadata_path, metadata_df)
    println("Metadata saved to: $metadata_path")

    # Print summary
    println("\n" * "="^70)
    println("PROCESSING COMPLETE")
    println("="^70)
    println("Total events processed: $nevents_to_process")
    println("Total single-track events: $total_tracks")
    println("Efficiency: $(round(100 * total_tracks / nevents_to_process, digits=1))%")
    println("\nOutput files:")
    println("  CSV: $csv_path")
    println("  Metadata: $metadata_path")
    println("  HDF5 files:")
    for r in all_results
        println("    $(r["output_file"]) ($(r["tracks_saved"]) tracks)")
    end
    println("="^70)

    return Dict(
        "thread_results" => all_results,
        "total_tracks" => total_tracks,
        "events_processed" => nevents_to_process,
        "csv_path" => csv_path,
        "metadata_path" => metadata_path
    )
end

#=============================================================================
# Command Line Interface
=============================================================================#

"""
    main()

Main function for command-line usage.
"""
function main()
    # Check minimum required arguments
    if length(ARGS) < 2
        println("Error: Missing required arguments")
        println("\nUsage: julia -t <nthreads> itaca_track_reco_mt.jl <cmdir> <input_file> [options]")
        println("\nRequired arguments:")
        println("  cmdir           Directory containing the input file")
        println("  input_file      Name of the HDF5 input file")
        println("\nOptional arguments:")
        println("  --outdir=PATH         Output directory (default: results)")
        println("  --outbase=NAME        Base name for output files (default: itaca)")
        println("  --ievt=N              First event to process (default: 1)")
        println("  --levt=N              Last event to process (default: -1, all)")
        println("  --nthreads=N          Number of threads to use (default: 1)")
        println("  --particle=TYPE       Particle type: 'ion' or 'electron' (default: ion)")
        println("  --ldrft=X             Drift length in cm (default: 100.0)")
        println("  --tK=X                Temperature in K (default: 297.0)")
        println("  --edrift=X            Drift field in V/cm (default: 500.0)")
        println("  --Pbar=X              Pressure in bar (default: 15.0)")
        println("  --dt=X                Transverse diffusion coeff mm/√cm (default: 3.5)")
        println("  --dl=X                Longitudinal diffusion coeff mm/√cm (default: 0.9)")
        println("  --eth-ion=X           Energy threshold for ions (default: 10.0)")
        println("  --nbins=N             Diffusion histogram bins (default: 100)")
        println("  --nsigma=X            Diffusion extent in sigma (default: 3.0)")
        println("  --nkde=N              Number of KDE evaluation points (default: 200)")
        println("  --Rb=X                Blob radius in mm (default: 10.0)")
        println("  --print_level=LEVEL   quiet, verbose, very_verbose (default: quiet)")
        println("\nExample:")
        println("  julia -t 8 itaca_track_reco_mt.jl /data/HD5t/itaca/ bb0nu.h5 --outdir=results/run01 --outbase=ion_100cm --nthreads=8")
        exit(1)
    end

    # Parse required arguments
    cmdir = ARGS[1]
    input_file = ARGS[2]

    # Default values
    outdir = "results"
    outbase = "itaca"
    ievt = 1
    levt = -1
    nthreads = 1
    particle_type = "ion"
    ldrft = 100.0
    tK = 297.0
    edrift = 500.0
    Pbar = 15.0
    dt = 3.5
    dl = 0.9
    energy_threshold_ions = 10.0
    nbins = 100
    nsigma = 3.0
    n_kde_eval = 200
    Rb = 10.0
    print_level = "quiet"

    # Parse optional arguments
    for arg in ARGS[3:end]
        if startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) != 2
                println("Warning: Ignoring malformed argument: $arg")
                continue
            end
            key, value = parts

            try
                if key == "outdir"
                    outdir = String(value)
                elseif key == "outbase"
                    outbase = String(value)
                elseif key == "ievt"
                    ievt = parse(Int, value)
                elseif key == "levt"
                    levt = parse(Int, value)
                elseif key == "nthreads"
                    nthreads = parse(Int, value)
                elseif key == "particle"
                    particle_type = String(value)
                elseif key == "ldrft"
                    ldrft = parse(Float64, value)
                elseif key == "tK"
                    tK = parse(Float64, value)
                elseif key == "edrift"
                    edrift = parse(Float64, value)
                elseif key == "Pbar"
                    Pbar = parse(Float64, value)
                elseif key == "dt"
                    dt = parse(Float64, value)
                elseif key == "dl"
                    dl = parse(Float64, value)
                elseif key == "eth-ion"
                    energy_threshold_ions = parse(Float64, value)
                elseif key == "nbins"
                    nbins = parse(Int, value)
                elseif key == "nsigma"
                    nsigma = parse(Float64, value)
                elseif key == "nkde"
                    n_kde_eval = parse(Int, value)
                elseif key == "Rb"
                    Rb = parse(Float64, value)
                elseif key == "print_level"
                    print_level = String(value)
                else
                    println("Warning: Unknown argument: --$key")
                end
            catch e
                println("Error parsing argument --$key=$value: $e")
                exit(1)
            end
        end
    end

    # Run
    event_loop_itaca_mt(cmdir;
                        input_file=input_file,
                        outdir=outdir,
                        outbase=outbase,
                        ievt=ievt,
                        levt=levt,
                        nthreads=nthreads,
                        particle_type=particle_type,
                        ldrft=ldrft,
                        tK=tK,
                        edrift=edrift,
                        Pbar=Pbar,
                        dt=dt,
                        dl=dl,
                        energy_threshold_ions=energy_threshold_ions,
                        nbins=nbins,
                        nsigma=nsigma,
                        n_kde_eval=n_kde_eval,
                        Rb=Rb,
                        print_level=print_level)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
