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

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit
# Helper functions (get_sigma, get_voxel_size_and_distance, get_energy_threshold)
# are now part of Petit module via itaca_functions.jl

#=============================================================================
# Helper Functions
=============================================================================#

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

#=============================================================================
# Data Structures
=============================================================================#

"""
    ItacaResult

Structure to hold ITACA analysis results for a single event.
Contains results for 4 analysis categories:
- MC (default/KDT method)
- MC-MST (MST method on MC track)
- RECO (default/KDT method)
- RECO-MST (MST method on RECO track)
"""
struct ItacaResult
    # Event identification
    event_id::Int
    thread_id::Int

    # Per-event diffusion parameters
    ldrft_cm::Float64
    sigma_t_mm::Float64
    sigma_l_mm::Float64
    voxel_size_mm::Float64
    max_distance_mm::Float64

    # MC track and paths
    mc_track::Petit.Tracks              # MC reconstructed track
    mc_path_true::DataFrame             # True MC path (from compute_mc_path)
    mc_path::DataFrame                  # MC reconstructed path (default/KDT)
    mc_path_mst::DataFrame              # MC reconstructed path (MST)
    mc_track_length::Float64            # MC track length (default)
    mc_track_length_mst::Float64        # MC track length (MST)

    # RECO track and paths
    track::Petit.Tracks                 # RECO track (default)
    path::DataFrame                     # RECO path (default/KDT)
    path_mst::DataFrame                 # RECO path (MST)
    track_length::Float64               # RECO track length (default)
    track_length_mst::Float64           # RECO track length (MST)

    # MC blob analysis (default/KDT method)
    Eb1_mc::Float64
    Eb2_mc::Float64
    asymmetry_mc::Float64
    blob1_pos_mc::NTuple{3,Float64}
    blob2_pos_mc::NTuple{3,Float64}

    # MC blob analysis (MST method)
    Eb1_mc_mst::Float64
    Eb2_mc_mst::Float64
    asymmetry_mc_mst::Float64
    blob1_pos_mc_mst::NTuple{3,Float64}
    blob2_pos_mc_mst::NTuple{3,Float64}

    # RECO blob analysis (default/KDT method)
    Eb1::Float64
    Eb2::Float64
    asymmetry::Float64
    blob1_pos::NTuple{3,Float64}
    blob2_pos::NTuple{3,Float64}

    # RECO blob analysis (MST method)
    Eb1_mst::Float64
    Eb2_mst::Float64
    asymmetry_mst::Float64
    blob1_pos_mst::NTuple{3,Float64}
    blob2_pos_mst::NTuple{3,Float64}

    # MC extreme distances (default/KDT method)
    d1_mc::Float64
    d2_mc::Float64

    # MC extreme distances (MST method)
    d1_mc_mst::Float64
    d2_mc_mst::Float64

    # RECO extreme distances (default/KDT method)
    d1::Float64
    d2::Float64

    # RECO extreme distances (MST method)
    d1_mst::Float64
    d2_mst::Float64

    # Blob radius used for analysis
    blob_radius::Float64
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
    save_track_to_hdf5(track, group, prefix)

Save a Tracks object to an HDF5 group with given prefix.
"""
function save_track_to_hdf5(track::Petit.Tracks, group, prefix::String)
    # Voxels data
    voxels_data = Matrix(track.voxels)
    group["$(prefix)_voxels"] = voxels_data
    group["$(prefix)_voxel_columns"] = String.(names(track.voxels))

    # Graph edges
    n_edges = ne(track.graph)
    if n_edges > 0
        edge_matrix = zeros(Int, n_edges, 2)
        for (i, edge) in enumerate(edges(track.graph))
            edge_matrix[i, 1] = src(edge)
            edge_matrix[i, 2] = dst(edge)
        end
        group["$(prefix)_graph_edges"] = edge_matrix
    else
        group["$(prefix)_graph_edges"] = zeros(Int, 0, 2)
    end
    group["$(prefix)_n_vertices"] = nv(track.graph)

    # Components
    if !isempty(track.components)
        max_comp_len = maximum(length.(track.components))
        comp_matrix = zeros(Int, length(track.components), max_comp_len)
        for (i, comp) in enumerate(track.components)
            for (j, v) in enumerate(comp)
                comp_matrix[i, j] = v
            end
        end
        group["$(prefix)_components"] = comp_matrix
    else
        group["$(prefix)_components"] = zeros(Int, 0, 0)
    end
end

"""
    save_path_to_hdf5_named(path, group, name)

Save a path DataFrame to an HDF5 group with given name.
"""
function save_path_to_hdf5_named(path::DataFrame, group, name::String)
    if nrow(path) == 0
        group["$(name)_data"] = zeros(Float64, 0, 0)
        group["$(name)_columns"] = String[]
        return
    end
    path_data = Matrix(path)
    group["$(name)_data"] = path_data
    group["$(name)_columns"] = String.(names(path))
end

"""
    save_itaca_results_to_hdf5(results, output_path, metadata)

Save ITACA analysis results to HDF5 file.
Saves all 4 analysis categories: MC-default, MC-MST, RECO-default, RECO-MST.
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

            # Event metadata
            g["event_id"] = result.event_id

            # Per-event diffusion parameters
            g["ldrft_cm"] = result.ldrft_cm
            g["sigma_t_mm"] = result.sigma_t_mm
            g["sigma_l_mm"] = result.sigma_l_mm
            g["voxel_size_mm"] = result.voxel_size_mm
            g["max_distance_mm"] = result.max_distance_mm

            # MC track
            save_track_to_hdf5(result.mc_track, g, "mc_track")

            # RECO track
            save_track_to_hdf5(result.track, g, "reco_track")

            # Paths
            save_path_to_hdf5_named(result.mc_path_true, g, "mc_path_true")
            save_path_to_hdf5_named(result.mc_path, g, "mc_path")
            save_path_to_hdf5_named(result.mc_path_mst, g, "mc_path_mst")
            save_path_to_hdf5_named(result.path, g, "reco_path")
            save_path_to_hdf5_named(result.path_mst, g, "reco_path_mst")

            # Track lengths
            g["mc_track_length"] = result.mc_track_length
            g["mc_track_length_mst"] = result.mc_track_length_mst
            g["track_length"] = result.track_length
            g["track_length_mst"] = result.track_length_mst

            # MC blob analysis (default)
            g["Eb1_mc_keV"] = result.Eb1_mc
            g["Eb2_mc_keV"] = result.Eb2_mc
            g["asymmetry_mc"] = result.asymmetry_mc
            g["blob1_pos_mc"] = collect(result.blob1_pos_mc)
            g["blob2_pos_mc"] = collect(result.blob2_pos_mc)

            # MC blob analysis (MST)
            g["Eb1_mc_mst_keV"] = result.Eb1_mc_mst
            g["Eb2_mc_mst_keV"] = result.Eb2_mc_mst
            g["asymmetry_mc_mst"] = result.asymmetry_mc_mst
            g["blob1_pos_mc_mst"] = collect(result.blob1_pos_mc_mst)
            g["blob2_pos_mc_mst"] = collect(result.blob2_pos_mc_mst)

            # RECO blob analysis (default)
            g["Eb1_keV"] = result.Eb1
            g["Eb2_keV"] = result.Eb2
            g["asymmetry"] = result.asymmetry
            g["blob1_pos"] = collect(result.blob1_pos)
            g["blob2_pos"] = collect(result.blob2_pos)

            # RECO blob analysis (MST)
            g["Eb1_mst_keV"] = result.Eb1_mst
            g["Eb2_mst_keV"] = result.Eb2_mst
            g["asymmetry_mst"] = result.asymmetry_mst
            g["blob1_pos_mst"] = collect(result.blob1_pos_mst)
            g["blob2_pos_mst"] = collect(result.blob2_pos_mst)

            # MC extreme distances (default)
            g["d1_mc"] = result.d1_mc
            g["d2_mc"] = result.d2_mc

            # MC extreme distances (MST)
            g["d1_mc_mst"] = result.d1_mc_mst
            g["d2_mc_mst"] = result.d2_mc_mst

            # RECO extreme distances (default)
            g["d1"] = result.d1
            g["d2"] = result.d2

            # RECO extreme distances (MST)
            g["d1_mst"] = result.d1_mst
            g["d2_mst"] = result.d2_mst

            # Blob radius
            g["blob_radius"] = result.blob_radius
        end
    end
end

#=============================================================================
# Event Processing (from itaca_single_track_analysis.jl)
=============================================================================#

"""
    process_single_event_itaca(hitsdf, event_id, thread_id, params)

Process a single event with full ITACA analysis pipeline.
Performs 4 analyses: MC-default, MC-MST, RECO-default, RECO-MST.
Returns ItacaResult if successful, nothing otherwise.
"""
function process_single_event_itaca(hitsdf::DataFrame, event_id::Int, thread_id::Int;
                                    ldrft_cm::Float64,
                                    σt::Float64, σl::Float64,
                                    voxel_size::Float64, max_distance::Float64,
                                    mcvox_size::Float64,
                                    energy_threshold_keV::Float64,
                                    dfpars::Petit.DiffusionParams,
                                    nbins::Int, nsigma::Float64,
                                    Rb::Float64,
                                    is_double_beta::Bool=false,
                                    graph_method::String="KDT",
                                    knn_k::Int=10,
                                    use_energy_weighting::Bool=true,
                                    use_edge_energy_weighting::Bool=true,
                                    use_mst_fallback::Bool=false)

    # Step 1: Get event hits
    event_df = Petit.get_event(hitsdf, event_id)
    if nrow(event_df) == 0
        return nothing
    end

    # Step 2: Compute true MC path
    mc_path_true = Petit.compute_mc_path(event_df, mcvox_size; is_double_beta=is_double_beta)
    if nrow(mc_path_true) == 0
        return nothing
    end

    #=========================================================================
    # MC TRACK RECONSTRUCTION
    =========================================================================#

    # Voxelize MC hits (without diffusion)
    mcvoxels = Petit.voxelize_event(event_df, mcvox_size)

    # Make MC tracks
    mctracks = Petit.make_tracks(mcvoxels;
                                 max_distance_mm=2.5 * mcvox_size,
                                 energy_threshold_kev=1.0,
                                 diffusion=dfpars)

    # Filter - need at least one MC track
    if length(mctracks) < 1
        return nothing
    end

    mctrack = mctracks[1]

    # MC Walk (default/KDT method)
    mcwalk = Petit.walk_track_from_extremes(mctrack;
        use_energy_weighting=use_energy_weighting,
        use_edge_energy_weighting=use_edge_energy_weighting)
    if isnothing(mcwalk.extremes[1]) || isempty(mcwalk.path_indices)
        return nothing
    end

    # MC path (default)
    mc_path = Petit.get_raw_path(mctrack, mcwalk.path_indices)
    if nrow(mc_path) < 2
        return nothing
    end
    mc_track_length = mc_path.s[end]

    # MC extreme distances (default vs true MC)
    xextreme_dists = Petit.compute_extreme_distances(mc_path, mc_path_true)

    # MC blob analysis (default)
    blobs_mc = Petit.find_blob_energies(mctrack, mc_path; radius=Rb)

    # MC Walk (MST method)
    mccoords = Petit.extract_coords(mctrack)
    mstmc = Petit.find_extremes_mst_diameter(mctrack, mccoords)
    mcwalk_mst = get_walk_mst(mstmc, mctrack)

    # MC path (MST)
    mc_path_mst = Petit.get_raw_path(mctrack, mcwalk_mst.path_indices)
    mc_track_length_mst = nrow(mc_path_mst) >= 2 ? mc_path_mst.s[end] : 0.0

    # MC extreme distances (MST vs true MC)
    xextreme_dists_mst = Petit.compute_extreme_distances(mc_path_mst, mc_path_true)

    # MC blob analysis (MST)
    blobs_mc_mst = Petit.find_blob_energies(mctrack, mc_path_mst; radius=Rb)

    #=========================================================================
    # RECO TRACK RECONSTRUCTION
    =========================================================================#

    # Transform hits
    event_mc = Petit.transform_hits_df(event_df)

    # Diffuse event
    diffused_df = Petit.diffuse_xyz_image_mc(event_mc;
                                              sigma_t_mm=σt,
                                              sigma_l_mm=σl,
                                              nbins=nbins,
                                              nsigma=nsigma)

    # Voxelize
    voxels = Petit.voxelize_event(diffused_df, voxel_size)

    # Make RECO tracks
    tracks = Petit.make_tracks(voxels;
                               max_distance_mm=max_distance,
                               energy_threshold_kev=energy_threshold_keV,
                               diffusion=dfpars,
                               method=graph_method,
                               k=knn_k)

    # Filter - only single-track events
    if length(tracks) != 1
        return nothing
    end

    track = tracks[1]

    # RECO Walk (default/KDT method)
    walk_result = Petit.walk_track_from_extremes(track;
        use_energy_weighting=use_energy_weighting,
        use_edge_energy_weighting=use_edge_energy_weighting,
        use_mst_fallback=use_mst_fallback)
    if isnothing(walk_result.extremes[1]) || isempty(walk_result.path_indices)
        return nothing
    end

    # RECO path (default)
    path = Petit.get_raw_path(track, walk_result.path_indices)
    if nrow(path) < 2
        return nothing
    end
    track_length = path.s[end]

    # RECO extreme distances (default vs true MC)
    extreme_dists = Petit.compute_extreme_distances(path, mc_path_true)

    # RECO blob analysis (default)
    blobs = Petit.find_blob_energies(track, path; radius=Rb)

    # RECO Walk (MST method)
    coords = Petit.extract_coords(track)
    mst = Petit.find_extremes_mst_diameter(track, coords)
    walk_mst = get_walk_mst(mst, track)

    # RECO path (MST)
    path_mst = Petit.get_raw_path(track, walk_mst.path_indices)
    track_length_mst = nrow(path_mst) >= 2 ? path_mst.s[end] : 0.0

    # RECO extreme distances (MST vs true MC)
    extreme_dists_mst = Petit.compute_extreme_distances(path_mst, mc_path_true)

    # RECO blob analysis (MST)
    blobs_mst = Petit.find_blob_energies(track, path_mst; radius=Rb)

    #=========================================================================
    # BUILD RESULT
    =========================================================================#

    return ItacaResult(
        event_id, thread_id,
        # Per-event diffusion parameters
        ldrft_cm, σt, σl, voxel_size, max_distance,
        # MC track and paths
        mctrack, mc_path_true, mc_path, mc_path_mst,
        mc_track_length, mc_track_length_mst,
        # RECO track and paths
        track, path, path_mst,
        track_length, track_length_mst,
        # MC blob analysis (default)
        blobs_mc.Eb1, blobs_mc.Eb2, blobs_mc.asymmetry,
        (blobs_mc.blob1.x, blobs_mc.blob1.y, blobs_mc.blob1.z),
        (blobs_mc.blob2.x, blobs_mc.blob2.y, blobs_mc.blob2.z),
        # MC blob analysis (MST)
        blobs_mc_mst.Eb1, blobs_mc_mst.Eb2, blobs_mc_mst.asymmetry,
        (blobs_mc_mst.blob1.x, blobs_mc_mst.blob1.y, blobs_mc_mst.blob1.z),
        (blobs_mc_mst.blob2.x, blobs_mc_mst.blob2.y, blobs_mc_mst.blob2.z),
        # RECO blob analysis (default)
        blobs.Eb1, blobs.Eb2, blobs.asymmetry,
        (blobs.blob1.x, blobs.blob1.y, blobs.blob1.z),
        (blobs.blob2.x, blobs.blob2.y, blobs.blob2.z),
        # RECO blob analysis (MST)
        blobs_mst.Eb1, blobs_mst.Eb2, blobs_mst.asymmetry,
        (blobs_mst.blob1.x, blobs_mst.blob1.y, blobs_mst.blob1.z),
        (blobs_mst.blob2.x, blobs_mst.blob2.y, blobs_mst.blob2.z),
        # MC extreme distances (default)
        xextreme_dists.d1, xextreme_dists.d2,
        # MC extreme distances (MST)
        xextreme_dists_mst.d1, xextreme_dists_mst.d2,
        # RECO extreme distances (default)
        extreme_dists.d1, extreme_dists.d2,
        # RECO extreme distances (MST)
        extreme_dists_mst.d1, extreme_dists_mst.d2,
        # Blob radius
        Rb
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
                                # Drift length range
                                lmin::Float64=100.0,
                                lmax::Float64=100.0,
                                use_random_ldrft::Bool=false,
                                # Fixed diffusion parameters (used when lmin == lmax)
                                ldrft_fixed::Float64=100.0,
                                σt_fixed::Float64=0.0, σl_fixed::Float64=0.0,
                                voxel_size_fixed::Float64=0.0, max_distance_fixed::Float64=0.0,
                                mcvox_size_fixed::Float64=0.0,
                                eth_fixed::Float64=10.0,
                                dfpars_fixed::Petit.DiffusionParams=Petit.DiffusionParams(0,0,0,0,0,0,0,0),
                                # Common parameters
                                nbins::Int=100, nsigma::Float64=3.0,
                                Rb::Float64=10.0,
                                # Parameters needed for recomputing diffusion
                                particle_type::String="ion",
                                dt::Float64=3.5, dl::Float64=0.9,
                                tK::Float64=297.0, edrift::Float64=500.0, Pbar::Float64=15.0,
                                energy_threshold_ions::Float64=10.0,
                                is_double_beta::Bool=false,
                                print_level::String="quiet",
                                # Graph construction method
                                graph_method::String="KDT",
                                knn_k::Int=10,
                                # Energy weighting options for walk_track_from_extremes
                                use_energy_weighting::Bool=true,
                                use_edge_energy_weighting::Bool=true,
                                use_mst_fallback::Bool=false)

    event_ids = Petit.get_events_to_process(hitsdf, events_to_run, initial_event)
    n_single_track = 0
    n_events_processed = 0
    n_failed = 0
    results = ItacaResult[]

    is_verbose = (print_level == "verbose")
    is_very_verbose = (print_level == "very_verbose")

    for event_id in event_ids
        n_events_processed += 1

        # Determine drift length and diffusion parameters for this event
        if use_random_ldrft
            # Sample random drift length and compute diffusion parameters
            event_ldrft = lmin + rand() * (lmax - lmin)
            event_σt, event_σl = Petit.get_sigma(particle_type, event_ldrft;
                                           dt=dt, dl=dl, tK=tK, edrift=edrift, Pbar=Pbar)
            event_voxel_size, event_mcvox_size, event_max_distance = Petit.get_voxel_size_and_distance(event_ldrft, event_σt)
            event_eth = Petit.get_energy_threshold(particle_type; energy_threshold_ions=energy_threshold_ions)
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
            event_dfpars = dfpars_fixed
        end

        # Print event header based on print_level
        if is_very_verbose
            println("Thread $thread_id: Event $event_id")
            println("  ldrft=$(round(event_ldrft, digits=2))cm, σt=$(round(event_σt, digits=3))mm, σl=$(round(event_σl, digits=3))mm, vox=$(round(event_voxel_size, digits=3))mm")
        elseif is_verbose
            print("T$thread_id E$event_id: l=$(round(event_ldrft, digits=1))cm, σt=$(round(event_σt, digits=2))mm ... ")
        end

        try
            result = process_single_event_itaca(hitsdf, event_id, thread_id;
                                                ldrft_cm=event_ldrft,
                                                σt=event_σt, σl=event_σl,
                                                voxel_size=event_voxel_size,
                                                max_distance=event_max_distance,
                                                mcvox_size=event_mcvox_size,
                                                energy_threshold_keV=event_eth,
                                                dfpars=event_dfpars,
                                                nbins=nbins, nsigma=nsigma,
                                                Rb=Rb,
                                                is_double_beta=is_double_beta,
                                                graph_method=graph_method,
                                                knn_k=knn_k,
                                                use_energy_weighting=use_energy_weighting,
                                                use_edge_energy_weighting=use_edge_energy_weighting,
                                                use_mst_fallback=use_mst_fallback)
            if !isnothing(result)
                push!(results, result)
                n_single_track += 1
                if is_verbose
                    println("1 track")
                elseif is_very_verbose
                    println("  → Single track found")
                end
            else
                if is_verbose
                    println("rejected")
                elseif is_very_verbose
                    println("  → Rejected (multi-track or no track)")
                end
            end
        catch e
            if is_verbose
                println("error: $e")
            elseif is_very_verbose
                println("  → Error: $e")
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
                             lmin::Float64=100.0,
                             lmax::Float64=100.0,
                             tK::Float64=297.0,
                             edrift::Float64=500.0,
                             Pbar::Float64=15.0,
                             dt::Float64=3.5,
                             dl::Float64=0.9,
                             energy_threshold_ions::Float64=10.0,
                             nbins::Int=100,
                             nsigma::Float64=3.0,
                             Rb::Float64=10.0,
                             print_level::String="quiet",
                             graph_method::String="KDT",
                             knn_k::Int=10,
                             use_energy_weighting::Bool=true,
                             use_edge_energy_weighting::Bool=true,
                             use_mst_fallback::Bool=false)

    # Validate lmin/lmax
    if lmin > lmax
        error("lmin ($lmin) must be less than or equal to lmax ($lmax)")
    end
    use_random_ldrft = (lmin < lmax)

    # Compute diffusion parameters for fixed drift case (lmin == lmax)
    # For random drift, these will be recomputed per event
    ldrft_fixed = lmin  # lmin == lmax in fixed mode
    σt, σl = Petit.get_sigma(particle_type, ldrft_fixed; dt=dt, dl=dl, tK=tK, edrift=edrift, Pbar=Pbar)
    voxel_size, mcvox_size, max_distance = Petit.get_voxel_size_and_distance(ldrft_fixed, σt)
    eth = Petit.get_energy_threshold(particle_type; energy_threshold_ions=energy_threshold_ions)

    dfpars = Petit.DiffusionParams(ldrft_fixed, σt, σl, voxel_size, max_distance, eth, nbins, nsigma)

    # Load and validate input (from track_reco_mt.jl)
    hitsdf, ntot, nevents_config, last_event = Petit.load_and_validate_input(cmdir, input_file, ievt, levt)
    nevents_to_process = last_event - ievt + 1
    optimal_nthreads = Petit.get_optimal_threads(nthreads)

    # Determine if double-beta (bb0nu) or single-electron (xe137) based on filename
    is_double_beta = occursin("bb", lowercase(input_file))
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
    println("  Event type:         $(is_double_beta ? "double-beta (bb0nu)" : "single-electron (xe137)")")
    if use_random_ldrft
        println("  Drift length:       RANDOM [$lmin, $lmax] cm")
        println("  (Diffusion params computed per event)")
    else
        println("  Drift length:       $lmin cm (fixed)")
        println("  σ_t:                $(round(σt, digits=3)) mm")
        println("  σ_l:                $(round(σl, digits=3)) mm")
        println("  Voxel size:         $(round(voxel_size, digits=3)) mm")
        println("  Max distance:       $(round(max_distance, digits=3)) mm")
    end
    println("  Energy threshold:   $(round(eth, digits=3)) keV")
    println("  Blob radius:        $(round(Rb, digits=3)) mm")
    println("  Graph method:       $graph_method$(graph_method in ["kNN", "kNN_mutual"] ? " (k=$knn_k)" : "")")
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
                                         lmin=lmin, lmax=lmax,
                                         use_random_ldrft=use_random_ldrft,
                                         # Fixed params (used when lmin == lmax)
                                         ldrft_fixed=ldrft_fixed,
                                         σt_fixed=σt, σl_fixed=σl,
                                         voxel_size_fixed=voxel_size,
                                         max_distance_fixed=max_distance,
                                         mcvox_size_fixed=mcvox_size,
                                         eth_fixed=eth,
                                         dfpars_fixed=dfpars,
                                         # Common params
                                         nbins=nbins, nsigma=nsigma,
                                         Rb=Rb,
                                         # Params for recomputing diffusion
                                         particle_type=particle_type,
                                         dt=dt, dl=dl,
                                         tK=tK, edrift=edrift, Pbar=Pbar,
                                         energy_threshold_ions=energy_threshold_ions,
                                         is_double_beta=is_double_beta,
                                         print_level=print_level,
                                         # Graph method
                                         graph_method=graph_method,
                                         knn_k=knn_k,
                                         # Energy weighting options
                                         use_energy_weighting=use_energy_weighting,
                                         use_edge_energy_weighting=use_edge_energy_weighting,
                                         use_mst_fallback=use_mst_fallback)

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
            "lmin_cm" => lmin,
            "lmax_cm" => lmax,
            "random_drift" => use_random_ldrft,
            "energy_threshold_kev" => eth,
            "blob_radius_mm" => Rb,
            "nbins" => nbins,
            "nsigma" => nsigma,
            "graph_method" => graph_method,
            "knn_k" => knn_k
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
        ldrft_cm = Float64[],
        sigma_t_mm = Float64[],
        sigma_l_mm = Float64[],
        voxel_size_mm = Float64[],
        max_distance_mm = Float64[],
        # Track lengths
        track_length_mm = Float64[],
        track_length_mst_mm = Float64[],
        # MC blob analysis (default/KDT)
        Eb1_mc_keV = Float64[],
        Eb2_mc_keV = Float64[],
        asymmetry_mc = Float64[],
        # MC blob analysis (MST)
        Eb1_mc_mst_keV = Float64[],
        Eb2_mc_mst_keV = Float64[],
        asymmetry_mc_mst = Float64[],
        # RECO blob analysis (default/KDT)
        Eb1_keV = Float64[],
        Eb2_keV = Float64[],
        asymmetry = Float64[],
        # RECO blob analysis (MST)
        Eb1_mst_keV = Float64[],
        Eb2_mst_keV = Float64[],
        asymmetry_mst = Float64[],
        # MC extreme distances (default/KDT)
        d1_mc_mm = Float64[],
        d2_mc_mm = Float64[],
        # MC extreme distances (MST)
        d1_mc_mst_mm = Float64[],
        d2_mc_mst_mm = Float64[],
        # RECO extreme distances (default/KDT)
        d1_mm = Float64[],
        d2_mm = Float64[],
        # RECO extreme distances (MST)
        d1_mst_mm = Float64[],
        d2_mst_mm = Float64[]
    )

    for thread_results in all_itaca_results
        for r in thread_results
            push!(csv_df, (
                event = r.event_id,
                thread_id = r.thread_id,
                ldrft_cm = r.ldrft_cm,
                sigma_t_mm = r.sigma_t_mm,
                sigma_l_mm = r.sigma_l_mm,
                voxel_size_mm = r.voxel_size_mm,
                max_distance_mm = r.max_distance_mm,
                track_length_mm = r.track_length,
                track_length_mst_mm = r.track_length_mst,
                Eb1_mc_keV = r.Eb1_mc,
                Eb2_mc_keV = r.Eb2_mc,
                asymmetry_mc = r.asymmetry_mc,
                Eb1_mc_mst_keV = r.Eb1_mc_mst,
                Eb2_mc_mst_keV = r.Eb2_mc_mst,
                asymmetry_mc_mst = r.asymmetry_mc_mst,
                Eb1_keV = r.Eb1,
                Eb2_keV = r.Eb2,
                asymmetry = r.asymmetry,
                Eb1_mst_keV = r.Eb1_mst,
                Eb2_mst_keV = r.Eb2_mst,
                asymmetry_mst = r.asymmetry_mst,
                d1_mc_mm = r.d1_mc,
                d2_mc_mm = r.d2_mc,
                d1_mc_mst_mm = r.d1_mc_mst,
                d2_mc_mst_mm = r.d2_mc_mst,
                d1_mm = r.d1,
                d2_mm = r.d2,
                d1_mst_mm = r.d1_mst,
                d2_mst_mm = r.d2_mst
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
    push!(metadata_df, (parameter="lmin_cm", value=string(lmin)))
    push!(metadata_df, (parameter="lmax_cm", value=string(lmax)))
    push!(metadata_df, (parameter="random_drift", value=string(use_random_ldrft)))
    push!(metadata_df, (parameter="energy_threshold_keV", value=string(round(eth, digits=4))))
    push!(metadata_df, (parameter="blob_radius_mm", value=string(Rb)))
    push!(metadata_df, (parameter="nbins", value=string(nbins)))
    push!(metadata_df, (parameter="nsigma", value=string(nsigma)))
    push!(metadata_df, (parameter="graph_method", value=graph_method))
    push!(metadata_df, (parameter="knn_k", value=string(knn_k)))
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
        println("  --lmin=X              Min drift length (cm, default: 100.0)")
        println("  --lmax=X              Max drift length (cm, default: 100.0)")
        println("                        If lmin==lmax: fixed drift; lmin<lmax: random per event")
        println("  --tK=X                Temperature in K (default: 297.0)")
        println("  --edrift=X            Drift field in V/cm (default: 500.0)")
        println("  --Pbar=X              Pressure in bar (default: 15.0)")
        println("  --dt=X                Transverse diffusion coeff mm/√cm (default: 3.5)")
        println("  --dl=X                Longitudinal diffusion coeff mm/√cm (default: 0.9)")
        println("  --eth-ion=X           Energy threshold for ions (default: 10.0)")
        println("  --nbins=N             Diffusion histogram bins (default: 100)")
        println("  --nsigma=X            Diffusion extent in sigma (default: 3.0)")
        println("  --Rb=X                Blob radius in mm (default: 10.0)")
        println("  --print_level=LEVEL   quiet, verbose, very_verbose (default: quiet)")
        println("  --graph-method=M      Graph method: KDT, kNN, kNN_mutual (default: KDT)")
        println("  --knn-k=N             Number of neighbors for kNN methods (default: 10)")
        println("  --energy-weight=BOOL  Enable energy-weighted extreme finding (default: true)")
        println("  --edge-energy-weight=BOOL Enable edge-energy-weighted method (default: true)")
        println("  --mst-fallback=BOOL   Enable MST-based fallback method (default: false)")
        println("\nExamples:")
        println("  # Fixed drift at 100 cm:")
        println("  julia -t 8 itaca_track_reco_mt.jl /data/HD5t/itaca/ bb0nu.h5 --lmin=100 --lmax=100")
        println("  # Random drift in [50, 150] cm:")
        println("  julia -t 8 itaca_track_reco_mt.jl /data/HD5t/itaca/ bb0nu.h5 --lmin=50 --lmax=150")
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
    lmin = 100.0
    lmax = 100.0
    tK = 297.0
    edrift = 500.0
    Pbar = 15.0
    dt = 3.5
    dl = 0.9
    energy_threshold_ions = 10.0
    nbins = 100
    nsigma = 3.0
    Rb = 10.0
    print_level = "quiet"
    graph_method = "KDT"
    knn_k = 10
    use_energy_weighting = true
    use_edge_energy_weighting = true
    use_mst_fallback = false

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
                elseif key == "lmin"
                    lmin = parse(Float64, value)
                elseif key == "lmax"
                    lmax = parse(Float64, value)
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
                elseif key == "Rb"
                    Rb = parse(Float64, value)
                elseif key == "print_level"
                    print_level = String(value)
                elseif key == "graph-method"
                    graph_method = String(value)
                elseif key == "knn-k"
                    knn_k = parse(Int, value)
                elseif key == "energy-weight"
                    use_energy_weighting = parse(Bool, lowercase(value))
                elseif key == "edge-energy-weight"
                    use_edge_energy_weighting = parse(Bool, lowercase(value))
                elseif key == "mst-fallback"
                    use_mst_fallback = parse(Bool, lowercase(value))
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
                        lmin=lmin,
                        lmax=lmax,
                        tK=tK,
                        edrift=edrift,
                        Pbar=Pbar,
                        dt=dt,
                        dl=dl,
                        energy_threshold_ions=energy_threshold_ions,
                        nbins=nbins,
                        nsigma=nsigma,
                        Rb=Rb,
                        print_level=print_level,
                        graph_method=graph_method,
                        knn_k=knn_k,
                        use_energy_weighting=use_energy_weighting,
                        use_edge_energy_weighting=use_edge_energy_weighting,
                        use_mst_fallback=use_mst_fallback)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
