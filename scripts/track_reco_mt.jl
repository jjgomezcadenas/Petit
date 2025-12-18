#!/usr/bin/env julia

"""
Multi-threaded track reconstruction script.

This script processes events from ievt to levt using multiple threads,
reconstructing tracks with central path computation and saving both
the track and central path to HDF5 files.

The output can then be used to run blob analysis with different radii
without having to redo the reconstruction.
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using HDF5
using DataFrames
using Statistics
using Graphs

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

# Constants
const f = 1e+5/2.5  # ions per MeV
const fkeV = f*1e-3  # ions per keV

"""
    transform_hits_df(df; energy_to_electrons)

Transform hits DataFrame to include electron counts.
"""
function transform_hits_df(df::DataFrame; energy_to_electrons::Float64=1e5/2.5)
    df2 = select(df, Not([:time, :label, :particle_id, :hit_id]))
    df2.electrons = round.(Int, df2.energy .* energy_to_electrons)
    return df2
end

"""
    RecoResult

Structure to hold reconstruction results for a single event.
Includes raw path, projected voxels, and KDE energy profiles for both RECO and MC.
"""
struct RecoResult
    event_id::Int
    track::Petit.Tracks
    path::DataFrame              # Raw path (x, y, z, s) - NOT smoothed
    track_length::Float64
    confidence::Float64
    mc_path::DataFrame           # MC path (x, y, z, energy, s, primary_electron)
    reco_s::Vector{Float64}      # Arc-length of each reco voxel on path
    kde_s::Vector{Float64}       # RECO KDE evaluation points (0 to track_length)
    reco_kde_f::Vector{Float64}  # RECO energy density f(s)
    mc_kde_s::Vector{Float64}    # MC KDE evaluation points (0 to mc_track_length)
    mc_kde_f::Vector{Float64}    # MC energy density f(s)
    kde_bandwidth::Float64       # RECO KDE bandwidth
    mc_kde_bandwidth::Float64    # MC KDE bandwidth
    d1::Float64                  # Distance: reco extreme 1 to matched MC extreme
    d2::Float64                  # Distance: reco extreme 2 to matched MC extreme
end

"""
    process_single_event(hitsdf, event_id, σt, voxel_size, max_distance,
                         energy_threshold_keV, dfpars, kde_bandwidth, mc_kde_bandwidth,
                         n_kde_eval, nbins, nsigma, mcvox_size)

Process a single event: diffuse, voxelize, make tracks, project onto path, compute KDE.
Computes KDE energy profiles for both RECO voxels and MC path.
Returns RecoResult if successful, nothing otherwise.
"""
function process_single_event(hitsdf::DataFrame, event_id::Int,
                              σt::Float64, σl::Float64,
                              voxel_size::Float64, max_distance::Float64,
                              energy_threshold_keV::Float64,
                              dfpars::Petit.DiffusionParams,
                              kde_bandwidth::Float64,
                              mc_kde_bandwidth::Float64,
                              n_kde_eval::Int,
                              nbins::Int, nsigma::Float64,
                              mcvox_size::Float64;
                              is_double_beta::Bool=false)

    # Get event hits
    event_df = Petit.get_event(hitsdf, event_id)
    if nrow(event_df) == 0
        return nothing
    end

    # Compute MC path (voxelized primary particle trajectory with arc-length)
    # First/last rows give MC truth extremes
    mc_path = Petit.compute_mc_path(event_df, mcvox_size; is_double_beta=is_double_beta)
    if nrow(mc_path) == 0
        return nothing  # No primary particle hits
    end

    # Transform hits (removes time, particle_id, etc.)
    event_mc = transform_hits_df(event_df)

    # Diffuse and voxelize
    diffused_df = Petit.diffuse_xyz_image_mc(event_mc;
                                              sigma_t_mm=σt,
                                              sigma_l_mm=σl,
                                              nbins=nbins,
                                              nsigma=nsigma)

    voxels = Petit.voxelize_event(diffused_df, voxel_size)

    # Make tracks
    tracks = Petit.make_tracks(voxels;
                               max_distance_mm=max_distance,
                               energy_threshold_kev=energy_threshold_keV,
                               diffusion=dfpars)

    # Only keep single-track events
    if length(tracks) != 1
        return nothing
    end

    track = tracks[1]

    # Walk track to get extremes and path
    walk_result = Petit.walk_track_from_extremes(track)

    if isnothing(walk_result.extremes[1]) || isempty(walk_result.path_indices)
        return nothing
    end

    # Get raw path (NO smoothing)
    path = Petit.get_raw_path(track, walk_result.path_indices)

    if nrow(path) < 2
        return nothing
    end

    track_length = path.s[end]

    # Project RECO voxels onto path
    reco_s = Petit.project_voxels_to_path(track.voxels, path)
    reco_E = Vector{Float64}(track.voxels.energy)

    # MC path data
    mc_s = Vector{Float64}(mc_path.s)
    mc_E = Vector{Float64}(mc_path.energy)
    mc_track_length = mc_path.s[end]

    # Create RECO evaluation grid (based on RECO track length only)
    kde_s = collect(range(0.0, track_length, length=n_kde_eval))

    # Create MC evaluation grid (based on MC track length)
    mc_kde_s = collect(range(0.0, mc_track_length, length=n_kde_eval))

    # Compute RECO KDE on RECO grid
    reco_kde_f, _ = Petit.energy_weighted_kde(reco_s, reco_E, kde_s; bandwidth=kde_bandwidth)

    # Compute MC KDE on MC grid
    mc_kde_f, _ = Petit.energy_weighted_kde(mc_s, mc_E, mc_kde_s; bandwidth=mc_kde_bandwidth)

    # Compute distances between reco and MC extremes
    extreme_dists = Petit.compute_extreme_distances(path, mc_path)

    return RecoResult(event_id, track, path,
                      track_length, walk_result.confidence,
                      mc_path, reco_s, kde_s, reco_kde_f,
                      mc_kde_s, mc_kde_f,
                      kde_bandwidth, mc_kde_bandwidth,
                      extreme_dists.d1, extreme_dists.d2)
end

"""
    analysis_loop_reco_mt(hitsdf, thread_id; kwargs...)

Multi-threaded analysis loop for track reconstruction with KDE.
Returns vector of RecoResult.
"""
function analysis_loop_reco_mt(hitsdf::DataFrame, thread_id::Int;
                               events_to_run=nothing,
                               initial_event::Int=1,
                               σt::Float64=3.0,
                               σl::Float64=0.0,
                               voxel_size::Float64=3.0,
                               max_distance::Float64=4.5,
                               energy_threshold_keV::Float64=0.25,
                               kde_bandwidth::Float64=5.0,
                               mc_kde_bandwidth::Float64=1.0,
                               n_kde_eval::Int=200,
                               nbins::Int=100,
                               nsigma::Float64=3.0,
                               mcvox_size::Float64=1.0,
                               dfpars::Petit.DiffusionParams=Petit.DiffusionParams(),
                               is_double_beta::Bool=false)

    event_ids = Petit.get_events_to_process(hitsdf, events_to_run, initial_event)
    n_single_track = 0
    n_events_processed = 0
    n_failed = 0
    results = RecoResult[]

    for event_id in event_ids
        n_events_processed += 1
        try
            result = process_single_event(hitsdf, event_id, σt, σl,
                                         voxel_size, max_distance,
                                         energy_threshold_keV, dfpars,
                                         kde_bandwidth, mc_kde_bandwidth,
                                         n_kde_eval, nbins, nsigma, mcvox_size;
                                         is_double_beta=is_double_beta)
            if !isnothing(result)
                push!(results, result)
                n_single_track += 1
            end
        catch e
            println("Thread $thread_id: Warning: Error processing event $event_id: $e")
            n_failed += 1
        end
    end

    println("Thread $thread_id: Completed! Processed $n_events_processed events, $n_single_track single-track.")
    return results
end

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

    # Save as matrix with column names
    path_data = Matrix(path)
    group["path_data"] = path_data
    group["path_columns"] = String.(names(path))
end

"""
    save_mc_path_to_hdf5(mc_path, group)

Save a MC path DataFrame to an HDF5 group.
MC path contains: x, y, z (energy-weighted position), energy (sum), s (arc-length)
"""
function save_mc_path_to_hdf5(mc_path::DataFrame, group)
    if nrow(mc_path) == 0
        group["mc_path_data"] = zeros(Float64, 0, 0)
        group["mc_path_columns"] = String[]
        return
    end

    # Save as matrix with column names
    mc_data = Matrix(mc_path)
    group["mc_path_data"] = mc_data
    group["mc_path_columns"] = String.(names(mc_path))
end

"""
    save_reco_results_to_hdf5(results, output_path, metadata)

Save reconstruction results (tracks + paths + KDE) to HDF5 file.
"""
function save_reco_results_to_hdf5(results::Vector{RecoResult}, output_path::String, metadata::Dict)
    h5open(output_path, "w") do fid
        # Save metadata as attributes
        for (key, val) in metadata
            attrs(fid)[key] = val
        end
        attrs(fid)["total_tracks_saved"] = length(results)

        if isempty(results)
            return
        end

        println("    Saving $(length(results)) reconstructed tracks to HDF5...")
        flush(stdout)

        for (idx, result) in enumerate(results)
            if idx % 10 == 0 || idx == length(results)
                print("\r    Saving track $idx/$(length(results)) to HDF5...")
                flush(stdout)
            end

            track_group_name = "batch_1/track_$(idx)"
            g = create_group(fid, track_group_name)

            track = result.track

            # Save event metadata
            g["event_id"] = result.event_id
            g["track_length"] = result.track_length
            g["confidence"] = result.confidence
            g["kde_bandwidth"] = result.kde_bandwidth
            g["mc_kde_bandwidth"] = result.mc_kde_bandwidth

            # Save voxels data
            voxels_data = Matrix(track.voxels)
            g["voxels"] = voxels_data
            g["voxel_columns"] = String.(names(track.voxels))

            # Save graph edges
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

            # Save components
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

            # Save raw path (unsmoothed)
            save_path_to_hdf5(result.path, g)

            # Save MC path
            save_mc_path_to_hdf5(result.mc_path, g)

            # Save projected voxel arc-lengths
            g["reco_s"] = result.reco_s

            # Save KDE results
            g["kde_s"] = result.kde_s
            g["reco_kde_f"] = result.reco_kde_f
            g["mc_kde_s"] = result.mc_kde_s
            g["mc_kde_f"] = result.mc_kde_f

            # Save extreme distances
            g["d1"] = result.d1
            g["d2"] = result.d2
        end
        println()  # New line after progress
    end
end

"""
    event_loop_reco_mt(cmdir, output_base; kwargs...)

Multi-threaded event processing for track reconstruction with KDE.

Produces HDF5 files containing:
- Reconstructed tracks (voxels, graph)
- Raw paths (unsmoothed track skeleton)
- KDE energy profiles (for RECO and MC)
- Track metadata (length, confidence, event_id)

# Diffusion parameters:
- σt: transverse diffusion (mm). If negative, computed from ion diffusion formula.
- σl: longitudinal diffusion (mm). Default 0.0.
"""
function event_loop_reco_mt(cmdir::String, output_base::String;
                            input_file::String="0nubb.next.h5",
                            ievt::Int=1,
                            levt::Int=-1,
                            nthreads::Int=1,
                            ldrft::Float64=100.0,
                            tK::Float64=297.0,
                            edrift::Float64=500.0,
                            σt::Float64=-1.0,
                            σl::Float64=0.0,
                            voxel_scale::Float64=1.0,
                            voxel_distance_scale::Float64=1.5,
                            energy_threshold_ions::Float64=10.0,
                            nbins::Int=100,
                            nsigma::Float64=3.0,
                            kde_bandwidth::Float64=5.0,
                            mc_kde_bandwidth::Float64=1.0,
                            n_kde_eval::Int=200,
                            mcvox_size::Float64=1.0)

    # Compute diffusion parameters
    # If σt < 0, compute from ion diffusion formula
    σt_use = σt < 0 ? Petit.sigma_t_ion_mm(tK, ldrft, edrift) : σt
    σl_use = σl

    voxel_size = σt_use * voxel_scale
    max_distance = voxel_size * voxel_distance_scale
    energy_threshold_keV = energy_threshold_ions / fkeV

    dfpars = Petit.DiffusionParams(ldrft, σt_use, σl_use,
                                   voxel_size, max_distance, energy_threshold_keV,
                                   nbins, nsigma)

    # Load and validate input
    hitsdf, ntot, nevents_config, last_event = Petit.load_and_validate_input(cmdir, input_file, ievt, levt)
    nevents_to_process = last_event - ievt + 1
    optimal_nthreads = Petit.get_optimal_threads(nthreads)
    thread_ranges = Petit.split_events_for_threads(nevents_to_process, optimal_nthreads, ievt)

    # Print configuration
    println("\n" * "="^60)
    println("TRACK RECONSTRUCTION + KDE (MT)")
    println("="^60)
    println("Configuration:")
    println("  Input file:         $input_file")
    println("  Output base:        $output_base")
    println("  Events:             $ievt to $last_event ($nevents_to_process total)")
    println("  Threads:            $optimal_nthreads")
    println("  Drift length:       $ldrft cm")
    println("  σ_t:                $(round(σt_use, digits=3)) mm")
    println("  σ_l:                $(round(σl_use, digits=3)) mm")
    println("  Voxel size:         $(round(voxel_size, digits=3)) mm")
    println("  Max distance:       $(round(max_distance, digits=3)) mm")
    println("  KDE bandwidth:      $(round(kde_bandwidth, digits=3)) mm (RECO)")
    println("  MC KDE bandwidth:   $(round(mc_kde_bandwidth, digits=3)) mm")
    println("  KDE eval points:    $n_kde_eval")
    println("  Energy threshold:   $(round(energy_threshold_keV, digits=3)) keV")
    println("  MC voxel size:      $(round(mcvox_size, digits=3)) mm")
    println("="^60)

    # Process in parallel
    println("\nStarting multi-threaded processing...")
    all_results = Vector{Any}(undef, length(thread_ranges))

    Threads.@threads for i in 1:length(thread_ranges)
        start_evt, num_evts = thread_ranges[i]

        results = analysis_loop_reco_mt(hitsdf, i;
                                        events_to_run=num_evts,
                                        initial_event=start_evt,
                                        σt=σt_use,
                                        σl=σl_use,
                                        voxel_size=voxel_size,
                                        max_distance=max_distance,
                                        energy_threshold_keV=energy_threshold_keV,
                                        kde_bandwidth=kde_bandwidth,
                                        mc_kde_bandwidth=mc_kde_bandwidth,
                                        n_kde_eval=n_kde_eval,
                                        nbins=nbins,
                                        nsigma=nsigma,
                                        mcvox_size=mcvox_size,
                                        dfpars=dfpars)

        output_file = "$(output_base)_th_$i.h5"
        output_path = joinpath(cmdir, output_file)

        metadata = Dict(
            "input_file" => input_file,
            "nevents_from_config" => nevents_config,
            "total_events_in_file" => ntot,
            "thread_id" => i,
            "first_event_processed" => start_evt,
            "last_event_processed" => start_evt + num_evts - 1,
            "events_processed" => num_evts,
            "ldrift_cm" => ldrft,
            "sigma_t_mm" => σt_use,
            "sigma_l_mm" => σl_use,
            "voxel_size_mm" => voxel_size,
            "max_distance_mm" => max_distance,
            "kde_bandwidth_mm" => kde_bandwidth,
            "mc_kde_bandwidth_mm" => mc_kde_bandwidth,
            "n_kde_eval" => n_kde_eval,
            "energy_threshold_kev" => energy_threshold_keV,
            "mcvox_size_mm" => mcvox_size,
            "nbins" => nbins,
            "nsigma" => nsigma
        )

        save_reco_results_to_hdf5(results, output_path, metadata)

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

    # Print summary
    total_tracks = sum(r["tracks_saved"] for r in all_results)
    println("\n" * "="^60)
    println("PROCESSING COMPLETE")
    println("="^60)
    println("Total tracks saved: $total_tracks")
    println("\nOutput files:")
    for r in all_results
        println("  $(r["output_file"]) ($(r["tracks_saved"]) tracks)")
    end

    return Dict(
        "thread_results" => all_results,
        "total_tracks" => total_tracks,
        "events_processed" => nevents_to_process
    )
end

"""
    main()

Main function for command-line usage.
"""
function main()
    # Default values
    ievt = 1
    levt = -1
    nthreads = 1
    ldrft = 100.0
    σt = -1.0  # Negative means compute from ion diffusion formula
    σl = 0.0
    voxel_scale = 1.0
    voxel_distance_scale = 1.5
    kde_bandwidth = 5.0     # RECO KDE bandwidth (mm)
    mc_kde_bandwidth = 1.0  # MC KDE bandwidth (mm)
    n_kde_eval = 200        # Number of KDE evaluation points
    mcvox_size = 1.0        # MC path voxel size (mm)

    # Check minimum required arguments
    if length(ARGS) < 3
        println("Error: Missing required arguments")
        println("\nUsage: julia -t <nthreads> track_reco_mt.jl <cmdir> <input_file> <output_base> [options]")
        println("\nRequired arguments:")
        println("  cmdir           Directory containing the input file")
        println("  input_file      Name of the HDF5 input file")
        println("  output_base     Base name for output files (no .h5 extension)")
        println("\nOptional arguments:")
        println("  --ievt=N                   First event to process (default: 1)")
        println("  --levt=N                   Last event to process (default: -1, all)")
        println("  --nthreads=N               Number of threads to use (default: 1)")
        println("  --ldrft=X                  Drift length in cm (default: 100.0)")
        println("  --sigmat=X                 Transverse diffusion in mm (default: -1, compute from ion)")
        println("  --sigmal=X                 Longitudinal diffusion in mm (default: 0.0)")
        println("  --voxel-scale=X            Voxel scale factor (default: 1.0)")
        println("  --voxel-distance-scale=X   Voxel distance scale (default: 1.5)")
        println("  --kde-bandwidth=X          RECO KDE bandwidth in mm (default: 5.0)")
        println("  --mc-kde-bandwidth=X       MC KDE bandwidth in mm (default: 1.0)")
        println("  --n-kde-eval=N             Number of KDE evaluation points (default: 200)")
        println("  --mcvox=X                  MC path voxel size in mm (default: 1.0)")
        println("\nExample:")
        println("  julia -t 8 track_reco_mt.jl /path/to/data/ input.h5 output --nthreads=8 --ldrft=100")
        println("  julia -t 8 track_reco_mt.jl /path/to/data/ input.h5 output --kde-bandwidth=3.0 --mcvox=1.0")
        exit(1)
    end

    # Parse required arguments
    cmdir = ARGS[1]
    input_file = ARGS[2]
    output_base = ARGS[3]

    # Parse optional arguments
    for arg in ARGS[4:end]
        if startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) != 2
                println("Warning: Ignoring malformed argument: $arg")
                continue
            end
            key, value = parts

            try
                if key == "ievt"
                    ievt = parse(Int, value)
                elseif key == "levt"
                    levt = parse(Int, value)
                elseif key == "nthreads"
                    nthreads = parse(Int, value)
                elseif key == "ldrft"
                    ldrft = parse(Float64, value)
                elseif key == "sigmat"
                    σt = parse(Float64, value)
                elseif key == "sigmal"
                    σl = parse(Float64, value)
                elseif key == "voxel-scale"
                    voxel_scale = parse(Float64, value)
                elseif key == "voxel-distance-scale"
                    voxel_distance_scale = parse(Float64, value)
                elseif key == "kde-bandwidth"
                    kde_bandwidth = parse(Float64, value)
                elseif key == "mc-kde-bandwidth"
                    mc_kde_bandwidth = parse(Float64, value)
                elseif key == "n-kde-eval"
                    n_kde_eval = parse(Int, value)
                elseif key == "mcvox"
                    mcvox_size = parse(Float64, value)
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
    event_loop_reco_mt(cmdir, output_base;
                       input_file=input_file,
                       ievt=ievt,
                       levt=levt,
                       nthreads=nthreads,
                       ldrft=ldrft,
                       σt=σt,
                       σl=σl,
                       voxel_scale=voxel_scale,
                       voxel_distance_scale=voxel_distance_scale,
                       kde_bandwidth=kde_bandwidth,
                       mc_kde_bandwidth=mc_kde_bandwidth,
                       n_kde_eval=n_kde_eval,
                       mcvox_size=mcvox_size)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
