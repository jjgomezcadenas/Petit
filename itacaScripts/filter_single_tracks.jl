#!/usr/bin/env julia

"""
Filter events by energy and single-track reconstruction.

Combines energy filtering with track reconstruction to select only events
that reconstruct as a single track (no diffusion applied - uses raw MC hits).

Features:
- Filters events by total energy in [emin, emax] keV
- Reconstructs each event with kNN track building
- Selects only events with exactly 1 reconstructed track
- Multi-threaded writing (one file per thread)
- Preserves HDF5 structure identical to input

Usage:
    julia -t <nthreads> filter_single_tracks.jl <cmdir> <input_file> [options]

Options:
    --emin=X            Minimum energy in keV (default: 2400)
    --emax=X            Maximum energy in keV (default: 2500)
    --voxel-size=X      Voxel size in mm (default: 2.0)
    --max-distance=X    Max distance for track building in mm (default: 2.0)
    --nthreads=N        Number of threads (default: 1)

Output:
    Creates directory: <cmdir>/<basename>_filtered_ecut_<emin>_<emax>_st/
    Contains:
    - <basename>_part_1.h5, <basename>_part_2.h5, ... (single-track events)
    - filter_summary.txt (statistics and parameters)
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using HDF5
using DataFrames
using Statistics
using Dates

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

# =============================================================================
# Helper Functions
# =============================================================================

"""
    compute_event_energies(hitsdf::DataFrame)

Compute total energy (in keV) for each event.
Returns a DataFrame with columns: event_id, total_energy_keV
"""
function compute_event_energies(hitsdf::DataFrame)
    event_energies = combine(groupby(hitsdf, :event_id)) do df
        (total_energy_keV = 1e3 * sum(df.energy),)  # Convert MeV to keV
    end
    return event_energies
end

"""
    filter_events_by_energy(event_energies::DataFrame, emin::Float64, emax::Float64)

Filter events by total energy in [emin, emax] keV.
Returns vector of event IDs that pass the cut.
"""
function filter_events_by_energy(event_energies::DataFrame, emin::Float64, emax::Float64)
    mask = (event_energies.total_energy_keV .>= emin) .& (event_energies.total_energy_keV .<= emax)
    return event_energies.event_id[mask]
end

"""
    is_single_track(hitsdf::DataFrame, event_id::Int, voxel_size::Float64, max_distance::Float64)

Check if an event reconstructs as exactly one track.
Returns true if n_tracks == 1, false otherwise.
"""
function is_single_track(hitsdf::DataFrame, event_id::Int, voxel_size::Float64, max_distance::Float64)
    # Get event hits
    event_df = Petit.get_event(hitsdf, event_id)

    if nrow(event_df) == 0
        return false
    end

    # Voxelize raw MC hits (no diffusion)
    voxels = Petit.voxelize_event(event_df, voxel_size)

    if nrow(voxels) == 0
        return false
    end

    # Build tracks using kNN
    tracks = Petit.make_tracks(voxels;
                               method="kNN",
                               k=10,
                               max_distance_mm=max_distance,
                               energy_threshold_kev=1.0)

    return length(tracks) == 1
end

"""
    create_output_dir(cmdir::String, base_name::String, emin::Float64, emax::Float64)

Create output directory for filtered single-track data.
Returns the path to the created directory.
"""
function create_output_dir(cmdir::String, base_name::String, emin::Float64, emax::Float64)
    # Format energy values (integer if whole number)
    emin_str = isinteger(emin) ? string(Int(emin)) : string(emin)
    emax_str = isinteger(emax) ? string(Int(emax)) : string(emax)

    output_dir = joinpath(cmdir, "$(base_name)_filtered_ecut_$(emin_str)_$(emax_str)_st")

    if !isdir(output_dir)
        mkpath(output_dir)
        println("Created output directory: $output_dir")
    else
        println("Output directory exists: $output_dir")
    end

    return output_dir
end

"""
    write_filtered_hits_hdf5(output_path::String, hits_df::DataFrame, metadata::Dict)

Write filtered hits to HDF5 file, preserving structure similar to input.
"""
function write_filtered_hits_hdf5(output_path::String, hits_df::DataFrame, metadata::Dict)
    h5open(output_path, "w") do fid
        # Write metadata as attributes
        for (key, val) in metadata
            attrs(fid)[key] = string(val)
        end

        # Create hits group and write data
        hits_group = create_group(fid, "MC/hits")

        if nrow(hits_df) > 0
            # Write each column
            for col in names(hits_df)
                data = hits_df[!, col]
                if eltype(data) <: AbstractString
                    hits_group[col] = String.(data)
                else
                    hits_group[col] = Vector(data)
                end
            end
        end

        # Store column names for reference
        attrs(hits_group)["columns"] = String.(names(hits_df))
        attrs(hits_group)["n_hits"] = nrow(hits_df)
        attrs(hits_group)["n_events"] = length(unique(hits_df.event_id))
    end
end

"""
    write_summary_file(output_dir, input_file, params, stats)

Write summary file with filter statistics.
"""
function write_summary_file(output_dir::String, input_file::String,
                            emin::Float64, emax::Float64,
                            voxel_size::Float64, max_distance::Float64,
                            n_total_events::Int, n_energy_pass::Int, n_single_track::Int,
                            energy_stats::NamedTuple, thread_results::Vector)

    summary_path = joinpath(output_dir, "filter_summary.txt")

    open(summary_path, "w") do f
        println(f, "=" ^ 70)
        println(f, "SINGLE-TRACK FILTER SUMMARY")
        println(f, "=" ^ 70)
        println(f)
        println(f, "Timestamp: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(f)
        println(f, "Input/Output:")
        println(f, "  Input file:       $input_file")
        println(f, "  Output directory: $output_dir")
        println(f)
        println(f, "Filter Parameters:")
        println(f, "  Energy min:       $emin keV")
        println(f, "  Energy max:       $emax keV")
        println(f, "  Voxel size:       $voxel_size mm")
        println(f, "  Max distance:     $max_distance mm")
        println(f, "  Track method:     kNN (k=10)")
        println(f, "  Diffusion:        DISABLED (raw MC hits)")
        println(f)
        println(f, "Event Statistics:")
        println(f, "  Total events read:        $n_total_events")
        println(f, "  Events passing energy:    $n_energy_pass ($(round(100 * n_energy_pass / n_total_events, digits=2))%)")
        println(f, "  Single-track events:      $n_single_track ($(round(100 * n_single_track / n_energy_pass, digits=2))% of energy-pass)")
        println(f, "  Overall efficiency:       $(round(100 * n_single_track / n_total_events, digits=4))%")
        println(f)
        println(f, "Energy Distribution (all events):")
        println(f, "  Mean:   $(round(energy_stats.mean, digits=2)) keV")
        println(f, "  Std:    $(round(energy_stats.std, digits=2)) keV")
        println(f, "  Min:    $(round(energy_stats.min, digits=2)) keV")
        println(f, "  Max:    $(round(energy_stats.max, digits=2)) keV")
        println(f)
        println(f, "Output Files:")
        for r in thread_results
            if r["n_events"] > 0
                println(f, "  $(r["output_file"]): $(r["n_events"]) events, $(r["n_hits"]) hits")
            end
        end
        println(f)
        println(f, "=" ^ 70)
    end

    println("Summary saved to: $summary_path")
    return summary_path
end

# =============================================================================
# Main Processing Functions
# =============================================================================

"""
    filter_and_write_thread(hitsdf, event_ids, thread_id, output_dir, base_name, metadata)

Process events for a single thread: filter hits and write to HDF5.
"""
function filter_and_write_thread(hitsdf::DataFrame, event_ids::Vector{Int},
                                  thread_id::Int, output_dir::String, base_name::String, metadata::Dict)
    # Filter hits for assigned events
    hits_filtered = filter(row -> row.event_id in event_ids, hitsdf)

    # Write to HDF5
    output_file = "$(base_name)_part_$thread_id.h5"
    output_path = joinpath(output_dir, output_file)

    thread_metadata = copy(metadata)
    thread_metadata["thread_id"] = thread_id
    thread_metadata["n_events_in_part"] = length(event_ids)
    thread_metadata["n_hits_in_part"] = nrow(hits_filtered)

    write_filtered_hits_hdf5(output_path, hits_filtered, thread_metadata)

    return Dict(
        "thread_id" => thread_id,
        "output_file" => output_file,
        "n_events" => length(event_ids),
        "n_hits" => nrow(hits_filtered)
    )
end

"""
    filter_single_tracks(cmdir, input_file; kwargs...)

Main function: filter events by energy and single-track reconstruction.
"""
function filter_single_tracks(cmdir::String, input_file::String;
                               emin::Float64=2400.0, emax::Float64=2500.0,
                               voxel_size::Float64=2.0, max_distance::Float64=2.0,
                               nthreads::Int=1)

    println("\n" * "=" ^ 70)
    println("SINGLE-TRACK FILTER")
    println("=" ^ 70)

    # Load data
    input_path = joinpath(cmdir, input_file)
    println("Loading data from: $input_path")

    if !isfile(input_path)
        error("File not found: $input_path")
    end

    dfs = Petit.get_dataset_dfs(input_path)
    hitsdf = dfs["hits"]

    # Compute event energies
    println("Computing event energies...")
    event_energies = compute_event_energies(hitsdf)
    n_total_events = nrow(event_energies)
    println("Total events: $n_total_events")

    # Energy statistics (for summary)
    energy_stats = (
        mean = mean(event_energies.total_energy_keV),
        std = std(event_energies.total_energy_keV),
        min = minimum(event_energies.total_energy_keV),
        max = maximum(event_energies.total_energy_keV)
    )

    # Filter events by energy
    println("\nFiltering events with energy in [$emin, $emax] keV...")
    energy_pass_ids = filter_events_by_energy(event_energies, emin, emax)
    n_energy_pass = length(energy_pass_ids)

    println("Events passing energy cut: $n_energy_pass / $n_total_events ($(round(100 * n_energy_pass / n_total_events, digits=2))%)")

    if n_energy_pass == 0
        println("\nWarning: No events pass the energy cut!")
        println("Energy range in data: [$(round(energy_stats.min, digits=2)), $(round(energy_stats.max, digits=2))] keV")
        return nothing
    end

    # Filter for single-track events
    println("\nChecking track reconstruction (voxel=$(voxel_size)mm, maxd=$(max_distance)mm)...")
    single_track_ids = Int[]

    for (i, event_id) in enumerate(energy_pass_ids)
        if i % 100 == 0 || i == n_energy_pass
            print("\r  Processing event $i / $n_energy_pass...")
        end

        if is_single_track(hitsdf, event_id, voxel_size, max_distance)
            push!(single_track_ids, event_id)
        end
    end
    println()

    n_single_track = length(single_track_ids)
    println("Single-track events: $n_single_track / $n_energy_pass ($(round(100 * n_single_track / n_energy_pass, digits=2))%)")

    if n_single_track == 0
        println("\nWarning: No single-track events found!")
        return nothing
    end

    # Get base name from input file (without extension)
    base_name = splitext(input_file)[1]

    # Create output directory
    output_dir = create_output_dir(cmdir, base_name, emin, emax)

    # Determine number of threads
    optimal_nthreads = Petit.get_optimal_threads(nthreads)
    if optimal_nthreads > n_single_track
        optimal_nthreads = n_single_track
        println("Reducing threads to $optimal_nthreads (one per event)")
    end

    # Split events among threads
    thread_ranges = Petit.split_events_for_threads(n_single_track, optimal_nthreads, 1)

    println("\nConfiguration:")
    println("  Energy cut:        [$emin, $emax] keV")
    println("  Voxel size:        $voxel_size mm")
    println("  Max distance:      $max_distance mm")
    println("  Events to write:   $n_single_track")
    println("  Threads:           $optimal_nthreads")
    println("  Output:            $output_dir")

    # Prepare metadata
    metadata = Dict(
        "source_file" => input_file,
        "emin_keV" => emin,
        "emax_keV" => emax,
        "voxel_size_mm" => voxel_size,
        "max_distance_mm" => max_distance,
        "total_events_in_source" => n_total_events,
        "events_passing_energy" => n_energy_pass,
        "single_track_events" => n_single_track,
        "filter_type" => "single_track",
        "timestamp" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    )

    # Process in parallel
    println("\nWriting filtered data...")
    results = Vector{Dict}(undef, length(thread_ranges))

    Threads.@threads for i in 1:length(thread_ranges)
        start_idx, num_events = thread_ranges[i]
        end_idx = start_idx + num_events - 1

        # Get event IDs for this thread
        thread_event_ids = single_track_ids[start_idx:end_idx]

        results[i] = filter_and_write_thread(hitsdf, thread_event_ids, i, output_dir, base_name, metadata)

        println("  Thread $i: wrote $(results[i]["n_events"]) events to $(results[i]["output_file"])")
    end

    # Write summary file
    summary_path = write_summary_file(output_dir, input_file,
                                       emin, emax, voxel_size, max_distance,
                                       n_total_events, n_energy_pass, n_single_track,
                                       energy_stats, results)

    # Print final summary
    println("\n" * "=" ^ 70)
    println("FILTER COMPLETE")
    println("=" ^ 70)
    println("Total events read:       $n_total_events")
    println("Events passing energy:   $n_energy_pass ($(round(100 * n_energy_pass / n_total_events, digits=2))%)")
    println("Single-track events:     $n_single_track ($(round(100 * n_single_track / n_energy_pass, digits=2))%)")
    println("Overall efficiency:      $(round(100 * n_single_track / n_total_events, digits=4))%")
    println("Output directory:        $output_dir")
    println("=" ^ 70)

    return Dict(
        "output_dir" => output_dir,
        "n_total_events" => n_total_events,
        "n_energy_pass" => n_energy_pass,
        "n_single_track" => n_single_track,
        "thread_results" => results,
        "summary_path" => summary_path
    )
end

# =============================================================================
# Command Line Interface
# =============================================================================

function main()
    # Check minimum required arguments
    if length(ARGS) < 2
        println("Error: Missing required arguments")
        println()
        println("Usage: julia -t <nthreads> filter_single_tracks.jl <cmdir> <input_file> [options]")
        println()
        println("Required arguments:")
        println("  cmdir           Directory containing the input file")
        println("  input_file      Name of the HDF5 input file")
        println()
        println("Optional arguments:")
        println("  --emin=X          Minimum energy in keV (default: 2400)")
        println("  --emax=X          Maximum energy in keV (default: 2500)")
        println("  --voxel-size=X    Voxel size in mm (default: 2.0)")
        println("  --max-distance=X  Max distance for track building in mm (default: 2.0)")
        println("  --nthreads=N      Number of threads (default: 1)")
        println()
        println("Output:")
        println("  Creates directory: <cmdir>/<basename>_filtered_ecut_<emin>_<emax>_st/")
        println("  Contains: <basename>_part_1.h5, ..., filter_summary.txt")
        println()
        println("Examples:")
        println("  # Default parameters (energy 2400-2500, voxel 2mm, maxd 2mm):")
        println("  julia -t 4 filter_single_tracks.jl /data/itaca/ bb0nu.h5")
        println()
        println("  # Custom parameters:")
        println("  julia -t 4 filter_single_tracks.jl /data/itaca/ bb0nu.h5 --voxel-size=1.5 --max-distance=2.0")
        exit(1)
    end

    # Parse required arguments
    cmdir = ARGS[1]
    input_file = ARGS[2]

    # Default values
    emin = 2400.0
    emax = 2500.0
    voxel_size = 2.0
    max_distance = 2.0
    nthreads = Threads.nthreads()  # Use Julia's thread count from -t flag

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
                if key == "emin"
                    emin = parse(Float64, value)
                elseif key == "emax"
                    emax = parse(Float64, value)
                elseif key == "voxel-size"
                    voxel_size = parse(Float64, value)
                elseif key == "max-distance"
                    max_distance = parse(Float64, value)
                elseif key == "nthreads"
                    nthreads = parse(Int, value)
                else
                    println("Warning: Unknown argument: --$key")
                end
            catch e
                println("Error parsing argument --$key=$value: $e")
                exit(1)
            end
        end
    end

    # Validate
    if emin >= emax
        println("Error: emin ($emin) must be less than emax ($emax)")
        exit(1)
    end

    if voxel_size <= 0
        println("Error: voxel-size must be positive")
        exit(1)
    end

    if max_distance <= 0
        println("Error: max-distance must be positive")
        exit(1)
    end

    # Run
    filter_single_tracks(cmdir, input_file;
                         emin=emin, emax=emax,
                         voxel_size=voxel_size, max_distance=max_distance,
                         nthreads=nthreads)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
