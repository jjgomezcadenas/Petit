#!/usr/bin/env julia

"""
Batch processing script for single track event analysis.

This script processes HDF5 files in batches, extracts single-track events,
and saves them to an output HDF5 file for later analysis.
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
#using .Petit

# Import Petit module functions
import Petit: event_loop_single_track, get_dataset_dfs, Tracks
import Petit: nof_events, count_events
import Petit: save_tracks_to_hdf5, read_tracks_from_hdf5

"""
    batch_process_tracks(cmdir::String, input_file::String, output_file::String;
                        nmax::Int=-1, nbatch::Int=100,
                        voxel_size_mm::Float64=5.0,
                        max_distance_mm::Float64=10.0,
                        energy_threshold_kev::Float64=10.0,
                        xyc::Float64=1950.0,
                        zc::Float64=10.0)

Process HDF5 file in batches and extract single-track events.

# Arguments
- `cmdir::String`: Directory containing the input file
- `input_file::String`: Name of the HDF5 input file
- `output_file::String`: Name of the HDF5 output file (will be created in cmdir)
- `nmax::Int=-1`: Maximum number of events to process (-1 means all events)
- `nbatch::Int=100`: Number of events to process per batch
- `voxel_size_mm::Float64=5.0`: Voxel size in mm
- `max_distance_mm::Float64=10.0`: Maximum distance for track clustering
- `energy_threshold_kev::Float64=10.0`: Energy threshold in keV
- `xyc::Float64=1950.0`: XY fiducial cut radius
- `zc::Float64=10.0`: Z fiducial cut distance

# Returns
- `Dict`: Summary statistics of the processing
"""
function batch_process_tracks(cmdir::String, input_file::String, output_file::String;
                              nmax::Int=-1, nbatch::Int=100,
                              voxel_size_mm::Float64=5.0,
                              max_distance_mm::Float64=10.0,
                              energy_threshold_kev::Float64=10.0,
                              xyc::Float64=1950.0,
                              zc::Float64=10.0)

    # Count total events
    println("Counting events in file...")
    input_path = joinpath(cmdir, input_file)

    # Get number of events from configuration
    nevents_config = nof_events(input_path)
    println("Number of events from config: $nevents_config")

    # Get number of events from actual data
    ntot = count_events(cmdir, input_file)
    println("Number of events with hits: $ntot")

    if nevents_config != ntot
        println("Warning: Config events ($nevents_config) != events with hits ($ntot)")
    end

    # Determine number of events to process
    nmx = nmax < 0 ? ntot : min(ntot, nmax)
    println("Events to process: $nmx")

    # Open output HDF5 file
    output_path = joinpath(cmdir, output_file)
    println("Creating output file: $output_path")

    # Statistics tracking
    total_tracks_saved = 0
    batches_processed = 0
    total_processing_time = 0.0
    total_writing_time = 0.0
    events_processed = 0

    h5open(output_path, "w") do fid
        # Store metadata
        attrs(fid)["input_file"] = input_file
        attrs(fid)["nevents_from_config"] = nevents_config
        attrs(fid)["total_events_in_file"] = ntot
        attrs(fid)["events_processed"] = nmx
        attrs(fid)["batch_size"] = nbatch
        attrs(fid)["voxel_size_mm"] = voxel_size_mm
        attrs(fid)["max_distance_mm"] = max_distance_mm
        attrs(fid)["energy_threshold_kev"] = energy_threshold_kev
        attrs(fid)["xyc"] = xyc
        attrs(fid)["zc"] = zc

        # Process in batches
        for ievent in 1:nbatch:nmx
            # Calculate actual batch size for this iteration
            current_batch_size = min(nbatch, nmx - ievent + 1)
            batches_processed += 1

            println("\n" * "="^60)
            println("Processing batch $batches_processed: events $ievent to $(ievent + current_batch_size - 1)")
            println("="^60)

            # Call event_loop_single_track with correct batch size
            proc_start = time()
            tracks = event_loop_single_track(cmdir;
                                           input_file=input_file,
                                           events_to_run=current_batch_size,
                                           voxel_size_mm=voxel_size_mm,
                                           max_distance_mm=max_distance_mm,
                                           energy_threshold_kev=energy_threshold_kev,
                                           xyc=xyc,
                                           zc=zc)
            proc_time = time() - proc_start
            total_processing_time += proc_time
            events_processed += current_batch_size

            # Flush output to ensure progress is visible
            flush(stdout)

            # Save tracks from this batch
            if !isempty(tracks)
                write_time = save_tracks_to_hdf5(tracks, fid, batches_processed)
                total_writing_time += write_time
                total_tracks_saved += length(tracks)
                println("    Total tracks saved so far: $total_tracks_saved")
            else
                println("    No single-track events found in this batch")
            end
            flush(stdout)
        end

        # Store final statistics
        attrs(fid)["total_tracks_saved"] = total_tracks_saved
        attrs(fid)["batches_processed"] = batches_processed
        attrs(fid)["total_processing_time"] = total_processing_time
        attrs(fid)["total_writing_time"] = total_writing_time
    end

    println("\n" * "="^60)
    println("PROCESSING COMPLETE")
    println("="^60)
    println("Total batches processed: $batches_processed")
    println("Total events processed: $events_processed")
    println("Total single tracks saved: $total_tracks_saved")
    println("Output file: $output_path")
    println("\n" * "-"^60)
    println("TIMING STATISTICS")
    println("-"^60)
    println("Total processing time: $(round(total_processing_time, digits=2))s")
    println("Total writing time: $(round(total_writing_time, digits=2))s")
    if events_processed > 0
        println("Average time per event (processing): $(round(1000*total_processing_time/events_processed, digits=2))ms")
    end
    if total_tracks_saved > 0
        println("Average time per track (writing): $(round(1000*total_writing_time/total_tracks_saved, digits=2))ms")
    end

    return Dict(
        "total_events" => ntot,
        "events_processed" => events_processed,
        "batches_processed" => batches_processed,
        "tracks_saved" => total_tracks_saved,
        "output_file" => output_path,
        "processing_time" => total_processing_time,
        "writing_time" => total_writing_time,
        "avg_time_per_event_ms" => events_processed > 0 ? 1000*total_processing_time/events_processed : 0.0,
        "avg_time_per_track_ms" => total_tracks_saved > 0 ? 1000*total_writing_time/total_tracks_saved : 0.0
    )
end

"""
    analyze_track_statistics(output_file::String)

Read tracks from output file and compute statistics.

# Arguments
- `output_file::String`: Path to the HDF5 output file

# Returns
- `DataFrame`: Statistics for each track
"""
function analyze_track_statistics(output_file::String)
    tracks, metadata = read_tracks_from_hdf5(output_file)

    println("\n" * "="^60)
    println("TRACK STATISTICS")
    println("="^60)
    println("Metadata:")
    for (key, val) in metadata
        println("  $key: $val")
    end

    if isempty(tracks)
        println("\nNo tracks found in file.")
        return DataFrame()
    end

    # Compute statistics for each track
    stats = DataFrame(
        track_id = Int[],
        n_voxels = Int[],
        n_vertices = Int[],
        n_edges = Int[],
        n_components = Int[],
        total_energy = Float64[],
        mean_energy = Float64[],
        max_energy = Float64[],
        x_min = Float64[],
        x_max = Float64[],
        y_min = Float64[],
        y_max = Float64[],
        z_min = Float64[],
        z_max = Float64[]
    )

    for (idx, track) in enumerate(tracks)
        push!(stats, (
            track_id = idx,
            n_voxels = nrow(track.voxels),
            n_vertices = nv(track.graph),
            n_edges = ne(track.graph),
            n_components = length(track.components),
            total_energy = sum(track.voxels.energy),
            mean_energy = mean(track.voxels.energy),
            max_energy = maximum(track.voxels.energy),
            x_min = minimum(track.voxels.x),
            x_max = maximum(track.voxels.x),
            y_min = minimum(track.voxels.y),
            y_max = maximum(track.voxels.y),
            z_min = minimum(track.voxels.z),
            z_max = maximum(track.voxels.z)
        ))
    end

    println("\n" * "-"^60)
    println("SUMMARY STATISTICS")
    println("-"^60)
    println("Total tracks: $(nrow(stats))")
    println("Voxels per track: $(round(mean(stats.n_voxels), digits=2)) ± $(round(std(stats.n_voxels), digits=2))")
    println("Energy per track (keV): $(round(mean(stats.total_energy), digits=2)) ± $(round(std(stats.total_energy), digits=2))")
    println("Vertices per track: $(round(mean(stats.n_vertices), digits=2)) ± $(round(std(stats.n_vertices), digits=2))")
    println("Edges per track: $(round(mean(stats.n_edges), digits=2)) ± $(round(std(stats.n_edges), digits=2))")

    return stats
end

"""
    parse_bool(s::String)

Parse a string to boolean value.
"""
function parse_bool(s::String)
    s_lower = lowercase(strip(s))
    if s_lower in ["true", "t", "yes", "y", "1"]
        return true
    elseif s_lower in ["false", "f", "no", "n", "0"]
        return false
    else
        error("Cannot parse '$s' as boolean")
    end
end

"""
    main()

Main function to run from command line.

Usage:
    julia batch_track_analysis.jl <cmdir> <input_file> <output_file> [options]

Required arguments:
    cmdir           Directory containing the input file
    input_file      Name of the HDF5 input file
    output_file     Name of the HDF5 output file

Optional arguments:
    --process       true/false - Run batch processing (default: true)
    --stats         true/false - Run statistics analysis (default: true)
    --nmax          Maximum number of events to process (default: -1, all events)
    --nbatch        Number of events per batch (default: 100)
    --voxel-size    Voxel size in mm (default: 5.0)
    --max-distance  Maximum distance for track clustering in mm (default: 10.0)
    --energy-threshold  Energy threshold in keV (default: 10.0)
    --xyc           XY fiducial cut radius (default: 1950.0)
    --zc            Z fiducial cut distance (default: 10.0)

Example:
    julia batch_track_analysis.jl /path/to/data/ input.h5 output.h5 --nmax=1000 --nbatch=50
"""
function main()
    # Default values
    process = true
    stats = true
    nmax = -1
    nbatch = 10000
    voxel_size_mm = 1.0
    max_distance_mm = 10.0
    energy_threshold_kev = 5.0
    xyc = 1990.0
    zc = 5.0

    # Check minimum required arguments
    if length(ARGS) < 3
        println("Error: Missing required arguments")
        println("\nUsage: julia batch_track_analysis.jl <cmdir> <input_file> <output_file> [options]")
        println("\nRequired arguments:")
        println("  cmdir           Directory containing the input file")
        println("  input_file      Name of the HDF5 input file")
        println("  output_file     Name of the HDF5 output file")
        println("\nOptional arguments:")
        println("  --process=true/false       Run batch processing (default: true)")
        println("  --stats=true/false         Run statistics analysis (default: true)")
        println("  --nmax=N                   Maximum events to process (default: -1)")
        println("  --nbatch=N                 Events per batch (default: 10000)")
        println("  --voxel-size=X             Voxel size in mm (default: 1.0)")
        println("  --max-distance=X           Max distance in mm (default: 10.0)")
        println("  --energy-threshold=X       Energy threshold in keV (default: 5.0)")
        println("  --xyc=X                    XY fiducial cut radius (default: 1990.0)")
        println("  --zc=X                     Z fiducial cut distance (default: 5.0)")
        exit(1)
    end

    # Parse required arguments
    cmdir = ARGS[1]
    input_file = ARGS[2]
    output_file = ARGS[3]

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
                if key == "process"
                    process = parse_bool(value)
                elseif key == "stats"
                    stats = parse_bool(value)
                elseif key == "nmax"
                    nmax = parse(Int, value)
                elseif key == "nbatch"
                    nbatch = parse(Int, value)
                elseif key == "voxel-size"
                    voxel_size_mm = parse(Float64, value)
                elseif key == "max-distance"
                    max_distance_mm = parse(Float64, value)
                elseif key == "energy-threshold"
                    energy_threshold_kev = parse(Float64, value)
                elseif key == "xyc"
                    xyc = parse(Float64, value)
                elseif key == "zc"
                    zc = parse(Float64, value)
                else
                    println("Warning: Unknown argument: --$key")
                end
            catch e
                println("Error parsing argument --$key=$value: $e")
                exit(1)
            end
        end
    end

    # Display configuration
    println("="^60)
    println("BATCH TRACK ANALYSIS")
    println("="^60)
    println("Configuration:")
    println("  Directory:          $cmdir")
    println("  Input file:         $input_file")
    println("  Output file:        $output_file")
    println("  Process data:       $process")
    println("  Run statistics:     $stats")
    println("  Max events:         $nmax")
    println("  Batch size:         $nbatch")
    println("  Voxel size:         $voxel_size_mm mm")
    println("  Max distance:       $max_distance_mm mm")
    println("  Energy threshold:   $energy_threshold_kev keV")
    println("  XY cut:             $xyc mm")
    println("  Z cut:              $zc mm")
    println("="^60)

    # Run processing if requested
    if process
        println("\nStarting batch processing...")
        result = batch_process_tracks(cmdir, input_file, output_file;
                                     nmax=nmax,
                                     nbatch=nbatch,
                                     voxel_size_mm=voxel_size_mm,
                                     max_distance_mm=max_distance_mm,
                                     energy_threshold_kev=energy_threshold_kev,
                                     xyc=xyc,
                                     zc=zc)
        println("\nProcessing result: $result")
    else
        println("\nSkipping batch processing (--process=false)")
    end

    # Run statistics if requested
    if stats
        println("\nStarting statistics analysis...")
        output_path = joinpath(cmdir, output_file)
        if isfile(output_path)
            stats_df = analyze_track_statistics(output_path)
            println("\nFirst 10 tracks:")
            println(first(stats_df, 10))
        else
            println("Error: Output file not found: $output_path")
            println("Run with --process=true first to generate the output file.")
        end
    else
        println("\nSkipping statistics analysis (--stats=false)")
    end

    println("\n" * "="^60)
    println("DONE")
    println("="^60)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# Example usage
"""
Example command line usage:

```bash
# Process all events in batches of 100
julia batch_track_analysis.jl /path/to/data/ input.h5 output.h5

# Process only first 1000 events in batches of 50
julia batch_track_analysis.jl /path/to/data/ input.h5 output.h5 --nmax=1000 --nbatch=50

# Only run statistics on existing output file
julia batch_track_analysis.jl /path/to/data/ input.h5 output.h5 --process=false --stats=true

# Process with custom parameters
julia batch_track_analysis.jl /path/to/data/ input.h5 output.h5 \
    --nmax=-1 --nbatch=100 --voxel-size=5.0 --max-distance=10.0 \
    --energy-threshold=10.0 --xyc=1950.0 --zc=10.0
```

Or use as a module:

```julia
using Petit

cmdir = "/Users/jjgomezcadenas/Data/HD5t/precdr/copperbkg/bi214/"
input_file = "bi214_copper_endcaps_1.next.h5"
output_file = "bi214_single_tracks.h5"

result = batch_process_tracks(cmdir, input_file, output_file;
                              nmax=-1, nbatch=100)

stats = analyze_track_statistics(joinpath(cmdir, output_file))
tracks, metadata = read_tracks_from_hdf5(joinpath(cmdir, output_file))
```
"""
