#!/usr/bin/env julia

"""
Event range processing script for single track event analysis.

This script processes events from ievt to levt (skipping the first ievt-1 events),
extracts single-track events, and saves them to an output HDF5 file.
Progress is shown with a progress bar.
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
using ProgressMeter

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))

# Import Petit module functions
import Petit: event_loop_single_track2, get_dataset_dfs, Tracks
import Petit: nof_events, count_events
import Petit: save_tracks_to_hdf5, read_tracks_from_hdf5

"""
    process_event_range(cmdir::String, input_file::String, output_file::String;
                       ievt::Int=1, levt::Int=-1,
                       voxel_size_mm::Float64=5.0,
                       max_distance_mm::Float64=10.0,
                       energy_threshold_kev::Float64=10.0,
                       emin::Float64=2400.0,
                       emax::Float64=2500.0,
                       xyc::Float64=1950.0,
                       zc::Float64=10.0)

Process HDF5 file from event ievt to levt and extract single-track events.

The script:
1. Loads data once for efficiency
2. Skips the first (ievt - 1) events
3. Processes events from ievt to levt (inclusive)
4. Writes all processed events to a single HDF5 file (no batching)
5. Shows progress with a progress bar

# Arguments
- `cmdir::String`: Directory containing the input file
- `input_file::String`: Name of the HDF5 input file
- `output_file::String`: Name of the HDF5 output file (will be created in cmdir)
- `ievt::Int=1`: First event to process (1-indexed)
- `levt::Int=-1`: Last event to process (-1 means all events)
- `voxel_size_mm::Float64=5.0`: Voxel size in mm
- `max_distance_mm::Float64=10.0`: Maximum distance for track clustering
- `energy_threshold_kev::Float64=10.0`: Energy threshold in keV
- `emin::Float64=2400.0`: Minimum event energy in keV
- `emax::Float64=2500.0`: Maximum event energy in keV
- `xyc::Float64=1950.0`: XY fiducial cut radius
- `zc::Float64=10.0`: Z fiducial cut distance

# Returns
- `Dict`: Summary statistics of the processing

# Example
```julia
# Process events 1000 to 2000 (skips first 999, processes 1001 events)
process_event_range("/path/to/data", "input.h5", "output.h5", ievt=1000, levt=2000)
```
"""
function process_event_range(cmdir::String, input_file::String, output_file::String;
                             ievt::Int=1, levt::Int=-1,
                             voxel_size_mm::Float64=5.0,
                             max_distance_mm::Float64=10.0,
                             energy_threshold_kev::Float64=10.0,
                             emin::Float64=2400.0,
                             emax::Float64=2500.0,
                             xyc::Float64=1950.0,
                             zc::Float64=10.0)

    input_path = joinpath(cmdir, input_file)

    # Get number of events from configuration
    nevents_config = nof_events(input_path)
    println("Number of events from config: $nevents_config")

    # Load data once (this is the slow operation)
    println(" Loading data from: $input_path")
    dfs = get_dataset_dfs(input_path)
    hitsdf = dfs["hits"]

    # Count events from already-loaded data (fast)
    println("Counting events from loaded data...")
    ntot = length(unique(hitsdf.event_id))
    println("Number of events with hits: $ntot")

    if nevents_config != ntot
        println("Warning: Config events ($nevents_config) != events with hits ($ntot)")
    end

    # Validate event range
    if ievt < 1
        error("ievt must be >= 1, got $ievt")
    end
    if ievt > ntot
        error("ievt ($ievt) exceeds total events ($ntot)")
    end

    # Determine last event to process
    last_event = levt < 0 ? ntot : min(levt, ntot)

    if last_event < ievt
        error("levt ($last_event) must be >= ievt ($ievt)")
    end

    # Calculate number of events to process
    nevents_to_process = last_event - ievt + 1

    println("\nEvent processing range:")
    println("  First event to process: $ievt")
    println("  Last event to process:  $last_event")
    println("  Events to skip:         $(ievt - 1)")
    println("  Events to process:      $nevents_to_process")

    # Open output HDF5 file
    output_path = joinpath(cmdir, output_file)
    println("\nCreating output file: $output_path")

    # Statistics tracking
    total_tracks_saved = 0
    total_processing_time = 0.0
    total_writing_time = 0.0

    # Open output file and store metadata
    h5open(output_path, "w") do fid
        # Store metadata
        attrs(fid)["input_file"] = input_file
        attrs(fid)["nevents_from_config"] = nevents_config
        attrs(fid)["total_events_in_file"] = ntot
        attrs(fid)["first_event_processed"] = ievt
        attrs(fid)["last_event_processed"] = last_event
        attrs(fid)["events_skipped"] = ievt - 1
        attrs(fid)["events_processed"] = nevents_to_process
        attrs(fid)["voxel_size_mm"] = voxel_size_mm
        attrs(fid)["max_distance_mm"] = max_distance_mm
        attrs(fid)["energy_threshold_kev"] = energy_threshold_kev
        attrs(fid)["emin"] = emin
        attrs(fid)["emax"] = emax
        attrs(fid)["xyc"] = xyc
        attrs(fid)["zc"] = zc

        println("\n" * "="^60)
        println("PROCESSING EVENTS $ievt TO $last_event")
        println("="^60)

        # Process all events in the range at once
        proc_start = time()

        # Call analysis_loop_single_track2 directly with already-loaded data
        # This avoids reloading the data which is very slow
        println(" Starting analysis from event $ievt, processing $nevents_to_process events")
        println(" Energy range filter: [$emin, $emax] keV")
        tracks = Petit.analysis_loop_single_track2(hitsdf;
                                                   events_to_run=nevents_to_process,
                                                   initial_event=ievt,
                                                   show_progress=true,
                                                   voxel_size_mm=voxel_size_mm,
                                                   max_distance_mm=max_distance_mm,
                                                   energy_threshold_kev=energy_threshold_kev,
                                                   emin=emin,
                                                   emax=emax)

        proc_time = time() - proc_start
        total_processing_time = proc_time

        # Save all tracks
        println("\nSaving tracks to output file...")
        if !isempty(tracks)
            write_start = time()
            write_time = save_tracks_to_hdf5(tracks, fid, 1)
            total_writing_time = write_time
            total_tracks_saved = length(tracks)
            println("Saved $total_tracks_saved tracks")
        else
            println("No single-track events found in the processed range")
        end

        # Store final statistics
        attrs(fid)["total_tracks_saved"] = total_tracks_saved
        attrs(fid)["total_processing_time"] = total_processing_time
        attrs(fid)["total_writing_time"] = total_writing_time
    end

    println("\n" * "="^60)
    println("PROCESSING COMPLETE")
    println("="^60)
    println("Event range: $ievt to $last_event")
    println("Events skipped: $(ievt - 1)")
    println("Events processed: $nevents_to_process")
    println("Single tracks saved: $total_tracks_saved")
    println("Output file: $output_path")
    println("\n" * "-"^60)
    println("TIMING STATISTICS")
    println("-"^60)
    println("Total processing time: $(round(total_processing_time, digits=2))s")
    println("Total writing time: $(round(total_writing_time, digits=2))s")
    if nevents_to_process > 0
        println("Average time per event: $(round(1000*total_processing_time/nevents_to_process, digits=2))ms")
    end
    if total_tracks_saved > 0
        println("Average time per track (writing): $(round(1000*total_writing_time/total_tracks_saved, digits=2))ms")
    end

    return Dict(
        "total_events_in_file" => ntot,
        "first_event" => ievt,
        "last_event" => last_event,
        "events_skipped" => ievt - 1,
        "events_processed" => nevents_to_process,
        "tracks_saved" => total_tracks_saved,
        "output_file" => output_path,
        "processing_time" => total_processing_time,
        "writing_time" => total_writing_time,
        "avg_time_per_event_ms" => nevents_to_process > 0 ? 1000*total_processing_time/nevents_to_process : 0.0,
        "avg_time_per_track_ms" => total_tracks_saved > 0 ? 1000*total_writing_time/total_tracks_saved : 0.0
    )
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
    julia batch_track_analysis2.jl <cmdir> <input_file> <output_file> [options]

Required arguments:
    cmdir           Directory containing the input file
    input_file      Name of the HDF5 input file
    output_file     Name of the HDF5 output file

Optional arguments:
    --ievt          First event to process (default: 1)
    --levt          Last event to process (default: -1, all events)
    --voxel-size    Voxel size in mm (default: 5.0)
    --max-distance  Maximum distance for track clustering in mm (default: 10.0)
    --energy-threshold  Energy threshold in keV (default: 10.0)
    --emin          Minimum event energy in keV (default: 2400.0)
    --emax          Maximum event energy in keV (default: 2500.0)
    --xyc           XY fiducial cut radius (default: 1950.0)
    --zc            Z fiducial cut distance (default: 10.0)

Example:
    # Process events 1000 to 2000 (skips first 999, processes 1001 events)
    julia batch_track_analysis2.jl /path/to/data/ input.h5 output.h5 --ievt=1000 --levt=2000

    # Process events 500 to end with energy window
    julia batch_track_analysis2.jl /path/to/data/ input.h5 output.h5 --ievt=500 --emin=2400 --emax=2500

    # Process all events from the beginning
    julia batch_track_analysis2.jl /path/to/data/ input.h5 output.h5
"""
function main()
    # Default values
    ievt = 1
    levt = -1
    voxel_size_mm = 1.0
    max_distance_mm = 10.0
    energy_threshold_kev = 5.0
    emin = 2400.0
    emax = 2500.0
    xyc = 1990.0
    zc = 5.0

    # Check minimum required arguments
    if length(ARGS) < 3
        println("Error: Missing required arguments")
        println("\nUsage: julia batch_track_analysis2.jl <cmdir> <input_file> <output_file> [options]")
        println("\nRequired arguments:")
        println("  cmdir           Directory containing the input file")
        println("  input_file      Name of the HDF5 input file")
        println("  output_file     Name of the HDF5 output file")
        println("\nOptional arguments:")
        println("  --ievt=N                   First event to process (default: 1)")
        println("  --levt=N                   Last event to process (default: -1, all)")
        println("  --voxel-size=X             Voxel size in mm (default: 1.0)")
        println("  --max-distance=X           Max distance in mm (default: 10.0)")
        println("  --energy-threshold=X       Energy threshold in keV (default: 5.0)")
        println("  --emin=X                   Minimum event energy in keV (default: 2400.0)")
        println("  --emax=X                   Maximum event energy in keV (default: 2500.0)")
        println("  --xyc=X                    XY fiducial cut radius (default: 1990.0)")
        println("  --zc=X                     Z fiducial cut distance (default: 5.0)")
        println("\nExample:")
        println("  # Process events 1000 to 2000 (skips first 999 events)")
        println("  julia batch_track_analysis2.jl /data/ input.h5 output.h5 --ievt=1000 --levt=2000")
        println("  # With energy window")
        println("  julia batch_track_analysis2.jl /data/ input.h5 output.h5 --emin=2400 --emax=2500")
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
                if key == "ievt"
                    ievt = parse(Int, value)
                elseif key == "levt"
                    levt = parse(Int, value)
                elseif key == "voxel-size"
                    voxel_size_mm = parse(Float64, value)
                elseif key == "max-distance"
                    max_distance_mm = parse(Float64, value)
                elseif key == "energy-threshold"
                    energy_threshold_kev = parse(Float64, value)
                elseif key == "emin"
                    emin = parse(Float64, value)
                elseif key == "emax"
                    emax = parse(Float64, value)
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
    println("EVENT RANGE TRACK ANALYSIS")
    println("="^60)
    println("Configuration:")
    println("  Directory:          $cmdir")
    println("  Input file:         $input_file")
    println("  Output file:        $output_file")
    println("  First event:        $ievt")
    println("  Last event:         $(levt < 0 ? "all" : levt)")
    println("  Voxel size:         $voxel_size_mm mm")
    println("  Max distance:       $max_distance_mm mm")
    println("  Energy threshold:   $energy_threshold_kev keV")
    println("  Energy range:       [$emin, $emax] keV")
    println("  XY cut:             $xyc mm")
    println("  Z cut:              $zc mm")
    println("="^60)

    # Run processing
    println("\nStarting event range processing...")
    result = process_event_range(cmdir, input_file, output_file;
                                 ievt=ievt,
                                 levt=levt,
                                 voxel_size_mm=voxel_size_mm,
                                 max_distance_mm=max_distance_mm,
                                 energy_threshold_kev=energy_threshold_kev,
                                 emin=emin,
                                 emax=emax,
                                 xyc=xyc,
                                 zc=zc)

    println("\nProcessing result:")
    for (key, val) in result
        println("  $key: $val")
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
# Process events 1000 to 2000 (skips first 999, processes 1001 events)
julia batch_track_analysis2.jl /path/to/data/ input.h5 output.h5 --ievt=1000 --levt=2000

# Process events 500 to end with energy window
julia batch_track_analysis2.jl /path/to/data/ input.h5 output.h5 --ievt=500 --emin=2400 --emax=2500

# Process all events from beginning with default energy window
julia batch_track_analysis2.jl /path/to/data/ input.h5 output.h5

# Process events 10000 to 15000 with custom parameters and energy window
julia batch_track_analysis2.jl /path/to/data/ input.h5 output.h5 \\
    --ievt=10000 --levt=15000 --voxel-size=5.0 --max-distance=10.0 \\
    --energy-threshold=10.0 --emin=2400 --emax=2500 --xyc=1950.0 --zc=10.0
```

Or use as a module:

```julia
using Petit

cmdir = "/Users/jjgomezcadenas/Data/HD5t/precdr/copperbkg/bi214/"
input_file = "bi214_copper_endcaps_1.next.h5"
output_file = "bi214_single_tracks_1000_2000.h5"

result = process_event_range(cmdir, input_file, output_file;
                             ievt=1000, levt=2000,
                             emin=2400.0, emax=2500.0)

tracks, metadata = read_tracks_from_hdf5(joinpath(cmdir, output_file))
```
"""
