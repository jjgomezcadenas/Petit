#!/usr/bin/env julia

"""
Multi-threaded event range processing script for single track event analysis.

This script processes events from ievt to levt using multiple threads,
with each thread writing to a separate output file.
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

# Import Petit module functions
import Petit: event_loop_single_track_mt, get_dataset_dfs, Tracks
import Petit: nof_events, count_events
import Petit: save_tracks_to_hdf5, read_tracks_from_hdf5

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
    julia batch_track_analysis_mt.jl <cmdir> <input_file> <output_base> [options]

Required arguments:
    cmdir           Directory containing the input file
    input_file      Name of the HDF5 input file
    output_base     Base name for output files (without .h5 extension)
                    Each thread will create {output_base}_th_{i}.h5

Optional arguments:
    --ievt          First event to process (default: 1)
    --levt          Last event to process (default: -1, all events)
    --nthreads      Number of threads to use (default: 1)
    --voxel-size    Voxel size in mm (default: 5.0)
    --max-distance  Maximum distance for track clustering in mm (default: 10.0)
    --energy-threshold  Energy threshold in keV (default: 10.0)
    --emin          Minimum event energy in keV (default: 2400.0)
    --emax          Maximum event energy in keV (default: 2500.0)
    --xyc           XY fiducial cut radius (default: 1950.0)
    --zc            Z fiducial cut distance (default: 10.0)

Example:
    # Process events 1000 to 2000 using 4 threads
    julia batch_track_analysis_mt.jl /path/to/data/ input.h5 output_base --ievt=1000 --levt=2000 --nthreads=4

    # Process all events with energy window using 8 threads
    julia batch_track_analysis_mt.jl /path/to/data/ input.h5 output_base --nthreads=8 --emin=2400 --emax=2500

Note:
    - You must start Julia with multiple threads: julia -t auto or julia -t 8
    - The script will check and cap threads at the system maximum
    - Each thread writes to a separate file: output_base_th_1.h5, output_base_th_2.h5, etc.
"""
function main()
    # Default values
    ievt = 1
    levt = -1
    nthreads = 1
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
        println("\nUsage: julia -t <nthreads> batch_track_analysis_mt.jl <cmdir> <input_file> <output_base> [options]")
        println("\nRequired arguments:")
        println("  cmdir           Directory containing the input file")
        println("  input_file      Name of the HDF5 input file")
        println("  output_base     Base name for output files (no .h5 extension)")
        println("\nOptional arguments:")
        println("  --ievt=N                   First event to process (default: 1)")
        println("  --levt=N                   Last event to process (default: -1, all)")
        println("  --nthreads=N               Number of threads to use (default: 1)")
        println("  --voxel-size=X             Voxel size in mm (default: 1.0)")
        println("  --max-distance=X           Max distance in mm (default: 10.0)")
        println("  --energy-threshold=X       Energy threshold in keV (default: 5.0)")
        println("  --emin=X                   Minimum event energy in keV (default: 2400.0)")
        println("  --emax=X                   Maximum event energy in keV (default: 2500.0)")
        println("  --xyc=X                    XY fiducial cut radius (default: 1990.0)")
        println("  --zc=X                     Z fiducial cut distance (default: 5.0)")
        println("\nExample:")
        println("  # Start Julia with 4 threads and process events 1000 to 2000")
        println("  julia -t 4 batch_track_analysis_mt.jl /data/ input.h5 output --ievt=1000 --levt=2000 --nthreads=4")
        println("\nImportant:")
        println("  - You MUST start Julia with multiple threads: julia -t auto or julia -t <N>")
        println("  - Each thread creates a separate file: output_base_th_1.h5, output_base_th_2.h5, etc.")
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
    println("MULTI-THREADED BATCH TRACK ANALYSIS")
    println("="^60)
    println("Configuration:")
    println("  Directory:          $cmdir")
    println("  Input file:         $input_file")
    println("  Output base:        $output_base")
    println("  First event:        $ievt")
    println("  Last event:         $(levt < 0 ? "all" : levt)")
    println("  Threads requested:  $nthreads")
    println("  Threads available:  $(Threads.nthreads())")
    println("  Voxel size:         $voxel_size_mm mm")
    println("  Max distance:       $max_distance_mm mm")
    println("  Energy threshold:   $energy_threshold_kev keV")
    println("  Energy range:       [$emin, $emax] keV")
    println("  XY cut:             $xyc mm")
    println("  Z cut:              $zc mm")
    println("="^60)

    # Check if Julia was started with enough threads
    if Threads.nthreads() < nthreads
        println("\nWarning: Julia started with only $(Threads.nthreads()) threads")
        println("Restart Julia with: julia -t $nthreads <script>")
        println("Proceeding with $(Threads.nthreads()) threads...")
    end

    # Run multi-threaded processing
    println("\nStarting multi-threaded event processing...")
    result = Petit.event_loop_single_track_mt(cmdir, output_base;
                                             input_file=input_file,
                                             ievt=ievt,
                                             levt=levt,
                                             nthreads=nthreads,
                                             voxel_size_mm=voxel_size_mm,
                                             max_distance_mm=max_distance_mm,
                                             energy_threshold_kev=energy_threshold_kev,
                                             emin=emin,
                                             emax=emax,
                                             xyc=xyc,
                                             zc=zc)

    println("\n" * "="^60)
    println("PROCESSING COMPLETE")
    println("="^60)
    println("\nOutput files created:")
    for thread_result in result["thread_results"]
        println("  $(thread_result["output_file"])")
    end

    println("\n" * "="^60)
    println("DONE")
    println("="^60)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

#=
Example usage:
Example command line usage:

```bash
# Start Julia with 4 threads and process events 1000 to 2000
julia -t 4 batch_track_analysis_mt.jl /path/to/data/ input.h5 output_base --ievt=1000 --levt=2000 --nthreads=4

# Use auto thread detection and process with energy window
julia -t auto batch_track_analysis_mt.jl /path/to/data/ input.h5 output_base --emin=2400 --emax=2500 --nthreads=8

# Process all events using 8 threads with custom parameters
julia -t 8 batch_track_analysis_mt.jl /path/to/data/ input.h5 output_base \\
    --nthreads=8 --voxel-size=5.0 --max-distance=10.0 \\
    --energy-threshold=10.0 --emin=2400 --emax=2500
```

Or use as a module:

```julia
# Start Julia with threads: julia -t 8
using Petit

cmdir = "/Users/jjgomezcadenas/Data/HD5t/precdr/copperbkg/bi214/"
input_file = "bi214_copper_endcaps_1.next.h5"
output_base = "bi214_single_tracks_mt"

result = event_loop_single_track_mt(cmdir, output_base;
                                    input_file=input_file,
                                    ievt=1000, levt=2000,
                                    nthreads=4,
                                    emin=2400.0, emax=2500.0)

# Read results from thread files
for i in 1:result["threads_used"]
    filename = "$(output_base)_th_$i.h5"
    tracks, metadata = read_tracks_from_hdf5(joinpath(cmdir, filename))
    println("Thread $i: $(length(tracks)) tracks")
end
```
=#
