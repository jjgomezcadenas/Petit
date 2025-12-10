#!/usr/bin/env julia

"""
Multi-threaded Batch Single Track Blob Analysis Script

This script performs blob analysis on single-track events using multiple threads.
Supports both fixed and variable radius blob analysis.

Usage:
    julia -t <N> scripts/batch_st_analysis_mt.jl --cmdir=<path> --outdir=<dir> [options]

Required Arguments:
    --cmdir       : Data directory (e.g., /path/to/HD5t/precdr)
    --outdir      : Output directory name (will be created)

Optional Arguments (General):
    --xedir       : "xenon" directory name (where the single tracks ) (default: "xe137r2")
    --tag         : File pattern tag, will be expanded to *tag* (default: "st3mm")
    --nmax        : Maximum number of tracks to process, -1 for all (default: -1)
    --nprint      : Progress print interval per thread (default: 100)
    --nthreads    : Number of threads to use (default: 1)
    --variable    : Use variable radius mode (true/false, default: false)

Fixed Radius Mode Arguments (when --variable=false):
    --rblob       : Fixed blob radius in mm (required if --variable=false)

Variable Radius Mode Arguments (when --variable=true):
    --seed        : Initial seed radius in mm (default: 3.0)
    --step        : Radius increment step in mm (default: 1.0)
    --maxrad      : Maximum radius in mm (default: 10.0)
    --threshold   : Relative energy change threshold (default: 0.05)

Examples:
    # Fixed radius mode
    julia -t 4 scripts/batch_st_analysis_mt.jl --cmdir=/data/HD5t/precdr --rblob=5.0 --outdir=xev1mm --xedir=xe137v1mm --tag=st1mm --nthreads=4

    # Variable radius mode
    julia -t 4 scripts/batch_st_analysis_mt.jl --cmdir=/data/HD5t/precdr --outdir=xev1mm_var --xedir=xe137v1mm --tag=st1mm --variable=true --seed=3.0 --step=1.0 --maxrad=10.0 --threshold=0.05 --nthreads=4
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using Printf
using DataFrames
using CSV
using Glob

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))

# Import Petit module functions
import Petit: Eff, get_xe137_tracks, blob_analysis, blob_analysis_variable
import Petit: get_optimal_threads, split_events_for_threads

#=============================================================================
# Main Script
=============================================================================#

function main()
    # Initialize variables with defaults
    cmdir = nothing
    rblob = nothing
    outdir = nothing
    xedir = "xe137r2"
    tag = "st3mm"
    nmax = -1
    nprint = 100
    nthreads = 1
    use_variable = false
    # Variable radius parameters
    seed_radius = 3.0
    step = 1.0
    max_radius = 10.0
    threshold = 0.05

    # Parse named command line arguments
    for arg in ARGS
        if startswith(arg, "--")
            parts = split(arg[3:end], '=', limit=2)
            if length(parts) != 2
                println("Error: Invalid argument format: $arg")
                println("Expected format: --name=value")
                exit(1)
            end

            key, value = parts

            if key == "cmdir"
                cmdir = value
            elseif key == "rblob"
                rblob = parse(Float64, value)
            elseif key == "outdir"
                outdir = value
            elseif key == "xedir"
                xedir = value
            elseif key == "tag"
                tag = value
            elseif key == "nmax"
                nmax = parse(Int, value)
            elseif key == "nprint"
                nprint = parse(Int, value)
            elseif key == "nthreads"
                nthreads = parse(Int, value)
            elseif key == "variable"
                use_variable = lowercase(value) in ["true", "t", "yes", "y", "1"]
            elseif key == "seed"
                seed_radius = parse(Float64, value)
            elseif key == "step"
                step = parse(Float64, value)
            elseif key == "maxrad"
                max_radius = parse(Float64, value)
            elseif key == "threshold"
                threshold = parse(Float64, value)
            else
                println("Warning: Unknown argument: --$key")
            end
        else
            println("Error: Arguments must start with --")
            println("Got: $arg")
            exit(1)
        end
    end

    # Check required arguments
    if isnothing(cmdir) || isnothing(outdir)
        println("Error: Missing required arguments")
        println()
        println("Usage: julia -t <N> batch_st_analysis_mt.jl --cmdir=<path> --outdir=<dir> [options]")
        println()
        println("Required:")
        println("  --cmdir       : Data directory")
        println("  --outdir      : Output directory name")
        println()
        println("Fixed radius mode (--variable=false, default):")
        println("  --rblob       : Blob radius in mm (required)")
        println()
        println("Variable radius mode (--variable=true):")
        println("  --seed        : Initial seed radius in mm (default: 3.0)")
        println("  --step        : Radius increment in mm (default: 1.0)")
        println("  --maxrad      : Maximum radius in mm (default: 10.0)")
        println("  --threshold   : Convergence threshold (default: 0.05)")
        println()
        println("General options:")
        println("  --xedir       : Xe directory name (default: xe137r2)")
        println("  --tag         : File pattern tag (default: st3mm)")
        println("  --nmax        : Max tracks to process, -1 for all (default: -1)")
        println("  --nprint      : Progress interval per thread (default: 100)")
        println("  --nthreads    : Number of threads to use (default: 1)")
        println()
        println("Examples:")
        println("  # Fixed radius")
        println("  julia -t 4 batch_st_analysis_mt.jl --cmdir=/data --rblob=5.0 --outdir=out --nthreads=4")
        println("  # Variable radius")
        println("  julia -t 4 batch_st_analysis_mt.jl --cmdir=/data --outdir=out --variable=true --nthreads=4")
        exit(1)
    end

    # Check mode-specific requirements
    if !use_variable && isnothing(rblob)
        println("Error: --rblob is required for fixed radius mode")
        println("Use --variable=true for variable radius mode, or specify --rblob=<value>")
        exit(1)
    end

    # Validate inputs
    if !isdir(cmdir)
        println("Error: Data directory does not exist: $cmdir")
        exit(1)
    end

    if !use_variable && rblob <= 0
        println("Error: Blob radius must be positive, got: $rblob")
        exit(1)
    end

    if use_variable
        if seed_radius <= 0 || step <= 0 || max_radius <= 0
            println("Error: seed, step, and maxrad must be positive")
            exit(1)
        end
        if seed_radius >= max_radius
            println("Error: seed radius must be less than max radius")
            exit(1)
        end
        if threshold <= 0 || threshold >= 1
            println("Error: threshold must be between 0 and 1")
            exit(1)
        end
    end

    if nprint < 0
        println("Error: Progress interval must be non-negative, got: $nprint")
        exit(1)
    end

    # Get optimal number of threads
    optimal_nthreads = get_optimal_threads(nthreads)

    # Check if Julia was started with enough threads
    if Threads.nthreads() < optimal_nthreads
        println("\nWarning: Julia started with only $(Threads.nthreads()) threads")
        println("Restart Julia with: julia -t $optimal_nthreads <script>")
        println("Proceeding with $(Threads.nthreads()) threads...")
        optimal_nthreads = Threads.nthreads()
    end

    # Create output directory
    if !isdir(outdir)
        mkdir(outdir)
        println("Created output directory: $outdir")
    else
        println("Using existing output directory: $outdir")
    end

    # Run analysis
    try
        println("="^60)
        println("Xe-137 BLOB ANALYSIS (MULTI-THREADED)")
        println("="^60)
        println("  Data directory: $cmdir")
        println("  Xe directory:   $xedir")
        println("  File pattern:   $tag")
        println("  Mode:           $(use_variable ? "Variable radius" : "Fixed radius")")
        if use_variable
            println("  Seed radius:    $seed_radius mm")
            println("  Step size:      $step mm")
            println("  Max radius:     $max_radius mm")
            println("  Threshold:      $(threshold*100)%")
        else
            println("  Blob radius:    $rblob mm")
        end
        println("  Output dir:     $outdir")
        println("  Max tracks:     $(nmax < 0 ? "all" : nmax)")
        println("  Progress every: $nprint tracks")
        println("  Threads:        $optimal_nthreads")
        println()

        # Load tracks (shared across all threads)
        println("Loading Xe-137 tracks...")
        trks = get_xe137_tracks(cmdir; xedir=xedir, tag=tag)

        # Calculate efficiency
        effXe = Eff(trks.n1trk/trks.ntot, 1.0, 1.0, 1.0)

        # Print statistics
        println()
        println("Statistics:")
        println("  Total events generated:    $(trks.ntot)")
        println("  Single track events:       $(trks.n1trk)")
        println("  Selection efficiency:      $(@sprintf("%.2e", effXe.eff1tr))")
        println()

        # Determine how many tracks to process
        total_tracks = length(trks.tracks)
        ntracks_to_process = nmax > 0 ? min(nmax, total_tracks) : total_tracks

        println("Processing $ntracks_to_process out of $total_tracks tracks")
        println()

        # Split tracks among threads
        thread_ranges = split_events_for_threads(ntracks_to_process, optimal_nthreads, 1)

        println("Track distribution across threads:")
        for (i, (start_track, num_tracks)) in enumerate(thread_ranges)
            end_track = start_track + num_tracks - 1
            println("  Thread $i: tracks $start_track to $end_track ($num_tracks tracks)")
        end
        println()

        # Run blob analysis in parallel
        println("Running multi-threaded blob analysis...")
        println()

        Threads.@threads for thread_id in 1:length(thread_ranges)
            start_track, num_tracks = thread_ranges[thread_id]
            end_track = start_track + num_tracks - 1

            println("Thread $thread_id: Starting analysis of tracks $start_track to $end_track")

            # Process this thread's subset of tracks
            thread_tracks = trks.tracks[start_track:end_track]

            if use_variable
                # Variable radius mode
                blobs = blob_analysis_variable(
                    thread_tracks;
                    seed_radius=seed_radius,
                    step=step,
                    max_radius=max_radius,
                    threshold=threshold,
                    nmax=-1,
                    nprint=nprint
                )

                # Write results with radius columns
                output_file = joinpath(outdir, "$(outdir)_th_$(thread_id).csv")
                df = DataFrame(
                    confidence = blobs.confidence,
                    trackLength = blobs.trackLength,
                    energyKeV = blobs.energyKeV,
                    eB1 = blobs.eB1,
                    eB2 = blobs.eB2,
                    rB1 = blobs.rB1,
                    rB2 = blobs.rB2
                )
                CSV.write(output_file, df)
            else
                # Fixed radius mode
                blobs = blob_analysis(thread_tracks, rblob; nmax=-1, nprint=nprint)

                # Write results without radius columns
                output_file = joinpath(outdir, "$(outdir)_th_$(thread_id).csv")
                df = DataFrame(
                    confidence = blobs.confidence,
                    trackLength = blobs.trackLength,
                    energyKeV = blobs.energyKeV,
                    eB1 = blobs.eB1,
                    eB2 = blobs.eB2
                )
                CSV.write(output_file, df)
            end

            println("Thread $thread_id: Completed! Written to $output_file")
        end

        println()
        println("="^60)
        println("All threads completed successfully!")
        println("Results written to directory: $outdir")
        println("  Files: $(outdir)_th_1.csv to $(outdir)_th_$(length(thread_ranges)).csv")
        println("="^60)

    catch e
        println()
        println("Error during analysis:")
        println(e)
        if isa(e, Exception)
            for (exc, bt) in Base.catch_stack()
                showerror(stdout, exc, bt)
                println()
            end
        end
        exit(1)
    end
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
