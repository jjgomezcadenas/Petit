#!/usr/bin/env julia

"""
Batch Single Track Blob Analysis Script

This script performs blob analysis on Xe-137 single-track events.

Usage:
    julia scripts/batch_st_analysis.jl --cmdir=<path> --rblob=<value> --output=<file> [options]

Required Arguments:
    --cmdir       : Data directory (e.g., /path/to/HD5t/precdr)
    --rblob       : Blob radius in mm (e.g., 5.0)
    --output      : Output CSV filename (e.g., xe137_blobs.csv)

Optional Arguments:
    --xedir       : Xe directory name (default: "xe137r2")
    --tag         : File pattern tag, will be expanded to *tag* (default: "st3mm")
    --nmax        : Maximum number of tracks to process, -1 for all (default: -1)
    --nprint      : Progress print interval (default: 100)

Example:
    julia scripts/batch_st_analysis.jl --cmdir=/data/HD5t/precdr --rblob=5.0 --output=xe137_blobs.csv --xedir=xe137v1mm --tag=st1mm --nmax=1000 --nprint=100
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
import Petit: Eff, get_xe137_tracks, blob_analysis, write_blobs_csv

#=============================================================================
# Main Script
=============================================================================#

function main()
    # Initialize variables with defaults
    cmdir = nothing
    rblob = nothing
    output_file = nothing
    xedir = "xe137r2"
    tag = "st3mm"
    nmax = -1
    nprint = 100

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
            elseif key == "output"
                output_file = value
            elseif key == "xedir"
                xedir = value
            elseif key == "tag"
                tag = value
            elseif key == "nmax"
                nmax = parse(Int, value)
            elseif key == "nprint"
                nprint = parse(Int, value)
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
    if isnothing(cmdir) || isnothing(rblob) || isnothing(output_file)
        println("Error: Missing required arguments")
        println()
        println("Usage: julia batch_st_analysis.jl --cmdir=<path> --rblob=<value> --output=<file> [options]")
        println()
        println("Required:")
        println("  --cmdir       : Data directory")
        println("  --rblob       : Blob radius in mm")
        println("  --output      : Output CSV filename")
        println()
        println("Optional:")
        println("  --xedir       : Xe directory name (default: xe137r2)")
        println("  --tag         : File pattern tag, expanded to *tag* (default: st3mm)")
        println("  --nmax        : Max tracks to process, -1 for all (default: -1)")
        println("  --nprint      : Progress interval (default: 100)")
        println()
        println("Example:")
        println("  julia batch_st_analysis.jl --cmdir=/data/HD5t/precdr --rblob=5.0 --output=xe137.csv --xedir=xe137v1mm --tag=st1mm --nmax=1000 --nprint=100")
        exit(1)
    end

    # Validate inputs
    if !isdir(cmdir)
        println("Error: Data directory does not exist: $cmdir")
        exit(1)
    end

    if rblob <= 0
        println("Error: Blob radius must be positive, got: $rblob")
        exit(1)
    end

    if nprint < 0
        println("Error: Progress interval must be non-negative, got: $nprint")
        exit(1)
    end

    # Run analysis
    try
        println("="^60)
        println("Xe-137 BLOB ANALYSIS")
        println("="^60)
        println("  Data directory: $cmdir")
        println("  Xe directory:   $xedir")
        println("  File pattern:   $tag")
        println("  Blob radius:    $rblob mm")
        println("  Output file:    $output_file")
        println("  Max tracks:     $(nmax < 0 ? "all" : nmax)")
        println("  Progress every: $nprint tracks")
        println()

        # Load tracks
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

        # Run blob analysis
        println("Running blob analysis...")
        blobs = blob_analysis(trks.tracks, rblob; nmax=nmax, nprint=nprint)

        # Write results
        df = DataFrame(
            confidence = blobs.confidence,
            trackLength = blobs.trackLength,
            energyKeV = blobs.energyKeV,
            eB1 = blobs.eB1,
            eB2 = blobs.eB2
        )
        CSV.write(output_file, df)
        println()
        println("Results written to: $output_file")
        println("="^60)
        println()
        println("Analysis completed successfully!")
    catch e
        println()
        println("Error during analysis:")
        println(e)
        exit(1)
    end
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
