#!/usr/bin/env julia

"""
Test batch analysis script
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
import Petit: Eff, get_xe137_tracks, blob_analysis

function main()
    println("Starting test...")

    # Just call Petit functions directly from main
    cmdir = "/Users/jjgomezcadenas/Data/HD5t/precdr"
    xedir = "xe137v1mm"
    tag = "st1mm"
    rblob = 10.0
    nmax = 5
    nprint = 0
    output_file = "test_output.csv"

    println("Loading tracks...")
    trks = get_xe137_tracks(cmdir; xedir=xedir, tag=tag)

    println("Running blob analysis...")
    blobs = blob_analysis(trks.tracks, rblob; nmax=nmax, nprint=nprint)

    println("Writing CSV...")
    df = DataFrame(
        confidence = blobs.confidence,
        trackLength = blobs.trackLength,
        energyKeV = blobs.energyKeV,
        eB1 = blobs.eB1,
        eB2 = blobs.eB2
    )
    CSV.write(output_file, df)

    println("Done! Written to $output_file")
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
