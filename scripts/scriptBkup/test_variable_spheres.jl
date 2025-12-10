#!/usr/bin/env julia

"""
Test script for energy_in_variable_spheres_around_extremes function.

Analyzes the expansion process for real track data to tune parameters.
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)

using Printf
using DataFrames
using Statistics

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))

# Import Petit module functions
import Petit: read_tracks_from_hdf5, walk_track_from_extremes
import Petit: energy_in_spheres_around_extremes, energy_in_variable_spheres_around_extremes

function analyze_sphere_expansion(track_file::String;
                                  nevents::Int=10,
                                  seed_radius::Float64=3.0,
                                  step::Float64=1.0,
                                  max_radius::Float64=10.0,
                                  threshold::Float64=0.05)

    println("="^70)
    println("VARIABLE SPHERE EXPANSION ANALYSIS")
    println("="^70)
    println("File: $track_file")
    println("Events to analyze: $nevents")
    println("Parameters:")
    println("  Seed radius:  $seed_radius mm")
    println("  Step size:    $step mm")
    println("  Max radius:   $max_radius mm")
    println("  Threshold:    $(threshold*100)%")
    println("="^70)
    println()

    # Read tracks
    println("Reading tracks from file...")
    tracks, metadata = read_tracks_from_hdf5(track_file)
    println("Loaded $(length(tracks)) tracks")
    println()

    # Limit to requested number of events
    ntracks = min(nevents, length(tracks))

    # Statistics for comparison
    fixed_radius = 5.0  # Compare with a fixed radius

    results = []

    for i in 1:ntracks
        track = tracks[i]

        println("\n" * "─"^70)
        println("EVENT $i / $ntracks")
        println("─"^70)
        println("Track has $(nrow(track.voxels)) voxels")

        # Walk track to find extremes
        walk_result = walk_track_from_extremes(track)

        if isnothing(walk_result.extremes[1])
            println("⚠ No valid extremes found, skipping")
            continue
        end

        println("Track length: $(round(walk_result.total_length, digits=2)) mm")
        println("Confidence: $(round(walk_result.confidence, digits=3))")
        println()

        # Calculate with variable spheres
        println("🔍 VARIABLE SPHERE ANALYSIS:")
        var_result = energy_in_variable_spheres_around_extremes(
            track, walk_result;
            seed_radius=seed_radius,
            step=step,
            max_radius=max_radius,
            threshold=threshold
        )

        # Calculate with fixed radius for comparison
        fixed_result = energy_in_spheres_around_extremes(track, walk_result, fixed_radius)

        # Print blob1 expansion history
        println("\n📊 BLOB 1 (higher energy) expansion:")
        println("  Final radius: $(round(var_result.blob1_radius, digits=2)) mm")
        println("  Final energy: $(round(var_result.blob1_energy*1000, digits=2)) keV")
        println("  Voxels: $(var_result.blob1_voxel_count)")
        println("\n  Expansion history:")
        println("  " * "-"^60)
        println("  Radius (mm)  Energy (keV)  ΔE (keV)  Rel. Change (%)")
        println("  " * "-"^60)

        prev_energy = 0.0
        for (radius, energy) in var_result.blob1_history
            energy_kev = energy * 1000
            delta_e = energy_kev - prev_energy
            rel_change = prev_energy > 0 ? abs(delta_e / energy_kev) * 100 : 0.0

            marker = ""
            if radius == var_result.blob1_radius
                marker = " ← STOPPED"
            end

            @printf("  %8.2f     %9.2f     %7.2f      %6.2f%s\n",
                    radius, energy_kev, delta_e, rel_change, marker)

            prev_energy = energy_kev
        end

        # Print blob2 expansion history
        println("\n📊 BLOB 2 (lower energy) expansion:")
        println("  Final radius: $(round(var_result.blob2_radius, digits=2)) mm")
        println("  Final energy: $(round(var_result.blob2_energy*1000, digits=2)) keV")
        println("  Voxels: $(var_result.blob2_voxel_count)")
        println("\n  Expansion history:")
        println("  " * "-"^60)
        println("  Radius (mm)  Energy (keV)  ΔE (keV)  Rel. Change (%)")
        println("  " * "-"^60)

        prev_energy = 0.0
        for (radius, energy) in var_result.blob2_history
            energy_kev = energy * 1000
            delta_e = energy_kev - prev_energy
            rel_change = prev_energy > 0 ? abs(delta_e / energy_kev) * 100 : 0.0

            marker = ""
            if radius == var_result.blob2_radius
                marker = " ← STOPPED"
            end

            @printf("  %8.2f     %9.2f     %7.2f      %6.2f%s\n",
                    radius, energy_kev, delta_e, rel_change, marker)

            prev_energy = energy_kev
        end

        # Comparison with fixed radius
        println("\n📏 COMPARISON (Variable vs Fixed radius = $fixed_radius mm):")
        blob1_var_kev = var_result.blob1_energy * 1000
        blob2_var_kev = var_result.blob2_energy * 1000
        blob1_fix_kev = fixed_result.blob1_energy * 1000
        blob2_fix_kev = fixed_result.blob2_energy * 1000

        blob1_diff = blob1_var_kev - blob1_fix_kev
        blob2_diff = blob2_var_kev - blob2_fix_kev
        blob1_diff_pct = (blob1_diff / blob1_fix_kev) * 100
        blob2_diff_pct = (blob2_diff / blob2_fix_kev) * 100

        println("  Blob1 - Variable: $(round(blob1_var_kev, digits=2)) keV at r=$(round(var_result.blob1_radius, digits=2)) mm")
        println("  Blob1 - Fixed:    $(round(blob1_fix_kev, digits=2)) keV at r=$fixed_radius mm")
        println("  Difference:       $(round(blob1_diff, digits=2)) keV ($(round(blob1_diff_pct, digits=1))%)")
        println()
        println("  Blob2 - Variable: $(round(blob2_var_kev, digits=2)) keV at r=$(round(var_result.blob2_radius, digits=2)) mm")
        println("  Blob2 - Fixed:    $(round(blob2_fix_kev, digits=2)) keV at r=$fixed_radius mm")
        println("  Difference:       $(round(blob2_diff, digits=2)) keV ($(round(blob2_diff_pct, digits=1))%)")

        # Store results for summary
        push!(results, Dict(
            "event" => i,
            "blob1_radius" => var_result.blob1_radius,
            "blob2_radius" => var_result.blob2_radius,
            "blob1_energy_kev" => blob1_var_kev,
            "blob2_energy_kev" => blob2_var_kev,
            "blob1_iterations" => length(var_result.blob1_history),
            "blob2_iterations" => length(var_result.blob2_history),
            "blob1_diff_pct" => blob1_diff_pct,
            "blob2_diff_pct" => blob2_diff_pct
        ))
    end

    # Summary statistics
    if !isempty(results)
        println("\n" * "="^70)
        println("SUMMARY STATISTICS")
        println("="^70)

        blob1_radii = [r["blob1_radius"] for r in results]
        blob2_radii = [r["blob2_radius"] for r in results]
        blob1_iters = [r["blob1_iterations"] for r in results]
        blob2_iters = [r["blob2_iterations"] for r in results]
        blob1_diffs = [r["blob1_diff_pct"] for r in results]
        blob2_diffs = [r["blob2_diff_pct"] for r in results]

        println("\nOptimal radii found:")
        println("  Blob1 - Mean: $(round(mean(blob1_radii), digits=2)) mm, Std: $(round(std(blob1_radii), digits=2)) mm")
        println("  Blob1 - Min: $(round(minimum(blob1_radii), digits=2)) mm, Max: $(round(maximum(blob1_radii), digits=2)) mm")
        println("  Blob2 - Mean: $(round(mean(blob2_radii), digits=2)) mm, Std: $(round(std(blob2_radii), digits=2)) mm")
        println("  Blob2 - Min: $(round(minimum(blob2_radii), digits=2)) mm, Max: $(round(maximum(blob2_radii), digits=2)) mm")

        println("\nIterations to convergence:")
        println("  Blob1 - Mean: $(round(mean(blob1_iters), digits=1)), Max: $(maximum(blob1_iters))")
        println("  Blob2 - Mean: $(round(mean(blob2_iters), digits=1)), Max: $(maximum(blob2_iters))")

        println("\nEnergy difference vs fixed radius ($fixed_radius mm):")
        println("  Blob1 - Mean: $(round(mean(blob1_diffs), digits=1))%, Max: $(round(maximum(abs.(blob1_diffs)), digits=1))%")
        println("  Blob2 - Mean: $(round(mean(blob2_diffs), digits=1))%, Max: $(round(maximum(abs.(blob2_diffs)), digits=1))%")

        println("\n" * "="^70)
        println("RECOMMENDATIONS:")
        println("="^70)

        avg_radius = mean([blob1_radii; blob2_radii])
        if avg_radius < seed_radius + step
            println("⚠ Most blobs converge quickly. Consider:")
            println("  - Smaller seed_radius (current: $seed_radius mm)")
            println("  - Smaller step size (current: $step mm)")
        elseif avg_radius > max_radius - step
            println("⚠ Many blobs reach max_radius. Consider:")
            println("  - Larger max_radius (current: $max_radius mm)")
        else
            println("✓ Parameters seem reasonable for this dataset")
        end

        avg_iters = mean([blob1_iters; blob2_iters])
        if avg_iters > 5
            println("  - Consider larger step size for faster convergence (current: $step mm)")
        elseif avg_iters < 3
            println("  - Consider smaller step size for finer resolution (current: $step mm)")
        end

        max_diff = maximum([abs.(blob1_diffs); abs.(blob2_diffs)])
        if max_diff > 10
            println("  - Variable spheres show >10% difference vs fixed radius")
            println("  - This suggests adaptive radius is capturing blob structure better")
        end
    end

    println("\n" * "="^70)
end

# Main execution
function main()
    if length(ARGS) < 1
        println("Usage: julia test_variable_spheres.jl <track_file> [nevents] [seed_radius] [step] [max_radius] [threshold]")
        println("\nExample:")
        println("  julia test_variable_spheres.jl /path/to/tracks.h5 10 3.0 1.0 10.0 0.05")
        println("\nDefaults: nevents=10, seed_radius=3.0, step=1.0, max_radius=10.0, threshold=0.05")
        exit(1)
    end

    track_file = ARGS[1]
    nevents = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 10
    seed_radius = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 3.0
    step = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0
    max_radius = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 10.0
    threshold = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 0.05

    if !isfile(track_file)
        println("Error: File not found: $track_file")
        exit(1)
    end

    analyze_sphere_expansion(track_file;
                            nevents=nevents,
                            seed_radius=seed_radius,
                            step=step,
                            max_radius=max_radius,
                            threshold=threshold)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
