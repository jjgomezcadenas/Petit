#!/usr/bin/env julia

"""
ITACA Event Statistics - Sequential Cut Analysis.

Reads event_stats.csv files from the analysis directory structure and applies
sequential cuts to compute efficiencies:

1. Total events: Number of events in the CSV (before any cuts)
2. Single-track cut: Events with n_tracks == 1
3. Energy ROI cut: Events with main_track_energy_keV in [2400, 2500] keV

For each isotope/zone/sigma/maxd combination, outputs:
- Event counts at each cut stage
- Cut efficiencies (relative and cumulative)
- Final filtered DataFrame saved as CSV

Directory structure expected:
    analysis/
    ├── 0nubb/
    │   ├── sigma1mm/   → contains results with event_stats.csv
    │   └── sigma10mm/
    ├── bi214/
    │   ├── sigma1mm/
    │   └── sigma10mm/
    └── tl208/

Usage:
    julia itaca_event_stats.jl [options]

Options:
    --cmdir=X       Base directory (default: analysis)
    --output-dir=X  Output directory for filtered CSVs (default: event_stats_filtered)
    --energy-min=X  Minimum energy for ROI cut (default: 2400 keV)
    --energy-max=X  Maximum energy for ROI cut (default: 2500 keV)
    --help, -h      Show this help message

Output:
    - Console summary table with event counts and efficiencies
    - CSV files with events passing all cuts
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using CSV
using DataFrames
using Printf

# =============================================================================
# Configuration
# =============================================================================

const ISOTOPES = ["0nubb", "bi214", "tl208", "xe137"]
const SIGMA_DIRS = ["sigma1mm", "sigma10mm"]

# Default energy ROI
const DEFAULT_ENERGY_MIN = 2400.0  # keV
const DEFAULT_ENERGY_MAX = 2500.0  # keV

# =============================================================================
# Directory Parsing (same as itaca_plot_stats.jl)
# =============================================================================

"""
    parse_results_dirname(dirname::String)

Parse a results directory name to extract isotope, zone, and maxd.
"""
function parse_results_dirname(dirname::String)
    # Extract maxd value: look for _maxd_Xmm_ pattern
    maxd_match = match(r"_maxd_(\d+)mm", dirname)
    if maxd_match === nothing
        return nothing
    end
    maxd = parse(Int, maxd_match.captures[1])

    # Determine isotope from the start of the name
    isotope = nothing
    for iso in ISOTOPES
        if startswith(dirname, iso)
            isotope = iso
            break
        end
    end
    if isotope === nothing
        return nothing
    end

    # Extract zone (region)
    zone = nothing
    if isotope != "0nubb"
        zone_pattern1 = Regex("^$(isotope)_(.+?)_filtered")
        zone_pattern2 = Regex("^$(isotope)_(.+?)_sigmat")

        zone_match = match(zone_pattern1, dirname)
        if zone_match === nothing
            zone_match = match(zone_pattern2, dirname)
        end

        if zone_match !== nothing
            zone = zone_match.captures[1]
        end
    end

    return (isotope=isotope, zone=zone, maxd=maxd)
end

"""
    discover_all_results_dirs(cmdir::String)

Walk the directory tree and collect all results directories with event_stats.csv.
Returns a vector of named tuples: (isotope, zone, sigma, maxd, dirname, path)
"""
function discover_all_results_dirs(cmdir::String)
    results = Vector{NamedTuple{(:isotope, :zone, :sigma, :maxd, :dirname, :path),
                                 Tuple{String, Union{String, Nothing}, String, Int, String, String}}}()

    for isotope in ISOTOPES
        isotope_dir = joinpath(cmdir, isotope)
        if !isdir(isotope_dir)
            continue
        end

        for sigma_dir in SIGMA_DIRS
            sigma_path = joinpath(isotope_dir, sigma_dir)
            if !isdir(sigma_path)
                continue
            end

            # Scan for results directories
            for entry in readdir(sigma_path)
                entry_path = joinpath(sigma_path, entry)
                if !isdir(entry_path)
                    continue
                end

                # Check if event_stats.csv exists
                csv_path = joinpath(entry_path, "event_stats.csv")
                if !isfile(csv_path)
                    continue
                end

                # Parse directory name
                parsed = parse_results_dirname(entry)
                if parsed === nothing
                    @warn "Could not parse directory name: $entry"
                    continue
                end

                push!(results, (
                    isotope = parsed.isotope,
                    zone = parsed.zone,
                    sigma = sigma_dir,
                    maxd = parsed.maxd,
                    dirname = entry,
                    path = entry_path
                ))
            end
        end
    end

    return results
end

# =============================================================================
# Cut Analysis Functions
# =============================================================================

"""
    CutResult

Holds the result of applying sequential cuts to a DataFrame.
"""
struct CutResult
    isotope::String
    zone::Union{String, Nothing}
    sigma::String
    maxd::Int
    dirname::String
    n_total::Int          # Events before any cut
    n_1track::Int         # Events after n_tracks == 1 cut
    n_roi::Int            # Events after energy ROI cut
    eff_1track::Float64   # n_1track / n_total
    eff_roi::Float64      # n_roi / n_1track (relative)
    eff_total::Float64    # n_roi / n_total (cumulative)
    df_final::DataFrame   # Final filtered DataFrame
end

"""
    apply_cuts(df::DataFrame, energy_min::Float64, energy_max::Float64)

Apply sequential cuts to the DataFrame:
1. n_tracks == 1
2. main_track_energy_keV in [energy_min, energy_max]

Returns (n_total, n_1track, n_roi, df_final)
"""
function apply_cuts(df::DataFrame, energy_min::Float64, energy_max::Float64)
    # Count total events
    n_total = nrow(df)

    # Cut 1: Single track
    df_1track = filter(row -> row.n_tracks == 1, df)
    n_1track = nrow(df_1track)

    # Cut 2: Energy ROI on main_track_energy_keV
    df_roi = filter(row -> energy_min <= row.main_track_energy_keV <= energy_max, df_1track)
    n_roi = nrow(df_roi)

    return (n_total, n_1track, n_roi, df_roi)
end

"""
    analyze_directory(result, energy_min::Float64, energy_max::Float64)

Analyze a single results directory: read CSV, apply cuts, compute efficiencies.
"""
function analyze_directory(result, energy_min::Float64, energy_max::Float64)
    csv_path = joinpath(result.path, "event_stats.csv")

    # Read CSV
    df = CSV.read(csv_path, DataFrame)

    # Apply cuts
    n_total, n_1track, n_roi, df_final = apply_cuts(df, energy_min, energy_max)

    # Compute efficiencies
    eff_1track = n_total > 0 ? n_1track / n_total : 0.0
    eff_roi = n_1track > 0 ? n_roi / n_1track : 0.0
    eff_total = n_total > 0 ? n_roi / n_total : 0.0

    return CutResult(
        result.isotope,
        result.zone,
        result.sigma,
        result.maxd,
        result.dirname,
        n_total,
        n_1track,
        n_roi,
        eff_1track,
        eff_roi,
        eff_total,
        df_final
    )
end

# =============================================================================
# Output Functions
# =============================================================================

"""
    format_zone(zone::Union{String, Nothing})

Format zone for display (handle nothing case).
"""
function format_zone(zone::Union{String, Nothing})
    return zone === nothing ? "-" : zone
end

"""
    format_sigma(sigma::String)

Format sigma directory name for display.
"""
function format_sigma(sigma::String)
    return sigma == "sigma1mm" ? "1mm" : "10mm"
end

"""
    print_results_table(results::Vector{CutResult}, energy_min::Float64, energy_max::Float64)

Print summary table of cut results.
"""
function print_results_table(results::Vector{CutResult}, energy_min::Float64, energy_max::Float64)
    println()
    println("=" ^ 120)
    println("ITACA Event Statistics - Sequential Cut Analysis")
    println("Energy ROI: [$(Int(energy_min)), $(Int(energy_max))] keV")
    println("=" ^ 120)
    println()

    # Header
    @printf("%-8s  %-15s  %-6s  %5s  %8s  %8s  %8s  %10s  %10s  %10s\n",
            "Isotope", "Zone", "Sigma", "maxd", "N_total", "N_1track", "N_ROI",
            "ε_1track", "ε_ROI|1t", "ε_total")
    println("-" ^ 120)

    # Group by isotope for cleaner output
    for isotope in ISOTOPES
        iso_results = filter(r -> r.isotope == isotope, results)
        if isempty(iso_results)
            continue
        end

        # Sort by zone, sigma, maxd
        sorted = sort(iso_results, by=r -> (format_zone(r.zone), r.sigma, r.maxd))

        for r in sorted
            @printf("%-8s  %-15s  %-6s  %5d  %8d  %8d  %8d  %10.4f  %10.4f  %10.6f\n",
                    r.isotope, format_zone(r.zone), format_sigma(r.sigma), r.maxd,
                    r.n_total, r.n_1track, r.n_roi,
                    r.eff_1track, r.eff_roi, r.eff_total)
        end
        println()  # Blank line between isotopes
    end

    println("=" ^ 120)
    println()
    println("Legend:")
    println("  N_total:   Total events in event_stats.csv")
    println("  N_1track:  Events with n_tracks == 1")
    println("  N_ROI:     Events in energy ROI after single-track cut")
    println("  ε_1track:  Single-track efficiency (N_1track / N_total)")
    println("  ε_ROI|1t:  ROI efficiency given single track (N_ROI / N_1track)")
    println("  ε_total:   Total efficiency (N_ROI / N_total)")
end

"""
    write_filtered_csv(result::CutResult, output_dir::String)

Write the filtered DataFrame to a CSV file.
"""
function write_filtered_csv(result::CutResult, output_dir::String)
    # Create filename
    filename = "$(result.dirname)_filtered.csv"
    output_path = joinpath(output_dir, filename)

    # Write CSV
    CSV.write(output_path, result.df_final)

    return output_path
end

"""
    write_summary_csv(results::Vector{CutResult}, output_path::String, energy_min::Float64, energy_max::Float64)

Write a summary CSV with all cut results.
"""
function write_summary_csv(results::Vector{CutResult}, output_path::String, energy_min::Float64, energy_max::Float64)
    rows = []

    for r in results
        push!(rows, (
            isotope = r.isotope,
            zone = format_zone(r.zone),
            sigma = format_sigma(r.sigma),
            maxd_mm = r.maxd,
            n_total = r.n_total,
            n_1track = r.n_1track,
            n_roi = r.n_roi,
            eff_1track = r.eff_1track,
            eff_roi_given_1track = r.eff_roi,
            eff_total = r.eff_total,
            energy_min_keV = energy_min,
            energy_max_keV = energy_max
        ))
    end

    df = DataFrame(rows)
    CSV.write(output_path, df)

    return output_path
end

# =============================================================================
# Main
# =============================================================================

function main()
    # Default values
    script_dir = dirname(@__FILE__)
    cmdir = joinpath(script_dir, "analysis")
    output_dir = joinpath(script_dir, "event_stats_filtered")
    energy_min = DEFAULT_ENERGY_MIN
    energy_max = DEFAULT_ENERGY_MAX

    # Parse arguments
    for arg in ARGS
        if startswith(arg, "--cmdir=")
            cmdir = String(split(arg, '=')[2])
        elseif startswith(arg, "--output-dir=")
            output_dir = String(split(arg, '=')[2])
        elseif startswith(arg, "--energy-min=")
            energy_min = parse(Float64, split(arg, '=')[2])
        elseif startswith(arg, "--energy-max=")
            energy_max = parse(Float64, split(arg, '=')[2])
        elseif arg == "--help" || arg == "-h"
            println("""
Usage: julia itaca_event_stats.jl [options]

Options:
    --cmdir=X       Base directory (default: analysis)
    --output-dir=X  Output directory for filtered CSVs (default: event_stats_filtered)
    --energy-min=X  Minimum energy for ROI cut (default: 2400 keV)
    --energy-max=X  Maximum energy for ROI cut (default: 2500 keV)
    --help, -h      Show this help message

Directory structure expected:
    <cmdir>/
    ├── 0nubb/sigma1mm/..., sigma10mm/...
    ├── bi214/sigma1mm/..., sigma10mm/...
    ├── tl208/sigma1mm/..., sigma10mm/...
    └── xe137/sigma1mm/..., sigma10mm/...

Output:
    Filtered CSV files and summary table in <output-dir>/

Sequential cuts applied:
    1. n_tracks == 1 (single-track selection)
    2. main_track_energy_keV in [energy_min, energy_max] (ROI selection)
""")
            return
        else
            @warn "Unknown argument: $arg"
        end
    end

    # Validate cmdir
    if !isdir(cmdir)
        error("Directory not found: $cmdir")
    end

    println("=" ^ 70)
    println("ITACA Event Statistics - Sequential Cut Analysis")
    println("=" ^ 70)
    println("Base directory:   $cmdir")
    println("Output directory: $output_dir")
    println("Energy ROI:       [$energy_min, $energy_max] keV")
    println()

    # Discover all results directories
    println("Scanning directory tree...")
    dirs = discover_all_results_dirs(cmdir)

    if isempty(dirs)
        println("\nNo results directories found with event_stats.csv!")
        return
    end

    println("Found $(length(dirs)) directories with event_stats.csv")
    println()

    # Create output directory
    if !isdir(output_dir)
        mkpath(output_dir)
        println("Created output directory: $output_dir")
    end

    # Process each directory
    results = CutResult[]
    n_success = 0
    n_failed = 0

    for (i, dir_info) in enumerate(dirs)
        print("[$i/$(length(dirs))] Processing: $(dir_info.dirname)... ")

        try
            result = analyze_directory(dir_info, energy_min, energy_max)
            push!(results, result)

            # Write filtered CSV
            if nrow(result.df_final) > 0
                output_path = write_filtered_csv(result, output_dir)
                println("OK ($(result.n_roi) events)")
            else
                println("OK (0 events, no CSV written)")
            end

            n_success += 1
        catch e
            println("ERROR: $e")
            n_failed += 1
        end
    end

    # Print results table
    if !isempty(results)
        print_results_table(results, energy_min, energy_max)

        # Write summary CSV
        summary_path = joinpath(output_dir, "cut_summary.csv")
        write_summary_csv(results, summary_path, energy_min, energy_max)
        println("Summary CSV: $summary_path")
    end

    # Final summary
    println()
    println("=" ^ 70)
    println("Complete: $n_success processed, $n_failed failed")
    println("=" ^ 70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
