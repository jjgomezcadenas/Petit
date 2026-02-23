#!/usr/bin/env julia

"""
Generate event statistics plots from analysis directory structure.

Scans a directory tree with structure:
    <cmdir>/
    ├── 0nubb/
    │   ├── sigma1mm/   → contains results with event_stats.csv
    │   └── sigma10mm/
    ├── bi214/
    │   ├── sigma1mm/
    │   └── sigma10mm/
    └── tl208/

Reads event_stats.csv from each results directory and creates 3x2 histogram plots.

Usage:
    julia itaca_plot_stats.jl [options]

Options:
    --cmdir=X       Base directory (default: analysis)
    --output-dir=X  Output directory for plots (default: statsplots)
    --dry-run       Show what would be done without generating plots
    --help, -h      Show this help message

Output:
    One PNG file per results directory, named after the directory.
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using CSV
using DataFrames
using Statistics
using Printf
using Plots

# =============================================================================
# Configuration
# =============================================================================

const ISOTOPES = ["0nubb", "bi214", "tl208", "xe137"]
const SIGMA_DIRS = ["sigma1mm", "sigma10mm"]

# Energy ROI limits for main track energy plot
const ENERGY_MIN = 2400.0  # keV
const ENERGY_MAX = 2500.0  # keV

# Variables to plot (6 total, 3x2 layout)
# Note: track3plus_energy_keV replaced by main_track_energy_keV in ROI
const PLOT_VARIABLES = [
    (col=:n_tracks, label="N Tracks", xlabel="Number of Tracks", xlim=nothing),
    (col=:total_energy_keV, label="Total Energy", xlabel="Energy (keV)", xlim=nothing),
    (col=:main_track_energy_keV, label="Main Track Energy", xlabel="Energy (keV)", xlim=nothing),
    (col=:track2_energy_keV, label="Track 2 Energy", xlabel="Energy (keV)", xlim=(0, 300)),
    (col=:main_track_energy_keV, label="Main Track Energy (ROI)", xlabel="Energy (keV)", xlim=(ENERGY_MIN, ENERGY_MAX)),
    (col=:track1_track2_distance_mm, label="Track 1-2 Distance", xlabel="Distance (mm)", xlim=(0, 500))
]

# =============================================================================
# Directory Parsing (adapted from itaca_tables.jl)
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
Returns a vector of named tuples: (isotope, zone, sigma, dirname, path)
"""
function discover_all_results_dirs(cmdir::String)
    results = Vector{NamedTuple{(:isotope, :zone, :sigma, :dirname, :path),
                                 Tuple{String, Union{String, Nothing}, String, String, String}}}()

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
                    dirname = entry,
                    path = entry_path
                ))
            end
        end
    end

    return results
end

# =============================================================================
# Plotting Functions
# =============================================================================

"""
    add_stats_legend!(p, data; digits=1)

Add μ and σ as legend entries using invisible scatter plots.
"""
function add_stats_legend!(p, data; digits=1)
    valid_data = filter(!isnan, data)
    if isempty(valid_data)
        return p
    end
    μ = mean(valid_data)
    σ = std(valid_data)
    μ_str = round(μ, digits=digits)
    σ_str = round(σ, digits=digits)
    plot!(p, [], [], label="μ = $μ_str", seriestype=:scatter, markersize=0)
    plot!(p, [], [], label="σ = $σ_str", seriestype=:scatter, markersize=0)
    return p
end

"""
    create_histogram(data, var_info; nbins=30)

Create a single histogram subplot for a variable.
If var_info.xlim is set, filter data to that range and set axis limits.
"""
function create_histogram(data, var_info; nbins=30)
    valid_data = filter(!isnan, data)

    # Apply xlim filter if specified
    if var_info.xlim !== nothing
        xmin, xmax = var_info.xlim
        valid_data = filter(x -> xmin <= x <= xmax, valid_data)
    end

    if isempty(valid_data)
        p = plot(title="$(var_info.label)\n(no data)", titlefontsize=9, legend=false)
        return p
    end

    # Build histogram with optional xlims
    if var_info.xlim !== nothing
        p = histogram(valid_data, bins=nbins,
                      xlabel=var_info.xlabel, ylabel="Events",
                      title=var_info.label, titlefontsize=9,
                      label="", color=:steelblue, alpha=0.7,
                      bottom_margin=3Plots.mm, left_margin=2Plots.mm,
                      xlims=var_info.xlim)
    else
        p = histogram(valid_data, bins=nbins,
                      xlabel=var_info.xlabel, ylabel="Events",
                      title=var_info.label, titlefontsize=9,
                      label="", color=:steelblue, alpha=0.7,
                      bottom_margin=3Plots.mm, left_margin=2Plots.mm)
    end

    # Add statistics
    add_stats_legend!(p, valid_data; digits=1)

    return p
end

"""
    create_event_stats_plot(df::DataFrame, title_str::String)

Create a 3x2 plot with histograms for the 6 variables.
"""
function create_event_stats_plot(df::DataFrame, title_str::String)
    n_events = nrow(df)

    # Create individual subplots
    subplots = []
    for var_info in PLOT_VARIABLES
        if hasproperty(df, var_info.col)
            data = df[!, var_info.col]
            p = create_histogram(data, var_info)
        else
            p = plot(title="$(var_info.label)\n(column not found)", titlefontsize=9, legend=false)
        end
        push!(subplots, p)
    end

    # Combine into 3x2 layout
    plt = plot(subplots..., layout=(3, 2), size=(900, 900),
               plot_title="Event Statistics (N=$n_events)\n$title_str",
               plot_titlefontsize=10)

    return plt
end

"""
    format_title(result)

Format a descriptive title from the result info.
"""
function format_title(result)
    isotope = result.isotope
    zone = result.zone === nothing ? "" : " ($(result.zone))"
    sigma = result.sigma == "sigma1mm" ? "σ=1mm" : "σ=10mm"
    return "$isotope$zone - $sigma"
end

# =============================================================================
# Main
# =============================================================================

function main()
    # Default values
    script_dir = dirname(@__FILE__)
    cmdir = joinpath(script_dir, "analysis")
    output_dir = joinpath(script_dir, "statsplots")
    dry_run = false

    # Parse arguments
    for arg in ARGS
        if startswith(arg, "--cmdir=")
            cmdir = String(split(arg, '=')[2])
        elseif startswith(arg, "--output-dir=")
            output_dir = String(split(arg, '=')[2])
        elseif arg == "--dry-run"
            dry_run = true
        elseif arg == "--help" || arg == "-h"
            println("""
Usage: julia itaca_plot_stats.jl [options]

Options:
    --cmdir=X       Base directory (default: analysis)
    --output-dir=X  Output directory for plots (default: statsplots)
    --dry-run       Show what would be done without generating plots
    --help, -h      Show this help message

Directory structure expected:
    <cmdir>/
    ├── 0nubb/sigma1mm/..., sigma10mm/...
    ├── bi214/sigma1mm/..., sigma10mm/...
    ├── tl208/sigma1mm/..., sigma10mm/...
    └── xe137/sigma1mm/..., sigma10mm/...

Output:
    PNG files in <output-dir>/ named after each results directory.
    Each plot contains 6 histograms in 3x2 layout:
    - n_tracks, total_energy_keV, main_track_energy_keV
    - track2_energy_keV (0-300 keV), main_track_energy_keV (ROI 2400-2500 keV), track1_track2_distance_mm (0-500 mm)
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
    println("ITACA Event Statistics Plotter")
    println("=" ^ 70)
    println("Base directory:   $cmdir")
    println("Output directory: $output_dir")
    println("Dry run:          $dry_run")
    println()

    # Discover all results directories
    println("Scanning directory tree...")
    results = discover_all_results_dirs(cmdir)

    if isempty(results)
        println("\nNo results directories found with event_stats.csv!")
        return
    end

    println("\nFound $(length(results)) directories with event_stats.csv")
    println()

    # Create output directory if needed
    if !dry_run && !isdir(output_dir)
        mkpath(output_dir)
        println("Created output directory: $output_dir")
    end

    # Process each directory
    n_success = 0
    n_failed = 0

    for (i, result) in enumerate(results)
        println("[$i/$(length(results))] Processing: $(result.dirname)")

        csv_path = joinpath(result.path, "event_stats.csv")
        output_path = joinpath(output_dir, "$(result.dirname).png")

        if dry_run
            println("  Would create: $output_path")
            continue
        end

        try
            # Read CSV
            df = CSV.read(csv_path, DataFrame)

            if nrow(df) == 0
                println("  Warning: Empty CSV file, skipping")
                n_failed += 1
                continue
            end

            # Create plot
            title_str = format_title(result)
            plt = create_event_stats_plot(df, title_str)

            # Save plot
            savefig(plt, output_path)
            println("  Saved: $output_path")
            n_success += 1

        catch e
            println("  Error: $e")
            n_failed += 1
        end
    end

    # Summary
    println()
    println("=" ^ 70)
    if dry_run
        println("DRY RUN complete - would create $(length(results)) plots")
    else
        println("Complete: $n_success plots created, $n_failed failed")
    end
    println("=" ^ 70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
