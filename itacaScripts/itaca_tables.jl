#!/usr/bin/env julia

"""
Generate efficiency tables from ITACA reconstruction results.

Generates tables with:
- Rows: ε_fid (fiducial efficiency), ε_1trk (single-track efficiency)
- Columns: maxd=2mm (σt=1mm), maxd=20mm (σt=10mm)

Tables generated:
1. 0nubb (signal)
2. bi214-copper_shell
3. bi214-copper_endcaps
4. bi214-weighted_mean
5. tl208-copper_shell
6. tl208-copper_endcaps
7. tl208-weighted_mean

Directory structure expected:
    itacaScripts/
    ├── <isotope>/sigma1mm/   → contains maxd=2mm, 4mm results
    ├── <isotope>/sigma10mm/  → contains maxd=40mm results
    └── data_info/<isotope>_<zone>_info.csv

Usage:
    julia itaca_tables.jl [--output=<dir>]

Options:
    --output=X    Output directory for tables (default: tables)
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using JSON
using Printf
using CSV
using DataFrames

#=============================================================================
# Configuration
=============================================================================#

const ISOTOPES = ["0nubb", "bi214", "tl208"]
const ZONES = Dict(
    "0nubb" => [nothing],  # signal, no zone
    "bi214" => ["copper_shell", "copper_endcaps"],
    "tl208" => ["copper_shell", "copper_endcaps"]
)

# maxd values to look for: 2 in sigma1mm; 20 in sigma10mm
const MAXD_VALUES = [2, 20]

# Mapping of maxd to sigma directory
const MAXD_TO_SIGMA = Dict(
    2 => "sigma1mm",
    20 => "sigma10mm"
)

# Pattern components for each sigma (voxel = 1.5 × σt)
const SIGMA_PATTERNS = Dict(
    "sigma1mm" => "sigmat_1mm_sigmal_1mm_voxel_1.5mm",
    "sigma10mm" => "sigmat_10mm_sigmal_3mm_voxel_15mm"
)

#=============================================================================
# Helper Functions
=============================================================================#

"""
    read_info_csv(script_dir::String, isotope::String, zone::Union{String,Nothing})

Read generation info from CSV file.
Returns (num_events, saved_events) or nothing if file not found.
"""
function read_info_csv(script_dir::String, isotope::String, zone::Union{String,Nothing})
    # Build CSV filename
    if isnothing(zone)
        csv_name = "$(isotope)_info.csv"
    else
        csv_name = "$(isotope)_$(zone)_info.csv"
    end
    csv_path = joinpath(script_dir, "data_info", csv_name)

    if !isfile(csv_path)
        println("  Warning: Info CSV not found: $csv_path")
        return nothing
    end

    df = CSV.read(csv_path, DataFrame)

    if nrow(df) == 0
        println("  Warning: Info CSV is empty: $csv_path")
        return nothing
    end

    return (num_events=df.num_events[1], saved_events=df.saved_events[1])
end

"""
    discover_results_dir(script_dir::String, isotope::String, zone::Union{String,Nothing}, maxd::Int)

Find the results directory matching the given parameters.
Returns the full path or nothing if not found.
"""
function discover_results_dir(script_dir::String, isotope::String, zone::Union{String,Nothing}, maxd::Int)
    sigma_dir = MAXD_TO_SIGMA[maxd]
    sigma_pattern = SIGMA_PATTERNS[sigma_dir]

    # Build the search directory
    search_dir = joinpath(script_dir, isotope, sigma_dir)

    if !isdir(search_dir)
        println("  Warning: Sigma directory not found: $search_dir")
        return nothing
    end

    # Build the expected directory name pattern
    if isnothing(zone)
        # 0nubb: 0nubb_filtered_ecut_2400_2500_sigmat_...
        prefix = "$(isotope)_filtered_ecut_2400_2500_$(sigma_pattern)_maxd_$(maxd)mm"
    else
        # bi214/tl208: bi214_copper_shell_filtered_ecut_2400_2500_sigmat_...
        prefix = "$(isotope)_$(zone)_filtered_ecut_2400_2500_$(sigma_pattern)_maxd_$(maxd)mm"
    end

    # Look for matching directory
    for entry in readdir(search_dir)
        entry_path = joinpath(search_dir, entry)
        if isdir(entry_path) && entry == prefix
            return entry_path
        end
    end

    println("  Warning: Results directory not found for $(isotope)/$(zone) maxd=$(maxd)mm in $sigma_dir")
    return nothing
end

"""
    read_statistics_json(results_dir::String)

Read reconstruction statistics from JSON file.
Returns (n_events, n_single) or nothing if file not found.
"""
function read_statistics_json(results_dir::String)
    stats_path = joinpath(results_dir, "statistics.json")

    if !isfile(stats_path)
        println("  Warning: statistics.json not found in: $results_dir")
        return nothing
    end

    stats = JSON.parsefile(stats_path)

    counts = stats["track_multiplicity"]["counts"]
    n_events = stats["metadata"]["n_events"]

    # Get count of single-track events (key "1")
    n_single = get(counts, "1", 0)

    return (n_events=n_events, n_single=n_single)
end

"""
    compute_efficiencies_for_table(info, stats)

Compute efficiencies for table:
- ε_fid = saved_events / num_events (fiducial efficiency from generation)
- ε_1trk = n_single / n_events (single-track efficiency)
"""
function compute_efficiencies_for_table(info, stats)
    ε_fid = info.saved_events / info.num_events
    ε_1trk = stats.n_single / stats.n_events

    return (ε_fid=ε_fid, ε_1trk=ε_1trk, saved_events=info.saved_events)
end

"""
    compute_weighted_mean(results_shell, results_endcaps)

Compute weighted mean of efficiencies using saved_events as weights.
"""
function compute_weighted_mean(results_shell::Dict, results_endcaps::Dict)
    weighted_results = Dict{Int, NamedTuple}()

    for maxd in MAXD_VALUES
        if haskey(results_shell, maxd) && haskey(results_endcaps, maxd)
            r_shell = results_shell[maxd]
            r_endcaps = results_endcaps[maxd]

            w_shell = r_shell.saved_events
            w_endcaps = r_endcaps.saved_events
            w_total = w_shell + w_endcaps

            ε_fid_weighted = (r_shell.ε_fid * w_shell + r_endcaps.ε_fid * w_endcaps) / w_total
            ε_1trk_weighted = (r_shell.ε_1trk * w_shell + r_endcaps.ε_1trk * w_endcaps) / w_total

            weighted_results[maxd] = (ε_fid=ε_fid_weighted, ε_1trk=ε_1trk_weighted, saved_events=w_total)
        end
    end

    return weighted_results
end

#=============================================================================
# Table Generation Functions
=============================================================================#

"""
    format_latex_scientific(x::Float64)

Format number in LaTeX scientific notation: mantissa × 10^{exp}
"""
function format_latex_scientific(x::Float64)
    if x == 0.0
        return "0"
    end
    exp = floor(Int, log10(abs(x)))
    mantissa = x / 10.0^exp
    return @sprintf("\$%.2f \\times 10^{%d}\$", mantissa, exp)
end

"""
    format_percentage(x::Float64)

Format efficiency as percentage with 2 decimal places.
"""
function format_percentage(x::Float64)
    return @sprintf("%.2f%%", x * 100)
end

"""
    write_efficiency_table_md(results::Dict, isotope::String, zone::Union{String,Nothing}, output_dir::String)

Write efficiency table in Markdown format.
"""
function write_efficiency_table_md(results::Dict{Int, NamedTuple}, isotope::String,
                                    zone::Union{String,Nothing}, output_dir::String)
    # Build filename
    if isnothing(zone)
        filename = "efficiency_$(isotope).md"
        title = uppercase(isotope)
    elseif zone == "weighted_mean"
        filename = "efficiency_$(isotope)_weighted_mean.md"
        title = "$(uppercase(isotope)) - Weighted Mean"
    else
        filename = "efficiency_$(isotope)_$(zone).md"
        title = "$(uppercase(isotope)) - $(replace(zone, "_" => " "))"
    end

    filepath = joinpath(output_dir, filename)

    open(filepath, "w") do f
        println(f, "# Efficiency Table: $title")
        println(f)
        println(f, "| Efficiency | maxd=2mm | maxd=20mm |")
        println(f, "|------------|----------|-----------|")

        # ε_fid row
        print(f, "| ε_fid |")
        for maxd in MAXD_VALUES
            if haskey(results, maxd)
                print(f, " $(@sprintf("%.2e", results[maxd].ε_fid)) |")
            else
                print(f, " - |")
            end
        end
        println(f)

        # ε_1trk row
        print(f, "| ε_1trk |")
        for maxd in MAXD_VALUES
            if haskey(results, maxd)
                print(f, " $(@sprintf("%.2e", results[maxd].ε_1trk)) |")
            else
                print(f, " - |")
            end
        end
        println(f)
    end

    println("  Saved: $filepath")
    return filepath
end

"""
    write_efficiency_table_tex(results::Dict, isotope::String, zone::Union{String,Nothing}, output_dir::String)

Write efficiency table in LaTeX format.
"""
function write_efficiency_table_tex(results::Dict{Int, NamedTuple}, isotope::String,
                                     zone::Union{String,Nothing}, output_dir::String)
    # Build filename
    if isnothing(zone)
        filename = "efficiency_$(isotope).tex"
        title = uppercase(isotope)
    elseif zone == "weighted_mean"
        filename = "efficiency_$(isotope)_weighted_mean.tex"
        title = "$(uppercase(isotope)) - Weighted Mean"
    else
        filename = "efficiency_$(isotope)_$(zone).tex"
        title = "$(uppercase(isotope)) - $(replace(zone, "_" => " "))"
    end

    filepath = joinpath(output_dir, filename)

    open(filepath, "w") do f
        # LaTeX preamble
        println(f, "\\documentclass{standalone}")
        println(f, "\\usepackage{booktabs}")
        println(f, "\\usepackage{amsmath}")
        println(f, "\\usepackage{caption}")
        println(f, "\\begin{document}")
        println(f, "\\footnotesize")
        println(f, "\\begin{minipage}{\\textwidth}")
        println(f, "\\captionof{table}{Efficiency table for $title}")
        println(f, "\\centering")
        println(f, "\\begin{tabular}{lcc}")
        println(f, "\\toprule")
        println(f, "Efficiency & maxd=2mm & maxd=20mm \\\\")
        println(f, "\\midrule")

        # ε_fid row
        print(f, "\$\\epsilon_{\\text{fid}}\$ ")
        for maxd in MAXD_VALUES
            if haskey(results, maxd)
                print(f, "& $(format_latex_scientific(results[maxd].ε_fid)) ")
            else
                print(f, "& -- ")
            end
        end
        println(f, "\\\\")

        # ε_1trk row
        print(f, "\$\\epsilon_{1\\text{trk}}\$ ")
        for maxd in MAXD_VALUES
            if haskey(results, maxd)
                print(f, "& $(format_latex_scientific(results[maxd].ε_1trk)) ")
            else
                print(f, "& -- ")
            end
        end
        println(f, "\\\\")

        # LaTeX closing
        println(f, "\\bottomrule")
        println(f, "\\end{tabular}")
        println(f, "\\end{minipage}")
        println(f, "\\end{document}")
    end

    println("  Saved: $filepath")
    return filepath
end

#=============================================================================
# Main Processing Functions
=============================================================================#

"""
    process_isotope_zone(script_dir::String, isotope::String, zone::Union{String,Nothing})

Process a single isotope/zone combination and collect efficiencies for all maxd values.
Returns Dict{Int, NamedTuple} mapping maxd -> efficiencies.
"""
function process_isotope_zone(script_dir::String, isotope::String, zone::Union{String,Nothing})
    results = Dict{Int, NamedTuple}()

    # Read info CSV (same for all maxd values)
    info = read_info_csv(script_dir, isotope, zone)
    if isnothing(info)
        return results
    end

    for maxd in MAXD_VALUES
        # Find results directory
        results_dir = discover_results_dir(script_dir, isotope, zone, maxd)
        if isnothing(results_dir)
            continue
        end

        # Read statistics
        stats = read_statistics_json(results_dir)
        if isnothing(stats)
            continue
        end

        # Compute efficiencies
        effs = compute_efficiencies_for_table(info, stats)
        results[maxd] = effs
    end

    return results
end

"""
    main()

Main function to generate all efficiency tables.
"""
function main()
    script_dir = dirname(@__FILE__)

    # Parse arguments
    output_dir = joinpath(script_dir, "tables")

    for arg in ARGS
        if startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) == 2
                key, value = parts
                if key == "output"
                    output_dir = String(value)
                else
                    println("Warning: Unknown argument: --$key")
                end
            end
        end
    end

    # Create output directory
    if !isdir(output_dir)
        mkpath(output_dir)
        println("Created output directory: $output_dir")
    end

    println("=" ^ 70)
    println("ITACA EFFICIENCY TABLES")
    println("=" ^ 70)
    println("Output directory: $output_dir")
    println()

    # Store results for weighted mean calculations
    all_results = Dict{String, Dict{String, Dict{Int, NamedTuple}}}()

    for isotope in ISOTOPES
        all_results[isotope] = Dict{String, Dict{Int, NamedTuple}}()

        for zone in ZONES[isotope]
            zone_name = isnothing(zone) ? "signal" : zone
            println("Processing: $isotope / $zone_name")

            results = process_isotope_zone(script_dir, isotope, zone)

            if !isempty(results)
                all_results[isotope][zone_name] = results

                # Write tables
                write_efficiency_table_md(results, isotope, zone, output_dir)
                write_efficiency_table_tex(results, isotope, zone, output_dir)
            else
                println("  No results found for $isotope / $zone_name")
            end
            println()
        end

        # Compute weighted mean for bi214 and tl208
        if isotope in ["bi214", "tl208"]
            if haskey(all_results[isotope], "copper_shell") && haskey(all_results[isotope], "copper_endcaps")
                println("Computing weighted mean for $isotope...")

                weighted = compute_weighted_mean(
                    all_results[isotope]["copper_shell"],
                    all_results[isotope]["copper_endcaps"]
                )

                if !isempty(weighted)
                    write_efficiency_table_md(weighted, isotope, "weighted_mean", output_dir)
                    write_efficiency_table_tex(weighted, isotope, "weighted_mean", output_dir)
                else
                    println("  Warning: Could not compute weighted mean for $isotope")
                end
                println()
            end
        end
    end

    println("=" ^ 70)
    println("Done! Tables saved to: $output_dir")
    println("=" ^ 70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
