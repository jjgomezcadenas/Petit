#!/usr/bin/env julia

"""
Generate efficiency tables from analysis directory structure.

Scans a directory tree with structure:
    <cmdir>/
    ├── 0nubb/
    │   ├── sigma1mm/   → contains results with pitch=2mm
    │   └── sigma10mm/  → contains results with pitch=20mm
    ├── bi214/
    │   ├── sigma1mm/
    │   └── sigma10mm/
    ├── tl208/
    └── xe137/

Extracts single-track efficiency from each results directory and generates tables.

Usage:
    julia itaca_tables.jl [options]

Options:
    --cmdir=X       Base directory (default: analysis)
    --output-dir=X  Output directory (default: tables)
    --output=X      Output file prefix (default: efficiency_table)
    --dry-run       Show what would be written without writing files
    --help, -h      Show this help message

Output:
    Generates .md, .tex, and .csv files with 4 columns:
    - Source: isotope and zone
    - σ=1mm, pitch=2mm: efficiency
    - σ=10mm, pitch=20mm: efficiency
    - Ratio: σ=10mm / σ=1mm
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using JSON
using Printf
using Dates

# =============================================================================
# Configuration
# =============================================================================

const ISOTOPES = ["0nubb", "bi214", "tl208", "xe137"]
const SIGMA_DIRS = ["sigma1mm", "sigma10mm"]

# =============================================================================
# Directory Parsing
# =============================================================================

"""
    parse_results_dirname(dirname::String)

Parse a results directory name to extract isotope, zone, and maxd.
Examples:
  - "0nubb_filtered_ecut_2400_2500_sigmat_1mm_sigmal_1mm_voxel_2mm_maxd_2mm_evt_1_1000"
    → (isotope="0nubb", zone=nothing, maxd=2)
  - "bi214_copper_shell_filtered_ecut_2400_2500_sigmat_1mm_sigmal_1mm_voxel_2mm_maxd_2mm_evt_1_end"
    → (isotope="bi214", zone="copper_shell", maxd=2)
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

    # Extract zone (region) - appears between isotope and _filtered or _sigmat
    # For 0nubb: no zone
    # For bi214/tl208: zone like "copper_shell", "copper_endcaps", "ptfe_barrel"
    zone = nothing
    if isotope != "0nubb"
        # Pattern: isotope_ZONE_filtered or isotope_ZONE_sigmat (without filtered)
        # Try _filtered first, then _sigmat
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
    read_single_track_efficiency(results_dir::String)

Read statistics.json from a results directory and extract single-track efficiency.
Returns the percentage of 1-track events from track_multiplicity.percentages["1"],
converted to a fraction (0-1). Returns nothing if file not found or invalid.
"""
function read_single_track_efficiency(results_dir::String)
    json_path = joinpath(results_dir, "statistics.json")
    if !isfile(json_path)
        return nothing
    end

    try
        data = JSON.parsefile(json_path)

        # Get track_multiplicity.percentages["1"]
        track_mult = get(data, "track_multiplicity", nothing)
        if track_mult === nothing
            @warn "No track_multiplicity in $json_path"
            return nothing
        end

        percentages = get(track_mult, "percentages", nothing)
        if percentages === nothing
            @warn "No percentages in track_multiplicity in $json_path"
            return nothing
        end

        # Get the "1" key (single-track percentage)
        single_pct = get(percentages, "1", nothing)
        if single_pct === nothing
            @warn "No '1' key in percentages in $json_path"
            return nothing
        end

        # Convert percentage to fraction
        return single_pct / 100.0
    catch e
        @warn "Error reading $json_path: $e"
        return nothing
    end
end

"""
    discover_all_results(cmdir::String)

Walk the directory tree and collect all efficiency data.
Returns a vector of named tuples: (isotope, zone, sigma, maxd, efficiency, path)
"""
function discover_all_results(cmdir::String)
    results = Vector{NamedTuple{(:isotope, :zone, :sigma, :maxd, :efficiency, :path),
                                 Tuple{String, Union{String, Nothing}, String, Int, Float64, String}}}()

    for isotope in ISOTOPES
        isotope_dir = joinpath(cmdir, isotope)
        if !isdir(isotope_dir)
            @info "Skipping missing isotope directory: $isotope"
            continue
        end

        for sigma_dir in SIGMA_DIRS
            sigma_path = joinpath(isotope_dir, sigma_dir)
            if !isdir(sigma_path)
                @info "Skipping missing sigma directory: $isotope/$sigma_dir"
                continue
            end

            # Scan for results directories
            for entry in readdir(sigma_path)
                entry_path = joinpath(sigma_path, entry)
                if !isdir(entry_path)
                    continue
                end

                # Parse directory name
                parsed = parse_results_dirname(entry)
                if parsed === nothing
                    @warn "Could not parse directory name: $entry"
                    continue
                end

                # Read efficiency
                efficiency = read_single_track_efficiency(entry_path)
                if efficiency === nothing
                    @warn "Could not read efficiency from: $entry_path"
                    continue
                end

                push!(results, (
                    isotope = parsed.isotope,
                    zone = parsed.zone,
                    sigma = sigma_dir,
                    maxd = parsed.maxd,
                    efficiency = efficiency,
                    path = entry_path
                ))

                println("  Found: $(parsed.isotope) $(parsed.zone === nothing ? "" : "($(parsed.zone))") " *
                        "$sigma_dir maxd=$(parsed.maxd)mm → ε_1trk = $(round(efficiency, digits=4))")
            end
        end
    end

    return results
end

# =============================================================================
# Table Generation
# =============================================================================

"""
    format_row_label(isotope::String, zone::Union{String, Nothing})

Format a row label for the table.
"""
function format_row_label(isotope::String, zone::Union{String, Nothing})
    if zone === nothing
        return isotope
    else
        return "$isotope ($zone)"
    end
end

"""
    format_efficiency(eff::Union{Float64, Nothing})

Format efficiency value for display.
"""
function format_efficiency(eff::Union{Float64, Nothing})
    if eff === nothing
        return "-"
    elseif eff >= 0.01
        return @sprintf("%.4f", eff)
    else
        # Scientific notation for small values
        exp = floor(Int, log10(abs(eff)))
        mantissa = eff / 10.0^exp
        return @sprintf("%.2f×10^%d", mantissa, exp)
    end
end

"""
    format_efficiency_latex(eff::Union{Float64, Nothing})

Format efficiency value in LaTeX.
"""
function format_efficiency_latex(eff::Union{Float64, Nothing})
    if eff === nothing
        return "-"
    elseif eff >= 0.01
        return @sprintf("%.4f", eff)
    else
        exp = floor(Int, log10(abs(eff)))
        mantissa = eff / 10.0^exp
        return @sprintf("\$%.2f \\times 10^{%d}\$", mantissa, exp)
    end
end

"""
    format_ratio(ratio::Union{Float64, Nothing})

Format ratio value for display.
"""
function format_ratio(ratio::Union{Float64, Nothing})
    if ratio === nothing
        return "-"
    else
        return @sprintf("%.3f", ratio)
    end
end

"""
    build_table_data(results::Vector)

Organize results into a simplified 4-column table structure.
Returns a vector of named tuples: (label, sigma_1mm, sigma_10mm, ratio)
where:
- label: row label (isotope + zone)
- sigma_1mm: efficiency for σ=1mm, pitch=2mm
- sigma_10mm: efficiency for σ=10mm, pitch=20mm
- ratio: sigma_10mm / sigma_1mm
"""
function build_table_data(results::Vector)
    # Collect unique rows (isotope, zone) - maintain order
    row_set = Set{Tuple{String, Union{String, Nothing}}}()
    for r in results
        push!(row_set, (r.isotope, r.zone))
    end

    # Sort rows: by isotope order, then by zone
    rows = sort(collect(row_set), by = x -> (findfirst(==(x[1]), ISOTOPES), something(x[2], "")))

    # Build table data
    table_rows = Vector{NamedTuple{(:label, :sigma_1mm, :sigma_10mm, :ratio),
                                    Tuple{String, Union{Float64, Nothing}, Union{Float64, Nothing}, Union{Float64, Nothing}}}}()

    for (isotope, zone) in rows
        label = format_row_label(isotope, zone)

        # Find efficiency for sigma1mm (expect pitch=2mm)
        sigma_1mm_eff = nothing
        for r in results
            if r.isotope == isotope && r.zone == zone && r.sigma == "sigma1mm"
                sigma_1mm_eff = r.efficiency
                break
            end
        end

        # Find efficiency for sigma10mm (expect pitch=20mm)
        sigma_10mm_eff = nothing
        for r in results
            if r.isotope == isotope && r.zone == zone && r.sigma == "sigma10mm"
                sigma_10mm_eff = r.efficiency
                break
            end
        end

        # Compute ratio
        ratio = nothing
        if sigma_1mm_eff !== nothing && sigma_10mm_eff !== nothing && sigma_1mm_eff > 0
            ratio = sigma_10mm_eff / sigma_1mm_eff
        end

        push!(table_rows, (label=label, sigma_1mm=sigma_1mm_eff, sigma_10mm=sigma_10mm_eff, ratio=ratio))
    end

    return table_rows
end

"""
    write_markdown_table(table_data, output_path::String)

Write the efficiency table in Markdown format.
"""
function write_markdown_table(table_data, output_path::String)
    open(output_path, "w") do f
        println(f, "# Single-Track Efficiency Table")
        println(f)
        println(f, "Generated: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))")
        println(f)

        # Header row
        println(f, "| Source | σ=1mm, pitch=2mm | σ=10mm, pitch=20mm | Ratio |")
        println(f, "|--------|------------------|--------------------| ------|")

        # Data rows
        for row in table_data
            eff_1mm = format_efficiency(row.sigma_1mm)
            eff_10mm = format_efficiency(row.sigma_10mm)
            ratio_str = format_ratio(row.ratio)
            println(f, "| $(row.label) | $eff_1mm | $eff_10mm | $ratio_str |")
        end
    end

    println("Saved Markdown table: $output_path")
end

"""
    write_latex_table(table_data, output_path::String)

Write the efficiency table in LaTeX format.
"""
function write_latex_table(table_data, output_path::String)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    open(output_path, "w") do f
        # Document header
        println(f, "% Generated: $timestamp")
        println(f, "\\documentclass{article}")
        println(f, "\\usepackage{booktabs}")
        println(f, "\\usepackage{siunitx}")
        println(f, "\\begin{document}")
        println(f)
        println(f, "\\begin{table}[htbp]")
        println(f, "\\centering")
        println(f, "\\caption{Single-Track Efficiency (Generated: $timestamp)}")

        # Table structure: 4 columns
        println(f, "\\begin{tabular}{lccc}")
        println(f, "\\toprule")

        # Header row
        println(f, "Source & \$\\sigma=1\\,\\text{mm}\$, pitch=2mm & \$\\sigma=10\\,\\text{mm}\$, pitch=20mm & Ratio \\\\")
        println(f, "\\midrule")

        # Data rows
        for row in table_data
            # Format isotope with superscript for LaTeX
            label = row.label
            # Replace isotope names with LaTeX formatting
            label = replace(label, "0nubb" => "\$0\\nu\\beta\\beta\$")
            label = replace(label, "bi214" => "\$^{214}\$Bi")
            label = replace(label, "tl208" => "\$^{208}\$Tl")
            label = replace(label, "xe137" => "\$^{137}\$Xe")
            # Escape underscores
            label = replace(label, "_" => "\\_")

            eff_1mm = format_efficiency_latex(row.sigma_1mm)
            eff_10mm = format_efficiency_latex(row.sigma_10mm)
            ratio_str = format_ratio(row.ratio)

            println(f, "$label & $eff_1mm & $eff_10mm & $ratio_str \\\\")
        end

        println(f, "\\bottomrule")
        println(f, "\\end{tabular}")
        println(f, "\\end{table}")
        println(f)
        println(f, "\\end{document}")
    end

    println("Saved LaTeX table: $output_path")
end

"""
    write_csv_table(table_data, output_path::String)

Write the efficiency table in CSV format.
"""
function write_csv_table(table_data, output_path::String)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    open(output_path, "w") do f
        # Timestamp and header row
        println(f, "# Generated: $timestamp")
        println(f, "source,sigma_1mm,sigma_10mm,ratio")

        # Data rows
        for row in table_data
            eff_1mm = row.sigma_1mm === nothing ? "" : @sprintf("%.6f", row.sigma_1mm)
            eff_10mm = row.sigma_10mm === nothing ? "" : @sprintf("%.6f", row.sigma_10mm)
            ratio_str = row.ratio === nothing ? "" : @sprintf("%.6f", row.ratio)
            println(f, "$(row.label),$eff_1mm,$eff_10mm,$ratio_str")
        end
    end

    println("Saved CSV table: $output_path")
end

# =============================================================================
# Main
# =============================================================================

function main()
    # Default values
    cmdir = "analysis"
    output_dir = "tables"
    output_prefix = "efficiency_table"
    dry_run = false

    # Parse arguments
    for arg in ARGS
        if startswith(arg, "--cmdir=")
            cmdir = String(split(arg, '=')[2])
        elseif startswith(arg, "--output-dir=")
            output_dir = String(split(arg, '=')[2])
        elseif startswith(arg, "--output=")
            output_prefix = String(split(arg, '=')[2])
        elseif arg == "--dry-run"
            dry_run = true
        elseif arg == "--help" || arg == "-h"
            println("""
Usage: julia itaca_tables.jl [options]

Options:
    --cmdir=X       Base directory (default: analysis)
    --output-dir=X  Output directory (default: tables)
    --output=X      Output file prefix (default: efficiency_table)
    --dry-run       Show what would be written without writing files
    --help, -h      Show this help message

Directory structure expected:
    <cmdir>/
    ├── 0nubb/sigma1mm/..., sigma10mm/...
    ├── bi214/sigma1mm/..., sigma10mm/...
    ├── tl208/sigma1mm/..., sigma10mm/...
    └── xe137/sigma1mm/..., sigma10mm/...

Output:
    Generates three files in <output-dir>/:
    - <prefix>.md   (Markdown format)
    - <prefix>.tex  (LaTeX format)
    - <prefix>.csv  (CSV format)

    Table columns:
    - Source: isotope and zone
    - σ=1mm, pitch=2mm: efficiency
    - σ=10mm, pitch=20mm: efficiency
    - Ratio: σ=10mm / σ=1mm
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

    println("=" ^ 60)
    println("ITACA Efficiency Table Generator")
    println("=" ^ 60)
    println("Base directory: $cmdir")
    println("Output directory: $output_dir")
    println("Output prefix: $output_prefix")
    println("Dry run: $dry_run")
    println()

    # Discover all results
    println("Scanning directory tree...")
    results = discover_all_results(cmdir)

    if isempty(results)
        println("\nNo results found!")
        return
    end

    println("\nFound $(length(results)) result directories")

    # Build table structure
    table_data = build_table_data(results)

    println("\nTable structure:")
    println("  Rows: $(length(table_data)) (isotope/zone combinations)")
    println("  Columns: Source, σ=1mm, σ=10mm, Ratio")

    # Preview table content
    println("\nTable preview:")
    println("-" ^ 60)
    @printf("%-30s %10s %10s %8s\n", "Source", "σ=1mm", "σ=10mm", "Ratio")
    println("-" ^ 60)
    for row in table_data
        eff_1mm = row.sigma_1mm === nothing ? "-" : @sprintf("%.4f", row.sigma_1mm)
        eff_10mm = row.sigma_10mm === nothing ? "-" : @sprintf("%.4f", row.sigma_10mm)
        ratio_str = row.ratio === nothing ? "-" : @sprintf("%.3f", row.ratio)
        @printf("%-30s %10s %10s %8s\n", row.label, eff_1mm, eff_10mm, ratio_str)
    end
    println("-" ^ 60)

    # Write output files
    println()
    if dry_run
        println("DRY RUN - would write the following files:")
        println("  - $(joinpath(output_dir, output_prefix)).md")
        println("  - $(joinpath(output_dir, output_prefix)).tex")
        println("  - $(joinpath(output_dir, output_prefix)).csv")
    else
        # Create output directory if needed
        script_dir = dirname(@__FILE__)
        full_output_dir = joinpath(script_dir, output_dir)
        if !isdir(full_output_dir)
            mkpath(full_output_dir)
            println("Created output directory: $full_output_dir")
        end

        # Write all three formats
        md_path = joinpath(full_output_dir, "$(output_prefix).md")
        tex_path = joinpath(full_output_dir, "$(output_prefix).tex")
        csv_path = joinpath(full_output_dir, "$(output_prefix).csv")

        write_markdown_table(table_data, md_path)
        write_latex_table(table_data, tex_path)
        write_csv_table(table_data, csv_path)
    end

    println()
    println("=" ^ 60)
    println("Done!")
    println("=" ^ 60)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
