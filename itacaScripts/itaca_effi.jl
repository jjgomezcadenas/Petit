#!/usr/bin/env julia

"""
Compute efficiencies from ITACA reconstruction results.

Computes three efficiencies:
1. interaction_efficiency = saved_events / num_events (from generation)
2. fiducial_efficiency = n_events (in ROI) / saved_events
3. eff_1trk = n_single_track / n_events

Reads data from:
- <isotope>_<zone>_info.csv (generation info)
- statistics.json (reconstruction results)

Usage:
    julia itaca_effi.jl --isotope=<isotope> --zone=<zone> [options]

Options:
    --isotope=X     Isotope: bi214, tl208, or all (required)
    --zone=X        Zone: copper_shell, copper_endcaps, etc. (required)
    --voxel=X       Voxel size in mm (default: 1)
    --maxd=X        Max distance in mm (default: 4)
    --writetable=X  Write efficiency table to files (default: false)

Output:
    Three efficiencies in scientific notation (2 decimals)
    If --writetable=true, creates directory with .md and .tex table files
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

"""
    read_info_csv(script_dir::String, isotope::String, zone::String)

Read generation info from CSV file.
Returns (num_events, saved_events) or nothing if file not found.
"""
function read_info_csv(script_dir::String, isotope::String, zone::String)
    csv_path = joinpath(script_dir, "data_info", "$(isotope)_$(zone)_info.csv")

    if !isfile(csv_path)
        println("Warning: Info CSV not found: $csv_path")
        return nothing
    end

    df = CSV.read(csv_path, DataFrame)

    if nrow(df) == 0
        println("Warning: Info CSV is empty: $csv_path")
        return nothing
    end

    return (num_events=df.num_events[1], saved_events=df.saved_events[1])
end

"""
    read_statistics_json(results_dir::String)

Read reconstruction statistics from JSON file.
Returns (n_events, n_single) or nothing if file not found.
"""
function read_statistics_json(results_dir::String)
    stats_path = joinpath(results_dir, "statistics.json")

    if !isfile(stats_path)
        println("Warning: statistics.json not found in: $results_dir")
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
    compute_efficiencies(num_events::Int, saved_events::Int, n_events::Int, n_single::Int)

Compute three efficiencies:
1. interaction_efficiency = saved_events / num_events
2. fiducial_efficiency = n_events / saved_events
3. eff_1trk = n_single / n_events

Returns NamedTuple with all efficiencies.
"""
function compute_efficiencies(num_events::Int, saved_events::Int, n_events::Int, n_single::Int)
    interaction_eff = saved_events / num_events
    fiducial_eff = n_events / saved_events
    eff_1trk = n_single / n_events

    return (
        interaction_efficiency = interaction_eff,
        fiducial_efficiency = fiducial_eff,
        eff_1trk = eff_1trk,
        # Also return the raw numbers for summary calculations
        num_events = num_events,
        saved_events = saved_events,
        n_events = n_events,
        n_single = n_single
    )
end

"""
    process_isotope_zone(isotope::String, zone::String, voxel::Float64, maxd::Float64)

Process a single isotope/zone combination and print all three efficiencies.
"""
function process_isotope_zone(isotope::String, zone::String, voxel::Float64, maxd::Float64)
    script_dir = dirname(@__FILE__)

    # Format voxel and maxd (integer if whole number)
    voxel_str = isinteger(voxel) ? string(Int(voxel)) : string(voxel)
    maxd_str = isinteger(maxd) ? string(Int(maxd)) : string(maxd)

    # Build directory name for results
    dir_name = "$(isotope)_$(zone)_filtered_ecut_2400_2500_voxel_$(voxel_str)mm_maxd_$(maxd_str)mm"
    results_dir = joinpath(script_dir, dir_name)

    if !isdir(results_dir)
        println("Warning: Directory not found: $dir_name")
        return nothing
    end

    # Read generation info from CSV
    info = read_info_csv(script_dir, isotope, zone)
    if isnothing(info)
        return nothing
    end

    # Read reconstruction statistics from JSON
    stats = read_statistics_json(results_dir)
    if isnothing(stats)
        return nothing
    end

    # Compute all three efficiencies
    effs = compute_efficiencies(info.num_events, info.saved_events, stats.n_events, stats.n_single)

    # Print results
    println("\n$isotope / $zone:")
    println("  Generation:   num_events=$(info.num_events), saved_events=$(info.saved_events)")
    println("  ROI filtered: n_events=$(stats.n_events), n_single=$(stats.n_single)")
    println("  ---")
    println("  interaction_efficiency = $(info.saved_events)/$(info.num_events) = $(@sprintf("%.2e", effs.interaction_efficiency))")
    println("  fiducial_efficiency    = $(stats.n_events)/$(info.saved_events) = $(@sprintf("%.2e", effs.fiducial_efficiency))")
    println("  eff_1trk               = $(stats.n_single)/$(stats.n_events) = $(@sprintf("%.2e", effs.eff_1trk))")

    return (
        isotope=isotope,
        zone=zone,
        num_events=info.num_events,
        saved_events=info.saved_events,
        n_events=stats.n_events,
        n_single=stats.n_single,
        interaction_efficiency=effs.interaction_efficiency,
        fiducial_efficiency=effs.fiducial_efficiency,
        eff_1trk=effs.eff_1trk
    )
end

#=============================================================================
# Table Generation Functions
=============================================================================#

"""
    get_table_basename(isotope::String, zone::String, voxel::Float64, maxd::Float64)

Get the base name for table files (without extension).
"""
function get_table_basename(isotope::String, zone::String, voxel::Float64, maxd::Float64)
    voxel_str = isinteger(voxel) ? string(Int(voxel)) : string(voxel)
    maxd_str = isinteger(maxd) ? string(Int(maxd)) : string(maxd)
    return "efficiency_table_$(isotope)_$(zone)_voxel_$(voxel_str)mm_maxd_$(maxd_str)mm"
end

"""
    get_table_caption(isotope::String, zone::String, voxel::Float64, maxd::Float64)

Generate a caption describing the table parameters.
"""
function get_table_caption(isotope::String, zone::String, voxel::Float64, maxd::Float64)
    voxel_str = isinteger(voxel) ? string(Int(voxel)) : string(voxel)
    maxd_str = isinteger(maxd) ? string(Int(maxd)) : string(maxd)

    iso_desc = isotope == "all" ? "all isotopes (Bi-214, Tl-208)" : isotope
    zone_desc = zone == "all" ? "all regions (copper shell, copper endcaps)" : replace(zone, "_" => " ")

    return "ITACA efficiency table for $iso_desc in $zone_desc. Voxel size: $(voxel_str) mm, max distance: $(maxd_str) mm."
end

"""
    create_table_dir(isotope::String, zone::String, voxel::Float64, maxd::Float64)

Create directory for efficiency tables.
Returns the path to the created directory.
"""
function create_table_dir(isotope::String, zone::String, voxel::Float64, maxd::Float64)
    script_dir = dirname(@__FILE__)

    # Format voxel and maxd (integer if whole number)
    voxel_str = isinteger(voxel) ? string(Int(voxel)) : string(voxel)
    maxd_str = isinteger(maxd) ? string(Int(maxd)) : string(maxd)

    dir_name = "tables_$(isotope)_$(zone)_voxel_$(voxel_str)mm_maxd_$(maxd_str)mm"
    table_dir = joinpath(script_dir, dir_name)

    if !isdir(table_dir)
        mkpath(table_dir)
        println("Created table directory: $table_dir")
    end

    return table_dir
end

"""
    write_markdown_table(results::Vector, totals::NamedTuple, output_dir::String,
                         isotope::String, zone::String, voxel::Float64, maxd::Float64)

Write efficiency table in Markdown format with caption.
"""
function write_markdown_table(results::Vector, totals::Union{NamedTuple, Nothing}, output_dir::String,
                               isotope::String, zone::String, voxel::Float64, maxd::Float64)
    basename = get_table_basename(isotope, zone, voxel, maxd)
    caption = get_table_caption(isotope, zone, voxel, maxd)
    md_path = joinpath(output_dir, "$(basename).md")

    open(md_path, "w") do f
        # Header with caption
        println(f, "# ITACA Efficiency Table")
        println(f)
        println(f, "**Caption:** $caption")
        println(f)
        println(f, "| Isotope | Zone | ε_int | ε_fid | ε_1trk | ε_total |")
        println(f, "|---------|------|-------|-------|--------|---------|")

        # Data rows
        for r in results
            total_eff = r.interaction_efficiency * r.fiducial_efficiency * r.eff_1trk
            println(f, "| $(r.isotope) | $(r.zone) | $(@sprintf("%.2e", r.interaction_efficiency)) | $(@sprintf("%.2e", r.fiducial_efficiency)) | $(@sprintf("%.2e", r.eff_1trk)) | $(@sprintf("%.2e", total_eff)) |")
        end

        # Total row if multiple results
        if !isnothing(totals)
            total_eff = totals.int_eff * totals.fid_eff * totals.trk_eff
            println(f, "| **TOTAL** | - | $(@sprintf("%.2e", totals.int_eff)) | $(@sprintf("%.2e", totals.fid_eff)) | $(@sprintf("%.2e", totals.trk_eff)) | $(@sprintf("%.2e", total_eff)) |")
        end
    end

    println("Saved Markdown table to: $md_path")
    return md_path
end

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
    write_latex_table(results::Vector, totals::NamedTuple, output_dir::String,
                      isotope::String, zone::String, voxel::Float64, maxd::Float64)

Write efficiency table in LaTeX format with standalone driver and caption.
"""
function write_latex_table(results::Vector, totals::Union{NamedTuple, Nothing}, output_dir::String,
                            isotope::String, zone::String, voxel::Float64, maxd::Float64)
    basename = get_table_basename(isotope, zone, voxel, maxd)
    caption = get_table_caption(isotope, zone, voxel, maxd)
    # Escape underscores in caption for LaTeX
    caption_tex = replace(caption, "_" => "\\_")
    tex_path = joinpath(output_dir, "$(basename).tex")

    open(tex_path, "w") do f
        # LaTeX preamble
        println(f, "\\documentclass{standalone}")
        println(f, "\\usepackage{booktabs}")
        println(f, "\\usepackage{amsmath}")
        println(f, "\\usepackage{caption}")
        println(f, "\\begin{document}")
        println(f, "\\footnotesize")  # Smaller font size
        println(f, "\\begin{minipage}{\\textwidth}")
        println(f, "\\captionof{table}{$caption_tex}")
        println(f, "\\centering")
        println(f, "\\begin{tabular}{llcccc}")
        println(f, "\\toprule")
        println(f, "Isotope & Zone & \$\\epsilon_{\\text{int}}\$ & \$\\epsilon_{\\text{fid}}\$ & \$\\epsilon_{1\\text{trk}}\$ & \$\\epsilon_{\\text{tot}}\$ \\\\")
        println(f, "\\midrule")

        # Data rows
        for r in results
            total_eff = r.interaction_efficiency * r.fiducial_efficiency * r.eff_1trk
            # Escape underscores in zone names for LaTeX
            zone_tex = replace(r.zone, "_" => "\\_")
            println(f, "$(r.isotope) & $zone_tex & $(format_latex_scientific(r.interaction_efficiency)) & $(format_latex_scientific(r.fiducial_efficiency)) & $(format_latex_scientific(r.eff_1trk)) & $(format_latex_scientific(total_eff)) \\\\")
        end

        # Total row if multiple results
        if !isnothing(totals)
            total_eff = totals.int_eff * totals.fid_eff * totals.trk_eff
            println(f, "\\midrule")
            println(f, "\\textbf{TOTAL} & -- & $(format_latex_scientific(totals.int_eff)) & $(format_latex_scientific(totals.fid_eff)) & $(format_latex_scientific(totals.trk_eff)) & $(format_latex_scientific(total_eff)) \\\\")
        end

        # LaTeX closing
        println(f, "\\bottomrule")
        println(f, "\\end{tabular}")
        println(f, "\\end{minipage}")
        println(f, "\\end{document}")
    end

    println("Saved LaTeX table to: $tex_path")
    return tex_path
end

"""
    write_efficiency_tables(results::Vector, totals::Union{NamedTuple, Nothing},
                            isotope::String, zone::String, voxel::Float64, maxd::Float64)

Create directory and write both Markdown and LaTeX efficiency tables.
"""
function write_efficiency_tables(results::Vector, totals::Union{NamedTuple, Nothing},
                                  isotope::String, zone::String, voxel::Float64, maxd::Float64)
    if isempty(results)
        println("Warning: No results to write to table")
        return nothing
    end

    # Create output directory
    output_dir = create_table_dir(isotope, zone, voxel, maxd)

    # Write both table formats
    write_markdown_table(results, totals, output_dir, isotope, zone, voxel, maxd)
    write_latex_table(results, totals, output_dir, isotope, zone, voxel, maxd)

    return output_dir
end

#=============================================================================
# Command Line Interface
=============================================================================#

function main()
    # Parse arguments
    isotope = nothing
    zone = nothing
    voxel = 1.0
    maxd = 4.0
    writetable = false

    for arg in ARGS
        if startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) == 2
                key, value = parts
                if key == "isotope"
                    isotope = String(value)
                elseif key == "zone"
                    zone = String(value)
                elseif key == "voxel"
                    voxel = parse(Float64, value)
                elseif key == "maxd"
                    maxd = parse(Float64, value)
                elseif key == "writetable"
                    writetable = lowercase(value) in ("true", "1", "yes")
                else
                    println("Warning: Unknown argument: --$key")
                end
            end
        end
    end

    # Validate required arguments
    if isnothing(isotope) || isnothing(zone)
        println("Usage: julia itaca_effi.jl --isotope=<isotope> --zone=<zone> [options]")
        println()
        println("Required:")
        println("  --isotope=X     Isotope: bi214, tl208, or all")
        println("  --zone=X        Zone: copper_shell, copper_endcaps, or all")
        println()
        println("Optional:")
        println("  --voxel=X       Voxel size in mm (default: 1)")
        println("  --maxd=X        Max distance in mm (default: 4)")
        println("  --writetable=X  Write efficiency table to files (default: false)")
        println()
        println("Examples:")
        println("  julia itaca_effi.jl --isotope=bi214 --zone=copper_shell")
        println("  julia itaca_effi.jl --isotope=all --zone=all --maxd=3")
        println("  julia itaca_effi.jl --isotope=all --zone=all --maxd=3 --writetable=true")
        exit(1)
    end

    println("=" ^ 70)
    println("ITACA EFFICIENCY ANALYSIS")
    println("=" ^ 70)
    println("Parameters: voxel=$(voxel)mm, maxd=$(maxd)mm")

    # Determine isotopes and zones to process
    isotopes = isotope == "all" ? ["bi214", "tl208"] : [isotope]
    zones = zone == "all" ? ["copper_shell", "copper_endcaps"] : [zone]

    results = []

    for iso in isotopes
        for z in zones
            result = process_isotope_zone(iso, z, voxel, maxd)
            if !isnothing(result)
                push!(results, result)
            end
        end
    end

    println("\n" * "-" ^ 70)

    # Compute totals if multiple results
    totals = nothing
    if length(results) > 1
        total_num_events = sum(r.num_events for r in results)
        total_saved = sum(r.saved_events for r in results)
        total_n_events = sum(r.n_events for r in results)
        total_n_single = sum(r.n_single for r in results)

        total_int_eff = total_saved / total_num_events
        total_fid_eff = total_n_events / total_saved
        total_1trk_eff = total_n_single / total_n_events

        totals = (int_eff=total_int_eff, fid_eff=total_fid_eff, trk_eff=total_1trk_eff)

        println("TOTAL (all isotopes/zones combined):")
        println("  interaction_efficiency = $total_saved/$total_num_events = $(@sprintf("%.2e", total_int_eff))")
        println("  fiducial_efficiency    = $total_n_events/$total_saved = $(@sprintf("%.2e", total_fid_eff))")
        println("  eff_1trk               = $total_n_single/$total_n_events = $(@sprintf("%.2e", total_1trk_eff))")
    end

    println("=" ^ 70)

    # Write tables if requested
    if writetable && !isempty(results)
        println("\nWriting efficiency tables...")
        write_efficiency_tables(results, totals, isotope, zone, voxel, maxd)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
