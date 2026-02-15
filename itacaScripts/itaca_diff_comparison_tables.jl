#!/usr/bin/env julia

"""
Generate diffusion comparison efficiency tables.

Reads statistics.json from reconstruction results and info CSV files from data_info/
to compute efficiencies for different diffusion configurations.
Outputs a single LaTeX document with 8 tables (one per isotope/zone/diffusion combo).

Hardwired configurations:
- Isotopes: tl208, bi214
- Zones: copper_shell, copper_endcaps
- Diffusion settings: (σt=1mm, σl=1mm) and (σt=10mm, σl=3mm)

Usage:
    julia itaca_diff_comparison_tables.jl

Output: tables_summary_diff_comparison/diffusion_comparison_tables.tex
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using JSON
using CSV
using DataFrames
using Printf

#=============================================================================
# Configuration
=============================================================================#

# Hardwired configurations
const ISOTOPES = ["tl208", "bi214"]
const ZONES = ["copper_shell", "copper_endcaps"]

# Additional configurations with custom directory naming
# Format: (isotope, zone, info_file_zone, config) -> custom_dir_name
const ADDITIONAL_CONFIGS = [
    (isotope="bi214", zone="ptfe_endcap_1.next", info_zone="ptfe_endcap_1.next",
     config=(sigmat=1.0, sigmal=1.0, voxel=2.0, maxd=4.0),
     dir_name="bi214_ptfe_endcap_1.next_sigmat_1mm_sigmal_1mm_voxel_2mm_maxd_4mm"),
    (isotope="bi214", zone="ptfe_endcap_1.next", info_zone="ptfe_endcap_1.next",
     config=(sigmat=10.0, sigmal=3.0, voxel=20.0, maxd=21.0),
     dir_name="bi214_ptfe_endcap_1.next_sigmat_10mm_sigmal_3mm_voxel_20mm_maxd_21mm")
]

# Diffusion settings: (sigmat, sigmal, voxel, maxd)
const DIFFUSION_CONFIGS = [
    (sigmat=1.0,  sigmal=1.0, voxel=2.0,  maxd=4.0),
    (sigmat=10.0, sigmal=3.0, voxel=20.0, maxd=21.0)
]

#=============================================================================
# Helper Functions
=============================================================================#

"""
    get_results_dir_name(isotope, zone, config)

Build the results directory name for a given configuration.
Format: <isotope>_<zone>_filtered_ecut_2400_2500_sigmat_Xmm_sigmal_Ymm_voxel_Zmm_maxd_Wmm
"""
function get_results_dir_name(isotope::String, zone::String, config)
    sigmat_str = isinteger(config.sigmat) ? string(Int(config.sigmat)) : string(config.sigmat)
    sigmal_str = isinteger(config.sigmal) ? string(Int(config.sigmal)) : string(config.sigmal)
    voxel_str = isinteger(config.voxel) ? string(Int(config.voxel)) : string(config.voxel)
    maxd_str = isinteger(config.maxd) ? string(Int(config.maxd)) : string(config.maxd)

    return "$(isotope)_$(zone)_filtered_ecut_2400_2500_sigmat_$(sigmat_str)mm_sigmal_$(sigmal_str)mm_voxel_$(voxel_str)mm_maxd_$(maxd_str)mm"
end

"""
    read_info_csv(script_dir, isotope, zone; info_zone=nothing)

Read generation info from CSV file in the data_info directory.
Returns num_events, interacting_events, saved_events.
If info_zone is provided, use that instead of zone for the filename.
"""
function read_info_csv(script_dir::String, isotope::String, zone::String; info_zone::Union{String,Nothing}=nothing)
    zone_for_file = info_zone !== nothing ? info_zone : zone
    csv_path = joinpath(script_dir, "data_info", "$(isotope)_$(zone_for_file)_info.csv")

    if !isfile(csv_path)
        error("Info CSV not found: $csv_path")
    end

    df = CSV.read(csv_path, DataFrame)

    if nrow(df) == 0
        error("Empty CSV file: $csv_path")
    end

    # Use first row (should only be one row per file)
    row = df[1, :]

    return (
        num_events = row.num_events,
        interacting_events = row.interacting_events,
        saved_events = row.saved_events
    )
end

"""
    read_statistics_json(script_dir, isotope, zone, config; custom_dir=nothing)

Read statistics.json from the results directory.
Returns track multiplicity data.
If custom_dir is provided, use that instead of constructing from isotope/zone/config.
"""
function read_statistics_json(script_dir::String, isotope::String, zone::String, config; custom_dir::Union{String,Nothing}=nothing)
    dir_name = custom_dir !== nothing ? custom_dir : get_results_dir_name(isotope, zone, config)
    json_path = joinpath(script_dir, dir_name, "statistics.json")

    if !isfile(json_path)
        error("Statistics JSON not found: $json_path")
    end

    data = JSON.parsefile(json_path)

    # Extract track multiplicity
    n_events = data["metadata"]["n_events"]
    track_mult = data["track_multiplicity"]["counts"]

    # Count 1-track events
    n_1track = get(track_mult, "1", 0)

    return (
        n_events = n_events,
        n_1track = n_1track,
        eff_1trk = n_1track / n_events
    )
end

"""
    compute_efficiencies(script_dir, isotope, zone, config; info_zone=nothing, custom_dir=nothing)

Compute all efficiencies for a given configuration.
Optional info_zone and custom_dir for non-standard directory naming.
"""
function compute_efficiencies(script_dir::String, isotope::String, zone::String, config;
                              info_zone::Union{String,Nothing}=nothing,
                              custom_dir::Union{String,Nothing}=nothing)
    # Read base info from CSV (independent of diffusion)
    info = read_info_csv(script_dir, isotope, zone; info_zone=info_zone)

    # Read track statistics (depends on diffusion)
    stats = read_statistics_json(script_dir, isotope, zone, config; custom_dir=custom_dir)

    # Compute efficiencies
    # eff_fid = fraction of generated events passing energy cut (2400-2500 keV)
    eff_fid = info.saved_events / info.num_events
    eff_1trk = stats.eff_1trk
    eff_total = eff_fid * eff_1trk

    return (
        isotope = isotope,
        zone = zone,
        sigmat = config.sigmat,
        sigmal = config.sigmal,
        eff_fid = eff_fid,
        eff_1trk = eff_1trk,
        eff_total = eff_total
    )
end

"""
    compute_total_efficiencies(script_dir, isotope, config)

Compute total efficiencies combining shell + endcaps for a given isotope and diffusion config.
Uses proper weighted averages based on raw counts.
"""
function compute_total_efficiencies(script_dir::String, isotope::String, config)
    # Read info for both zones
    info_shell = read_info_csv(script_dir, isotope, "copper_shell")
    info_endcaps = read_info_csv(script_dir, isotope, "copper_endcaps")

    # Read statistics for both zones
    stats_shell = read_statistics_json(script_dir, isotope, "copper_shell", config)
    stats_endcaps = read_statistics_json(script_dir, isotope, "copper_endcaps", config)

    # Combine raw counts
    total_num_events = info_shell.num_events + info_endcaps.num_events
    total_saved = info_shell.saved_events + info_endcaps.saved_events
    total_n_events_reco = stats_shell.n_events + stats_endcaps.n_events
    total_n_1track = stats_shell.n_1track + stats_endcaps.n_1track

    # Compute combined efficiencies
    # eff_fid = fraction of generated events passing energy cut (2400-2500 keV)
    eff_fid = total_saved / total_num_events
    eff_1trk = total_n_1track / total_n_events_reco
    eff_total = eff_fid * eff_1trk

    return (
        isotope = isotope,
        zone = "total",
        sigmat = config.sigmat,
        sigmal = config.sigmal,
        eff_fid = eff_fid,
        eff_1trk = eff_1trk,
        eff_total = eff_total,
        # Store raw counts for reference
        num_events = total_num_events,
        saved_events = total_saved,
        n_events_reco = total_n_events_reco,
        n_1track = total_n_1track
    )
end

"""
    format_latex_scientific(x::Float64)

Format a number in LaTeX scientific notation.
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
    format_isotope_name(isotope::String)

Format isotope name for display (e.g., "tl208" -> "Tl-208").
"""
function format_isotope_name(isotope::String)
    if isotope == "tl208"
        return "Tl-208"
    elseif isotope == "bi214"
        return "Bi-214"
    else
        return isotope
    end
end

"""
    format_zone_name(zone::String)

Format zone name for display (e.g., "copper_shell" -> "copper shell").
"""
function format_zone_name(zone::String)
    return replace(zone, "_" => " ")
end

#=============================================================================
# LaTeX Generation
=============================================================================#

"""
    write_latex_document(results, totals, output_path)

Write a complete LaTeX document with all efficiency tables.
Includes individual zone tables, total tables, and ratio comparison.
"""
function write_latex_document(results::Vector, totals::Vector, output_path::String)
    open(output_path, "w") do f
        # LaTeX preamble
        println(f, "\\documentclass{article}")
        println(f, "\\usepackage{booktabs}")
        println(f, "\\usepackage{amsmath}")
        println(f, "\\usepackage{caption}")
        println(f, "\\usepackage[margin=1in]{geometry}")
        println(f, "")
        println(f, "\\begin{document}")
        println(f, "")
        println(f, "\\section{Efficiency Calculations}")
        println(f, "")

        for isotope in ISOTOPES
            for zone in ZONES
                for config in DIFFUSION_CONFIGS
                    # Find matching result
                    r = nothing
                    for res in results
                        if res.isotope == isotope && res.zone == zone &&
                           res.sigmat == config.sigmat && res.sigmal == config.sigmal
                            r = res
                            break
                        end
                    end

                    if r === nothing
                        println("Warning: No result found for $isotope / $zone / σt=$(config.sigmat)")
                        continue
                    end

                    # Generate table
                    write_single_table(f, r)
                    println(f, "")
                end
            end
        end

        # Tables for additional configurations
        for add_cfg in ADDITIONAL_CONFIGS
            # Find matching result
            r = nothing
            for res in results
                if res.isotope == add_cfg.isotope && res.zone == add_cfg.zone &&
                   res.sigmat == add_cfg.config.sigmat && res.sigmal == add_cfg.config.sigmal
                    r = res
                    break
                end
            end

            if r === nothing
                println("Warning: No result found for $(add_cfg.isotope) / $(add_cfg.zone)")
                continue
            end

            # Generate table
            write_single_table(f, r)
            println(f, "")
        end

        # Total tables (4 tables)

        for isotope in ISOTOPES
            for config in DIFFUSION_CONFIGS
                # Find matching total
                t = nothing
                for tot in totals
                    if tot.isotope == isotope &&
                       tot.sigmat == config.sigmat && tot.sigmal == config.sigmal
                        t = tot
                        break
                    end
                end

                if t === nothing
                    println("Warning: No total found for $isotope / σt=$(config.sigmat)")
                    continue
                end

                write_total_table(f, t)
                println(f, "")
            end
        end

        # Ratio comparison table
        write_ratio_table(f, totals)
        println(f, "")

        # Close document
        println(f, "\\end{document}")
    end

    println("Saved LaTeX document to: $output_path")
end

"""
    write_single_table(f, r)

Write a single efficiency table to the file handle.
"""
function write_single_table(f::IO, r)
    iso_name = format_isotope_name(r.isotope)
    zone_name = format_zone_name(r.zone)
    zone_tex = replace(r.zone, "_" => "\\_")

    sigmat_str = isinteger(r.sigmat) ? string(Int(r.sigmat)) : string(r.sigmat)
    sigmal_str = isinteger(r.sigmal) ? string(Int(r.sigmal)) : string(r.sigmal)

    caption = "Efficiency table for $iso_name in $zone_name. Diffusion: \$\\sigma_t=$(sigmat_str)\$ mm, \$\\sigma_l=$(sigmal_str)\$ mm."

    println(f, "\\begin{table}[h!]")
    println(f, "\\centering")
    println(f, "\\caption{$caption}")
    println(f, "\\begin{tabular}{lccc}")
    println(f, "\\toprule")
    println(f, "Isotope & \$\\epsilon_{\\text{fid}}\$ & \$\\epsilon_{1\\text{trk}}\$ & \$\\epsilon_{\\text{tot}}\$ \\\\")
    println(f, "\\midrule")

    # Data row
    println(f, "$(r.isotope) & $(format_latex_scientific(r.eff_fid)) & $(format_latex_scientific(r.eff_1trk)) & $(format_latex_scientific(r.eff_total)) \\\\")

    println(f, "\\bottomrule")
    println(f, "\\end{tabular}")
    println(f, "\\end{table}")
end

"""
    write_total_table(f, r)

Write a total efficiency table (shell + endcaps combined) to the file handle.
"""
function write_total_table(f::IO, r)
    iso_name = format_isotope_name(r.isotope)

    sigmat_str = isinteger(r.sigmat) ? string(Int(r.sigmat)) : string(r.sigmat)
    sigmal_str = isinteger(r.sigmal) ? string(Int(r.sigmal)) : string(r.sigmal)

    caption = "Total efficiency table for $iso_name (shell + endcaps). Diffusion: \$\\sigma_t=$(sigmat_str)\$ mm, \$\\sigma_l=$(sigmal_str)\$ mm."

    println(f, "\\begin{table}[h!]")
    println(f, "\\centering")
    println(f, "\\caption{$caption}")
    println(f, "\\begin{tabular}{lccc}")
    println(f, "\\toprule")
    println(f, "Isotope & \$\\epsilon_{\\text{fid}}\$ & \$\\epsilon_{1\\text{trk}}\$ & \$\\epsilon_{\\text{tot}}\$ \\\\")
    println(f, "\\midrule")

    # Data row
    println(f, "$(r.isotope) & $(format_latex_scientific(r.eff_fid)) & $(format_latex_scientific(r.eff_1trk)) & $(format_latex_scientific(r.eff_total)) \\\\")

    println(f, "\\bottomrule")
    println(f, "\\end{tabular}")
    println(f, "\\end{table}")
end

"""
    write_ratio_table(f, totals)

Write a ratio comparison table (large diffusion / small diffusion).
"""
function write_ratio_table(f::IO, totals::Vector)
    println(f, "\\begin{table}[h!]")
    println(f, "\\centering")
    println(f, "\\caption{Efficiency ratios: large diffusion (\$\\sigma_t=10\$ mm) / small diffusion (\$\\sigma_t=1\$ mm).}")
    println(f, "\\begin{tabular}{lccc}")
    println(f, "\\toprule")
    println(f, "Isotope & \$\\epsilon_{1\\text{trk}}\$ ratio & \$\\epsilon_{\\text{tot}}\$ ratio \\\\")
    println(f, "\\midrule")

    for isotope in ISOTOPES
        # Find small and large diffusion results for this isotope
        small_diff = nothing
        large_diff = nothing
        for t in totals
            if t.isotope == isotope
                if t.sigmat == 1.0
                    small_diff = t
                elseif t.sigmat == 10.0
                    large_diff = t
                end
            end
        end

        if small_diff !== nothing && large_diff !== nothing
            ratio_1trk = large_diff.eff_1trk / small_diff.eff_1trk
            ratio_total = large_diff.eff_total / small_diff.eff_total
            iso_name = format_isotope_name(isotope)
            println(f, "$iso_name & $(@sprintf("%.2f", ratio_1trk)) & $(@sprintf("%.2f", ratio_total)) \\\\")
        end
    end

    println(f, "\\bottomrule")
    println(f, "\\end{tabular}")
    println(f, "\\end{table}")
end

#=============================================================================
# Main
=============================================================================#

function main()
    script_dir = dirname(@__FILE__)
    output_dir = joinpath(script_dir, "tables_summary_diff_comparison")

    # Create output directory
    if !isdir(output_dir)
        mkpath(output_dir)
        println("Created output directory: $output_dir")
    end

    println("=" ^ 70)
    println("DIFFUSION COMPARISON EFFICIENCY TABLES")
    println("=" ^ 70)
    println()
    println("Script directory: $script_dir")
    println("Info files: $(joinpath(script_dir, "data_info"))")
    println()

    # Collect all results
    results = []

    for isotope in ISOTOPES
        for zone in ZONES
            println("Processing: $isotope / $zone")

            for config in DIFFUSION_CONFIGS
                println("  σt=$(config.sigmat) mm, σl=$(config.sigmal) mm...")

                try
                    r = compute_efficiencies(script_dir, isotope, zone, config)
                    push!(results, r)

                    println("    ε_fid = $(@sprintf("%.2e", r.eff_fid))")
                    println("    ε_1trk = $(@sprintf("%.2e", r.eff_1trk))")
                    println("    ε_total = $(@sprintf("%.2e", r.eff_total))")
                catch e
                    println("    ERROR: $e")
                end
            end
        end
    end

    # Process additional configurations with custom directory names
    println("\nProcessing additional configurations...")
    for add_cfg in ADDITIONAL_CONFIGS
        println("Processing: $(add_cfg.isotope) / $(add_cfg.zone)")
        println("  σt=$(add_cfg.config.sigmat) mm, σl=$(add_cfg.config.sigmal) mm...")
        println("  Directory: $(add_cfg.dir_name)")

        try
            r = compute_efficiencies(script_dir, add_cfg.isotope, add_cfg.zone, add_cfg.config;
                                     info_zone=add_cfg.info_zone, custom_dir=add_cfg.dir_name)
            push!(results, r)

            println("    ε_fid = $(@sprintf("%.2e", r.eff_fid))")
            println("    ε_1trk = $(@sprintf("%.2e", r.eff_1trk))")
            println("    ε_total = $(@sprintf("%.2e", r.eff_total))")
        catch e
            println("    ERROR: $e")
        end
    end

    println()
    println("=" ^ 70)

    if isempty(results)
        println("No results to write!")
        return
    end

    # Compute total efficiencies (shell + endcaps combined)
    println("\nComputing total efficiencies (shell + endcaps)...")
    totals = []

    for isotope in ISOTOPES
        for config in DIFFUSION_CONFIGS
            println("  Total $isotope / σt=$(config.sigmat) mm...")

            try
                t = compute_total_efficiencies(script_dir, isotope, config)
                push!(totals, t)

                println("    ε_fid = $(@sprintf("%.2e", t.eff_fid))")
                println("    ε_1trk = $(@sprintf("%.2e", t.eff_1trk))")
                println("    ε_total = $(@sprintf("%.2e", t.eff_total))")
            catch e
                println("    ERROR: $e")
            end
        end
    end

    # Print ratio comparison (large/small)
    println("\nDiffusion comparison ratios (large/small):")
    for isotope in ISOTOPES
        small_diff = nothing
        large_diff = nothing
        for t in totals
            if t.isotope == isotope
                if t.sigmat == 1.0
                    small_diff = t
                elseif t.sigmat == 10.0
                    large_diff = t
                end
            end
        end

        if small_diff !== nothing && large_diff !== nothing
            ratio_1trk = large_diff.eff_1trk / small_diff.eff_1trk
            ratio_total = large_diff.eff_total / small_diff.eff_total
            println("  $(format_isotope_name(isotope)): ε_1trk ratio = $(@sprintf("%.2f", ratio_1trk)), ε_total ratio = $(@sprintf("%.2f", ratio_total))")
        end
    end

    println()
    println("=" ^ 70)

    # Write LaTeX document
    output_path = joinpath(output_dir, "diffusion_comparison_tables.tex")
    write_latex_document(results, totals, output_path)

    println()
    println("Generated $(length(results)) individual tables + $(length(totals)) total tables + 1 ratio table.")
    println("=" ^ 70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
