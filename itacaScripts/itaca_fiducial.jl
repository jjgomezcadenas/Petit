#!/usr/bin/env julia

"""
Fiducial efficiency calculator.

Reads CSV info files (geometry cut) and TXT info files (energy cut) from the
datainfo directory and computes the combined fiducial efficiency.

Usage:
    julia itaca_fiducial.jl [options]

Options:
    --datadir=X     Directory containing info files (default: datainfo)
    --output=X      Output CSV file name (default: fiducial_efficiency.csv)
    --help, -h      Show this help message

Output:
    CSV file with columns:
    - source: source name (e.g., bi214_copper_shell)
    - num_events_tot: total generated events
    - events_fid_geom: events after geometry cut
    - fiducial_geom_efficiency: geometry cut efficiency
    - num_events_cut: events after energy cut
    - fiducial_cut: energy cut efficiency
    - fiducial_tot: total fiducial efficiency
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
# File Parsing Functions
# =============================================================================

"""
    parse_energy_cut_txt(filepath::String)

Parse an energy cut info TXT file and extract event statistics.
Returns (total_events_read, events_passing_cut) or nothing if parsing fails.
"""
function parse_energy_cut_txt(filepath::String)
    if !isfile(filepath)
        return nothing
    end

    total_read = nothing
    passing_cut = nothing

    try
        for line in eachline(filepath)
            # Look for "Total events read:" line
            m = match(r"Total events read:\s+(\d+)", line)
            if m !== nothing
                total_read = parse(Int, m.captures[1])
            end

            # Look for "Events passing cut:" line
            m = match(r"Events passing cut:\s+(\d+)", line)
            if m !== nothing
                passing_cut = parse(Int, m.captures[1])
            end
        end

        if total_read !== nothing && passing_cut !== nothing
            return (total_read=total_read, passing_cut=passing_cut)
        else
            return nothing
        end
    catch e
        @warn "Error parsing $filepath: $e"
        return nothing
    end
end

"""
    extract_source_from_csv_filename(filename::String)

Extract source name from CSV info filename.
Examples:
  - "bi214_copper_shell_info.csv" → "bi214_copper_shell"
  - "bi214_ptfe_barrel_1.next_info.csv" → "bi214_ptfe_barrel_1.next"
"""
function extract_source_from_csv_filename(filename::String)
    source = replace(filename, r"_info\.csv$" => "")
    source = replace(source, r"\.next_info\.csv$" => "")
    return source
end

"""
    extract_source_from_txt_filename(filename::String)

Extract source name from TXT info filename.
Examples:
  - "bi214_copper_shell_filtered_ecut_2400_2500_info.txt" → "bi214_copper_shell"
  - "bi214_ptfe_barrel_1.next_filtered_ecut_2400_2500_info.txt" → "bi214_ptfe_barrel_1.next"
"""
function extract_source_from_txt_filename(filename::String)
    # Remove the _filtered_ecut_XXXX_XXXX_info.txt suffix
    source = replace(filename, r"_filtered_ecut_\d+_\d+_info\.txt$" => "")
    return source
end

"""
    normalize_source_name(source::String)

Normalize source name for matching (remove .next suffix variations).
"""
function normalize_source_name(source::String)
    # For matching purposes, treat "bi214_ptfe_barrel_1.next" same as "bi214_ptfe_barrel"
    # But keep the original for output
    normalized = replace(source, r"_\d+\.next$" => "")
    normalized = replace(normalized, r"\.next$" => "")
    return normalized
end

"""
    read_geometry_cut_files(datadir::String)

Read all *_info.csv files from the specified directory.
Returns a Dict mapping source name to (num_events, saved_events).
"""
function read_geometry_cut_files(datadir::String)
    all_files = readdir(datadir)
    csv_files = filter(f -> endswith(f, "_info.csv") || endswith(f, ".next_info.csv"), all_files)

    results = Dict{String, NamedTuple{(:num_events, :saved_events), Tuple{Int, Int}}}()

    for f in csv_files
        filepath = joinpath(datadir, f)
        df = CSV.read(filepath, DataFrame)

        if nrow(df) > 0
            source = extract_source_from_csv_filename(f)
            num_events = df.num_events[1]
            saved_events = df.saved_events[1]
            results[source] = (num_events=num_events, saved_events=saved_events)
        end
    end

    return results
end

"""
    read_energy_cut_files(datadir::String)

Read all *_filtered_ecut_*_info.txt files from the specified directory.
Returns a Dict mapping source name to (total_read, passing_cut).
"""
function read_energy_cut_files(datadir::String)
    all_files = readdir(datadir)
    txt_files = filter(f -> occursin(r"_filtered_ecut_\d+_\d+_info\.txt$", f), all_files)

    results = Dict{String, NamedTuple{(:total_read, :passing_cut), Tuple{Int, Int}}}()

    for f in txt_files
        filepath = joinpath(datadir, f)
        parsed = parse_energy_cut_txt(filepath)

        if parsed !== nothing
            source = extract_source_from_txt_filename(f)
            results[source] = (total_read=parsed.total_read, passing_cut=parsed.passing_cut)
        end
    end

    return results
end

"""
    combine_efficiencies(geom_data::Dict, ecut_data::Dict)

Combine geometry and energy cut data into a results DataFrame.
"""
function combine_efficiencies(geom_data::Dict, ecut_data::Dict)
    results = DataFrame(
        source = String[],
        num_events_tot = Int[],
        events_fid_geom = Int[],
        fiducial_geom_efficiency = Float64[],
        num_events_cut = Int[],
        fiducial_cut = Float64[],
        fiducial_tot = Float64[]
    )

    # Process all sources that have geometry data
    for (source, geom) in sort(collect(geom_data), by=x->x[1])
        num_events_tot = geom.num_events
        events_fid_geom = geom.saved_events
        fiducial_geom_eff = events_fid_geom / num_events_tot

        # Try to find matching energy cut data
        # First try exact match, then try normalized match
        ecut = get(ecut_data, source, nothing)
        if ecut === nothing
            # Try normalized matching
            norm_source = normalize_source_name(source)
            for (ecut_source, ecut_val) in ecut_data
                if normalize_source_name(ecut_source) == norm_source
                    ecut = ecut_val
                    break
                end
            end
        end

        if ecut !== nothing
            # Validate: total_read should match saved_events (geometry cut output)
            if ecut.total_read != events_fid_geom
                @warn "Mismatch for $source: geometry saved_events=$(events_fid_geom), energy total_read=$(ecut.total_read)"
            end

            num_events_cut = ecut.passing_cut
            fiducial_cut = num_events_cut / events_fid_geom
            fiducial_tot = num_events_cut / num_events_tot
        else
            # No energy cut data available
            num_events_cut = 0
            fiducial_cut = 0.0
            fiducial_tot = 0.0
            @warn "No energy cut data found for: $source"
        end

        push!(results, (
            source,
            num_events_tot,
            events_fid_geom,
            fiducial_geom_eff,
            num_events_cut,
            fiducial_cut,
            fiducial_tot
        ))
    end

    return results
end

# =============================================================================
# Main
# =============================================================================

function main()
    # Default values
    datadir = joinpath(dirname(@__FILE__), "datainfo")
    output_file = "fiducial_efficiency.csv"

    # Parse arguments
    for arg in ARGS
        if startswith(arg, "--datadir=")
            datadir = String(split(arg, '=')[2])
        elseif startswith(arg, "--output=")
            output_file = String(split(arg, '=')[2])
        elseif arg == "--help" || arg == "-h"
            println("""
Usage: julia itaca_fiducial.jl [options]

Options:
    --datadir=X     Directory containing info files (default: datainfo)
    --output=X      Output CSV file name (default: fiducial_efficiency.csv)
    --help, -h      Show this help message

Output:
    CSV file with columns:
    - source: source name
    - num_events_tot: total generated events
    - events_fid_geom: events after geometry cut
    - fiducial_geom_efficiency: geometry cut efficiency
    - num_events_cut: events after energy cut (2400-2500 keV)
    - fiducial_cut: energy cut efficiency
    - fiducial_tot: total fiducial efficiency (geometry × energy)
""")
            return
        else
            println("Warning: Unknown argument: $arg")
        end
    end

    if !isdir(datadir)
        error("Directory not found: $datadir")
    end

    println("=" ^ 80)
    println("FIDUCIAL EFFICIENCY CALCULATOR (Geometry + Energy Cut)")
    println("=" ^ 80)
    println("Data directory: $datadir")
    println("Output file:    $output_file")
    println("=" ^ 80)
    println()

    # Read geometry cut data (CSV files)
    println("Reading geometry cut data (CSV files)...")
    geom_data = read_geometry_cut_files(datadir)
    println("  Found $(length(geom_data)) sources with geometry cut data")
    println()

    # Read energy cut data (TXT files)
    println("Reading energy cut data (TXT files)...")
    ecut_data = read_energy_cut_files(datadir)
    println("  Found $(length(ecut_data)) sources with energy cut data")
    println()

    # Combine efficiencies
    println("Combining efficiencies...")
    results = combine_efficiencies(geom_data, ecut_data)

    if nrow(results) == 0
        println("No data found.")
        return
    end

    # Save results
    output_path = joinpath(dirname(@__FILE__), output_file)
    CSV.write(output_path, results)
    println()
    println("Saved results to: $output_path")

    # Print summary table
    println()
    println("=" ^ 120)
    println("SUMMARY")
    println("=" ^ 120)
    println()
    @printf("%-30s  %12s  %12s  %12s  %12s  %12s  %12s\n",
            "Source", "N_tot", "N_geom", "ε_geom", "N_cut", "ε_cut", "ε_tot")
    println("-" ^ 120)

    for row in eachrow(results)
        @printf("%-30s  %12d  %12d  %12.6f  %12d  %12.6f  %12.6f\n",
                row.source,
                row.num_events_tot,
                row.events_fid_geom,
                row.fiducial_geom_efficiency,
                row.num_events_cut,
                row.fiducial_cut,
                row.fiducial_tot)
    end

    println("=" ^ 120)
    println()
    println("Legend:")
    println("  N_tot:   Total generated events")
    println("  N_geom:  Events after geometry (fiducial) cut")
    println("  ε_geom:  Geometry cut efficiency (N_geom / N_tot)")
    println("  N_cut:   Events after energy cut (2400-2500 keV)")
    println("  ε_cut:   Energy cut efficiency (N_cut / N_geom)")
    println("  ε_tot:   Total fiducial efficiency (N_cut / N_tot)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
