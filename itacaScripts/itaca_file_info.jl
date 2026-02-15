#!/usr/bin/env julia

"""
Display event generation information from ITACA HDF5 files.

Returns the number of generated events, saved events, and energy cut parameters
stored in the MC/configuration dataset.

Usage:
    julia itaca_file_info.jl <cmdir> <input_file> [options]

Options:
    --writeinfo=X      Write info to CSV file (default: true)

Output:
    num_events:        Number of events generated
    saved_events:      Number of events saved (after cuts)
    interacting_events: Number of events with interactions
    min_energy:        Minimum energy cut (keV)
    max_energy:        Maximum energy cut (keV)

    If --writeinfo=true, creates <basename>_info.csv in itacaScripts/
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using HDF5
using CSV
using DataFrames

"""
    get_file_info(filepath::String)

Extract generation info from HDF5 file configuration.
Returns a Dict with num_events, saved_events, interacting_events, min_energy, max_energy.
"""
function get_file_info(filepath::String)
    info = Dict{String, Any}()

    h5open(filepath, "r") do f
        config = read(f["MC/configuration"])

        for row in config
            key = strip(String(row.param_key))
            val = strip(String(row.param_value))

            if key == "num_events"
                info["num_events"] = parse(Int, val)
            elseif key == "saved_events"
                info["saved_events"] = parse(Int, val)
            elseif key == "interacting_events"
                info["interacting_events"] = parse(Int, val)
            elseif key == "/Actions/DefaultEventAction/min_energy"
                # Parse "2300 keV" -> 2300.0
                info["min_energy_keV"] = parse(Float64, split(val)[1])
            elseif key == "/Actions/DefaultEventAction/max_energy"
                info["max_energy_keV"] = parse(Float64, split(val)[1])
            end
        end
    end

    return info
end

"""
    write_info_csv(info::Dict, input_file::String)

Write info to CSV file in itacaScripts directory.
Returns the path to the created file.
"""
function write_info_csv(info::Dict, input_file::String)
    script_dir = dirname(@__FILE__)
    base_name = splitext(input_file)[1]
    csv_path = joinpath(script_dir, "$(base_name)_info.csv")

    # Create DataFrame with one row
    df = DataFrame(
        file = input_file,
        num_events = get(info, "num_events", missing),
        saved_events = get(info, "saved_events", missing),
        interacting_events = get(info, "interacting_events", missing),
        min_energy_keV = get(info, "min_energy_keV", missing),
        max_energy_keV = get(info, "max_energy_keV", missing)
    )

    CSV.write(csv_path, df)
    println("Saved info to: $csv_path")

    return csv_path
end

"""
    print_file_info(cmdir::String, input_file::String; writeinfo::Bool=true)

Print generation info for an HDF5 file.
If writeinfo=true, also saves to CSV in itacaScripts directory.
"""
function print_file_info(cmdir::String, input_file::String; writeinfo::Bool=true)
    filepath = joinpath(cmdir, input_file)

    if !isfile(filepath)
        error("File not found: $filepath")
    end

    info = get_file_info(filepath)

    println("=" ^ 50)
    println("File: $input_file")
    println("=" ^ 50)
    println("  num_events:         $(get(info, "num_events", "N/A"))")
    println("  saved_events:       $(get(info, "saved_events", "N/A"))")
    println("  interacting_events: $(get(info, "interacting_events", "N/A"))")

    if haskey(info, "min_energy_keV") && haskey(info, "max_energy_keV")
        println("  energy_cut:         [$(info["min_energy_keV"]), $(info["max_energy_keV"])] keV")
    else
        println("  energy_cut:         N/A")
    end
    println("=" ^ 50)

    if writeinfo
        write_info_csv(info, input_file)
    end

    return info
end

#=============================================================================
# Command Line Interface
=============================================================================#

function main()
    if length(ARGS) < 2
        println("Usage: julia itaca_file_info.jl <cmdir> <input_file> [options]")
        println()
        println("Options:")
        println("  --writeinfo=X    Write info to CSV file (default: true)")
        println()
        println("Example:")
        println("  julia itaca_file_info.jl /Users/jjgomezcadenas/Data/HD5t/precdr/bb0nu/ 0nubb.h5")
        println("  julia itaca_file_info.jl /Users/jjgomezcadenas/Data/HD5t/precdr/bb0nu/ 0nubb.h5 --writeinfo=false")
        exit(1)
    end

    cmdir = ARGS[1]
    input_file = ARGS[2]

    # Default values
    writeinfo = true

    # Parse optional arguments
    for arg in ARGS[3:end]
        if startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) == 2
                key, value = parts
                if key == "writeinfo"
                    writeinfo = lowercase(value) in ("true", "1", "yes")
                else
                    println("Warning: Unknown argument: --$key")
                end
            end
        end
    end

    print_file_info(cmdir, input_file; writeinfo=writeinfo)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
