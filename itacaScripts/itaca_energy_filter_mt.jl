#!/usr/bin/env julia

"""
Multi-threaded energy filter for ITACA HDF5 files.

Filters events by total energy and writes filtered data to HDF5 files.

Features:
- Filters events where total energy is in [emin, emax] keV
- Multi-threaded writing (one file per thread)
- Preserves HDF5 structure identical to input
- Creates summary file with filter statistics

Usage:
    julia -t <nthreads> itaca_energy_filter_mt.jl <cmdir> <input_file> [options]

Options:
    --emin=X            Minimum energy in keV (default: 2400)
    --emax=X            Maximum energy in keV (default: 2500)
    --nthreads=N        Number of threads (default: 1)

Output:
    Creates directory: <cmdir>/<basename>_filtered_ecut_<emin>_<emax>/
    Contains:
    - <basename>_part_1.h5, <basename>_part_2.h5, ... (filtered data, one per thread)
    - filter_summary.txt (statistics and parameters)
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using HDF5
using DataFrames
using Statistics
using Dates

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

#=============================================================================
# Helper Functions
=============================================================================#

"""
    compute_event_energies(hitsdf::DataFrame)

Compute total energy (in keV) for each event.
Returns a DataFrame with columns: event_id, total_energy_keV
"""
function compute_event_energies(hitsdf::DataFrame)
    event_energies = combine(groupby(hitsdf, :event_id)) do df
        (total_energy_keV = 1e3 * sum(df.energy),)  # Convert MeV to keV
    end
    return event_energies
end

"""
    filter_events_by_energy(event_energies::DataFrame, emin::Float64, emax::Float64)

Filter events by total energy in [emin, emax] keV.
Returns vector of event IDs that pass the cut.
"""
function filter_events_by_energy(event_energies::DataFrame, emin::Float64, emax::Float64)
    mask = (event_energies.total_energy_keV .>= emin) .& (event_energies.total_energy_keV .<= emax)
    return event_energies.event_id[mask]
end

"""
    create_output_dir(cmdir::String, base_name::String, emin::Float64, emax::Float64)

Create output directory for filtered data.
Returns the path to the created directory.
"""
function create_output_dir(cmdir::String, base_name::String, emin::Float64, emax::Float64)
    # Format energy values (integer if whole number)
    emin_str = isinteger(emin) ? string(Int(emin)) : string(emin)
    emax_str = isinteger(emax) ? string(Int(emax)) : string(emax)

    output_dir = joinpath(cmdir, "$(base_name)_filtered_ecut_$(emin_str)_$(emax_str)")

    if !isdir(output_dir)
        mkpath(output_dir)
        println("Created output directory: $output_dir")
    else
        println("Output directory exists: $output_dir")
    end

    return output_dir
end

"""
    write_filtered_hits_hdf5(output_path::String, hits_df::DataFrame, metadata::Dict)

Write filtered hits to HDF5 file, preserving structure similar to input.
"""
function write_filtered_hits_hdf5(output_path::String, hits_df::DataFrame, metadata::Dict)
    h5open(output_path, "w") do fid
        # Write metadata as attributes
        for (key, val) in metadata
            attrs(fid)[key] = string(val)
        end

        # Create hits group and write data
        hits_group = create_group(fid, "MC/hits")

        if nrow(hits_df) > 0
            # Write each column
            for col in names(hits_df)
                data = hits_df[!, col]
                if eltype(data) <: AbstractString
                    hits_group[col] = String.(data)
                else
                    hits_group[col] = Vector(data)
                end
            end
        end

        # Store column names for reference
        attrs(hits_group)["columns"] = String.(names(hits_df))
        attrs(hits_group)["n_hits"] = nrow(hits_df)
        attrs(hits_group)["n_events"] = length(unique(hits_df.event_id))
    end
end

"""
    write_summary_file(output_dir::String, input_file::String, emin::Float64, emax::Float64,
                       n_total_events::Int, n_filtered_events::Int,
                       energy_stats::NamedTuple, thread_results::Vector)

Write summary file with filter statistics.
"""
function write_summary_file(output_dir::String, input_file::String,
                            emin::Float64, emax::Float64,
                            n_total_events::Int, n_filtered_events::Int,
                            energy_stats::NamedTuple, thread_results::Vector)

    summary_path = joinpath(output_dir, "filter_summary.txt")

    open(summary_path, "w") do f
        println(f, "=" ^ 70)
        println(f, "ITACA ENERGY FILTER SUMMARY")
        println(f, "=" ^ 70)
        println(f)
        println(f, "Timestamp: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(f)
        println(f, "Input/Output:")
        println(f, "  Input file:      $input_file")
        println(f, "  Output directory: $output_dir")
        println(f)
        println(f, "Filter Parameters:")
        println(f, "  Energy min:      $emin keV")
        println(f, "  Energy max:      $emax keV")
        println(f)
        println(f, "Event Statistics:")
        println(f, "  Total events read:     $n_total_events")
        println(f, "  Events passing cut:    $n_filtered_events")
        println(f, "  Events rejected:       $(n_total_events - n_filtered_events)")
        println(f, "  Filter efficiency:     $(round(100 * n_filtered_events / n_total_events, digits=2))%")
        println(f)
        println(f, "Energy Distribution (all events):")
        println(f, "  Mean:   $(round(energy_stats.mean, digits=2)) keV")
        println(f, "  Std:    $(round(energy_stats.std, digits=2)) keV")
        println(f, "  Min:    $(round(energy_stats.min, digits=2)) keV")
        println(f, "  Max:    $(round(energy_stats.max, digits=2)) keV")
        println(f)
        println(f, "Output Files:")
        for r in thread_results
            println(f, "  $(r["output_file"]): $(r["n_events"]) events, $(r["n_hits"]) hits")
        end
        println(f)
        println(f, "=" ^ 70)
    end

    println("Summary saved to: $summary_path")
    return summary_path
end

#=============================================================================
# Main Processing Functions
=============================================================================#

"""
    filter_and_write_thread(hitsdf::DataFrame, event_ids::Vector{Int},
                            thread_id::Int, output_dir::String, base_name::String, metadata::Dict)

Process events for a single thread: filter hits and write to HDF5.
"""
function filter_and_write_thread(hitsdf::DataFrame, event_ids::Vector{Int},
                                  thread_id::Int, output_dir::String, base_name::String, metadata::Dict)
    # Filter hits for assigned events
    hits_filtered = filter(row -> row.event_id in event_ids, hitsdf)

    # Write to HDF5
    output_file = "$(base_name)_part_$thread_id.h5"
    output_path = joinpath(output_dir, output_file)

    thread_metadata = copy(metadata)
    thread_metadata["thread_id"] = thread_id
    thread_metadata["n_events_in_part"] = length(event_ids)
    thread_metadata["n_hits_in_part"] = nrow(hits_filtered)

    write_filtered_hits_hdf5(output_path, hits_filtered, thread_metadata)

    return Dict(
        "thread_id" => thread_id,
        "output_file" => output_file,
        "n_events" => length(event_ids),
        "n_hits" => nrow(hits_filtered)
    )
end

"""
    energy_filter_mt(cmdir::String, input_file::String;
                     emin::Float64=2400.0, emax::Float64=2500.0,
                     nthreads::Int=1)

Main function: filter events by energy and write to HDF5 files.
"""
function energy_filter_mt(cmdir::String, input_file::String;
                          emin::Float64=2400.0, emax::Float64=2500.0,
                          nthreads::Int=1)

    println("\n" * "=" ^ 70)
    println("ITACA ENERGY FILTER (MT)")
    println("=" ^ 70)

    # Load data
    input_path = joinpath(cmdir, input_file)
    println("Loading data from: $input_path")

    if !isfile(input_path)
        error("File not found: $input_path")
    end

    dfs = Petit.get_dataset_dfs(input_path)
    hitsdf = dfs["hits"]

    # Compute event energies
    println("Computing event energies...")
    event_energies = compute_event_energies(hitsdf)
    n_total_events = nrow(event_energies)
    println("Total events: $n_total_events")

    # Energy statistics (for summary)
    energy_stats = (
        mean = mean(event_energies.total_energy_keV),
        std = std(event_energies.total_energy_keV),
        min = minimum(event_energies.total_energy_keV),
        max = maximum(event_energies.total_energy_keV)
    )

    # Filter events
    println("\nFiltering events with energy in [$emin, $emax] keV...")
    filtered_event_ids = filter_events_by_energy(event_energies, emin, emax)
    n_filtered_events = length(filtered_event_ids)

    println("Events passing cut: $n_filtered_events / $n_total_events ($(round(100 * n_filtered_events / n_total_events, digits=2))%)")

    if n_filtered_events == 0
        println("\nWarning: No events pass the energy cut!")
        println("Energy range in data: [$(round(energy_stats.min, digits=2)), $(round(energy_stats.max, digits=2))] keV")
        return nothing
    end

    # Get base name from input file (without extension)
    base_name = splitext(input_file)[1]

    # Create output directory
    output_dir = create_output_dir(cmdir, base_name, emin, emax)

    # Determine number of threads
    optimal_nthreads = Petit.get_optimal_threads(nthreads)
    if optimal_nthreads > n_filtered_events
        optimal_nthreads = n_filtered_events
        println("Reducing threads to $optimal_nthreads (one per event)")
    end

    # Split events among threads
    thread_ranges = Petit.split_events_for_threads(n_filtered_events, optimal_nthreads, 1)

    println("\nConfiguration:")
    println("  Energy cut:    [$emin, $emax] keV")
    println("  Events to write: $n_filtered_events")
    println("  Threads:       $optimal_nthreads")
    println("  Output:        $output_dir")

    # Prepare metadata
    metadata = Dict(
        "source_file" => input_file,
        "emin_keV" => emin,
        "emax_keV" => emax,
        "total_events_in_source" => n_total_events,
        "total_events_filtered" => n_filtered_events,
        "filter_efficiency_pct" => round(100 * n_filtered_events / n_total_events, digits=2),
        "timestamp" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    )

    # Process in parallel
    println("\nWriting filtered data...")
    results = Vector{Dict}(undef, length(thread_ranges))

    Threads.@threads for i in 1:length(thread_ranges)
        start_idx, num_events = thread_ranges[i]
        end_idx = start_idx + num_events - 1

        # Get event IDs for this thread
        thread_event_ids = filtered_event_ids[start_idx:end_idx]

        results[i] = filter_and_write_thread(hitsdf, thread_event_ids, i, output_dir, base_name, metadata)

        println("  Thread $i: wrote $(results[i]["n_events"]) events to $(results[i]["output_file"])")
    end

    # Write summary file
    summary_path = write_summary_file(output_dir, input_file, emin, emax,
                                       n_total_events, n_filtered_events,
                                       energy_stats, results)

    # Print final summary
    println("\n" * "=" ^ 70)
    println("FILTER COMPLETE")
    println("=" ^ 70)
    println("Total events read:    $n_total_events")
    println("Events passing cut:   $n_filtered_events ($(round(100 * n_filtered_events / n_total_events, digits=2))%)")
    println("Output directory:     $output_dir")
    println("Summary file:         $summary_path")
    println("=" ^ 70)

    return Dict(
        "output_dir" => output_dir,
        "n_total_events" => n_total_events,
        "n_filtered_events" => n_filtered_events,
        "thread_results" => results,
        "summary_path" => summary_path
    )
end

#=============================================================================
# Command Line Interface
=============================================================================#

function main()
    # Check minimum required arguments
    if length(ARGS) < 2
        println("Error: Missing required arguments")
        println()
        println("Usage: julia -t <nthreads> itaca_energy_filter_mt.jl <cmdir> <input_file> [options]")
        println()
        println("Required arguments:")
        println("  cmdir           Directory containing the input file")
        println("  input_file      Name of the HDF5 input file")
        println()
        println("Optional arguments:")
        println("  --emin=X        Minimum energy in keV (default: 2400)")
        println("  --emax=X        Maximum energy in keV (default: 2500)")
        println("  --nthreads=N    Number of threads (default: 1)")
        println()
        println("Output:")
        println("  Creates directory: <cmdir>/<basename>_filtered_ecut_<emin>_<emax>/")
        println("  Contains: <basename>_part_1.h5, <basename>_part_2.h5, ..., filter_summary.txt")
        println()
        println("Examples:")
        println("  # Default energy cut [2400, 2500] keV:")
        println("  julia -t 4 itaca_energy_filter_mt.jl /data/itaca/ bb0nu.h5")
        println()
        println("  # Custom energy cut:")
        println("  julia -t 4 itaca_energy_filter_mt.jl /data/itaca/ bb0nu.h5 --emin=2450 --emax=2480")
        exit(1)
    end

    # Parse required arguments
    cmdir = ARGS[1]
    input_file = ARGS[2]

    # Default values
    emin = 2400.0
    emax = 2500.0
    nthreads = Threads.nthreads()  # Use Julia's thread count from -t flag

    # Parse optional arguments
    for arg in ARGS[3:end]
        if startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) != 2
                println("Warning: Ignoring malformed argument: $arg")
                continue
            end
            key, value = parts

            try
                if key == "emin"
                    emin = parse(Float64, value)
                elseif key == "emax"
                    emax = parse(Float64, value)
                elseif key == "nthreads"
                    nthreads = parse(Int, value)
                else
                    println("Warning: Unknown argument: --$key")
                end
            catch e
                println("Error parsing argument --$key=$value: $e")
                exit(1)
            end
        end
    end

    # Validate
    if emin >= emax
        println("Error: emin ($emin) must be less than emax ($emax)")
        exit(1)
    end

    # Run
    energy_filter_mt(cmdir, input_file;
                     emin=emin, emax=emax, nthreads=nthreads)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
