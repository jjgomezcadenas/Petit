#!/usr/bin/env julia

"""
Select events from H5 files for faster testing.

Usage:
    julia select_events.jl <input_file> --events 1,5,12,100 --output selected.h5
    julia select_events.jl <input_file> --range 1:100 --output selected.h5
    julia select_events.jl <input_file> --file event_list.txt --output selected.h5
    julia select_events.jl <input_file> --list  # Just list available events
"""

const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)

using HDF5
using DataFrames
using ArgParse
using Printf
using Statistics

include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

const DEFAULT_CMDIR = joinpath(ENV["DATA"], "HD5t/itaca")

"""
    get_available_events(h5file)

Get list of available event IDs from H5 file.
"""
function get_available_events(h5file::String)
    events = Int[]
    h5open(h5file, "r") do fid
        if haskey(fid, "MC/hits")
            hits = read(fid, "MC/hits")
            # hits is a Vector of NamedTuples
            event_ids = [h.event_id for h in hits]
            events = sort(unique(event_ids))
        end
    end
    return events
end

"""
    count_hits_per_event(h5file)

Count hits per event and return DataFrame with event_id and nhits.
"""
function count_hits_per_event(h5file::String)
    h5open(h5file, "r") do fid
        hits = read(fid, "MC/hits")
        # hits is a Vector of NamedTuples
        event_ids = [h.event_id for h in hits]

        # Count hits per event
        counts = Dict{Int, Int}()
        for eid in event_ids
            counts[eid] = get(counts, eid, 0) + 1
        end

        df = DataFrame(event_id = collect(keys(counts)), nhits = collect(values(counts)))
        sort!(df, :event_id)
        return df
    end
end

"""
    extract_events(input_file, output_file, event_ids)

Extract selected events from input H5 file to output H5 file.
Uses h5copy to preserve the exact compound dataset structure.
"""
function extract_events(input_file::String, output_file::String, event_ids::Vector{Int})
    println("Extracting $(length(event_ids)) events...")
    println("  Input: $input_file")
    println("  Output: $output_file")

    # Read all hits and filter
    hits_df = h5open(input_file, "r") do fin
        hits = read(fin, "MC/hits")
        event_set = Set(event_ids)
        selected = [h for h in hits if h.event_id in event_set]
        DataFrame(selected)
    end

    n_selected = nrow(hits_df)
    if n_selected == 0
        error("No hits found for selected events!")
    end

    println("  Selected $n_selected hits")

    # Write as simple columnar format that load_data can handle
    h5open(output_file, "w") do fout
        create_group(fout, "MC")

        # Write each column as separate dataset under MC/hits group
        # This is a workaround - we'll need to modify load_data to handle this
        hits_grp = create_group(fout, "MC/hits")
        hits_grp["event_id"] = Vector{Int64}(hits_df.event_id)
        hits_grp["x"] = Vector{Float32}(hits_df.x)
        hits_grp["y"] = Vector{Float32}(hits_df.y)
        hits_grp["z"] = Vector{Float32}(hits_df.z)
        hits_grp["time"] = Vector{Float32}(hits_df.time)
        hits_grp["energy"] = Vector{Float32}(hits_df.energy)
        hits_grp["label"] = Vector{String}(hits_df.label)
        hits_grp["particle_id"] = Vector{Int32}(hits_df.particle_id)
        hits_grp["hit_id"] = Vector{Int32}(hits_df.hit_id)

        # Store metadata
        attrs(fout)["selected_events"] = collect(event_ids)
        attrs(fout)["source_file"] = input_file
        attrs(fout)["n_events_selected"] = length(event_ids)
        attrs(fout)["format"] = "columnar"  # Mark as columnar format
    end

    println("  Done!")
end


"""
    load_selected_data(filename)

Load data from a selected events file (columnar format).
Returns a DataFrame compatible with Petit functions.
"""
function load_selected_data(filename::String)
    h5open(filename, "r") do fid
        hits_grp = fid["MC/hits"]
        DataFrame(
            event_id = read(hits_grp["event_id"]),
            x = read(hits_grp["x"]),
            y = read(hits_grp["y"]),
            z = read(hits_grp["z"]),
            time = read(hits_grp["time"]),
            energy = read(hits_grp["energy"]),
            label = read(hits_grp["label"]),
            particle_id = read(hits_grp["particle_id"]),
            hit_id = read(hits_grp["hit_id"])
        )
    end
end

"""
    parse_event_spec(spec)

Parse event specification string into list of event IDs.
Supports: "1,5,12" or "1:10" or "1:2:10" (start:step:end)
"""
function parse_event_spec(spec::String)
    spec = strip(spec)

    if contains(spec, ",")
        # Comma-separated list
        return [parse(Int, s) for s in split(spec, ",")]
    elseif contains(spec, ":")
        # Range specification
        parts = split(spec, ":")
        if length(parts) == 2
            start, stop = parse(Int, parts[1]), parse(Int, parts[2])
            return collect(start:stop)
        elseif length(parts) == 3
            start, step, stop = parse(Int, parts[1]), parse(Int, parts[2]), parse(Int, parts[3])
            return collect(start:step:stop)
        else
            error("Invalid range specification: $spec")
        end
    else
        # Single event
        return [parse(Int, spec)]
    end
end

"""
    read_events_from_file(filename)

Read event IDs from a text file (one per line or comma-separated).
"""
function read_events_from_file(filename::String)
    events = Int[]
    for line in eachline(filename)
        line = strip(line)
        isempty(line) && continue
        startswith(line, "#") && continue  # Skip comments

        if contains(line, ",")
            append!(events, [parse(Int, s) for s in split(line, ",")])
        else
            push!(events, parse(Int, line))
        end
    end
    return events
end

function main()
    s = ArgParseSettings(description="Select events from H5 files for faster testing")

    @add_arg_table! s begin
        "input"
            help = "Input H5 file (relative to DATA/HD5t/itaca or absolute path)"
            required = true
        "--events", "-e"
            help = "Event IDs to select (comma-separated or range like 1:10)"
            arg_type = String
            default = ""
        "--file", "-f"
            help = "File containing event IDs (one per line)"
            arg_type = String
            default = ""
        "--output", "-o"
            help = "Output H5 file"
            arg_type = String
            default = ""
        "--list", "-l"
            help = "List available events and exit"
            action = :store_true
        "--stats"
            help = "Show hit count statistics per event"
            action = :store_true
        "--cmdir"
            help = "Base directory for input files"
            arg_type = String
            default = DEFAULT_CMDIR
    end

    args = parse_args(s)

    # Resolve input file path
    input_file = args["input"]
    if !isabspath(input_file)
        input_file = joinpath(args["cmdir"], input_file)
    end

    if !isfile(input_file)
        error("Input file not found: $input_file")
    end

    println("Input file: $input_file")

    # List mode
    if args["list"]
        events = get_available_events(input_file)
        println("\nAvailable events: $(length(events))")
        println("Event ID range: $(minimum(events)) - $(maximum(events))")
        println("\nFirst 20 events: $(events[1:min(20, length(events))])")
        if length(events) > 20
            println("Last 20 events: $(events[end-19:end])")
        end
        return
    end

    # Stats mode
    if args["stats"]
        println("\nComputing hit statistics...")
        df = count_hits_per_event(input_file)
        println("\nHit count statistics:")
        println("  Total events: $(nrow(df))")
        println("  Min hits: $(minimum(df.nhits))")
        println("  Max hits: $(maximum(df.nhits))")
        println("  Mean hits: $(round(mean(df.nhits), digits=1))")
        println("  Median hits: $(round(median(df.nhits), digits=1))")

        println("\nTop 10 events by hit count:")
        sorted_df = sort(df, :nhits, rev=true)
        for i in 1:min(10, nrow(sorted_df))
            println("  Event $(sorted_df.event_id[i]): $(sorted_df.nhits[i]) hits")
        end
        return
    end

    # Get event IDs to select
    event_ids = Int[]

    if !isempty(args["events"])
        event_ids = parse_event_spec(args["events"])
    elseif !isempty(args["file"])
        event_ids = read_events_from_file(args["file"])
    else
        error("Must specify --events or --file")
    end

    # Validate events exist
    available = Set(get_available_events(input_file))
    valid_events = [e for e in event_ids if e in available]
    invalid_events = [e for e in event_ids if !(e in available)]

    if !isempty(invalid_events)
        println("Warning: $(length(invalid_events)) events not found: $(invalid_events[1:min(5, length(invalid_events))])...")
    end

    if isempty(valid_events)
        error("No valid events to select!")
    end

    println("Selected $(length(valid_events)) valid events")

    # Output file
    output_file = args["output"]
    if isempty(output_file)
        base = basename(input_file)
        name, ext = splitext(base)
        output_file = "$(name)_selected$(ext)"
    end

    # Extract events
    extract_events(input_file, output_file, valid_events)

    println("\nOutput file: $output_file")
    println("Use with: --input $output_file")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
