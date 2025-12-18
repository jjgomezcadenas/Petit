#!/usr/bin/env julia

"""
Event selection script based on blob energy cut.

Reads reconstructed tracks (single file or multiple files via pattern),
computes blob energies at given radius, and saves events passing the
Eb2 cut to a new HDF5 file.

Eb1 is the higher energy blob, Eb2 is the lower energy blob.
Events with Eb2 >= eb2cut are selected (signal-like).

The output HDF5 file has the same format as the input.

Usage:
  # Single file:
  julia select_events_rb.jl /path/to/file.h5 --rb=12 --eb2cut=250 --outfile=selected.h5

  # Multiple files with pattern:
  julia select_events_rb.jl /path/to/dir "pattern*.h5" --rb=12 --eb2cut=250 --outfile=selected.h5
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using Printf
using HDF5
using DataFrames
using Graphs
using Glob

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

"""
    save_central_path_to_hdf5(central_path, group)

Save a central path DataFrame to an HDF5 group.
"""
function save_central_path_to_hdf5(central_path::DataFrame, group)
    if nrow(central_path) == 0
        group["central_path_data"] = zeros(Float64, 0, 0)
        group["central_path_columns"] = String[]
        return
    end

    # Save as matrix with column names
    cp_data = Matrix(central_path)
    group["central_path_data"] = cp_data
    group["central_path_columns"] = String.(names(central_path))
end

"""
    save_selected_results_to_hdf5(selected_results, output_path, metadata)

Save selected reconstruction results to HDF5 file in the same format as input.
"""
function save_selected_results_to_hdf5(selected_results::Vector, output_path::String, metadata::Dict)
    h5open(output_path, "w") do fid
        # Update total_tracks_saved in metadata before saving
        metadata["total_tracks_saved"] = length(selected_results)

        # Save metadata as attributes
        for (key, val) in metadata
            HDF5.attributes(fid)[key] = val
        end

        if isempty(selected_results)
            return
        end

        println("    Saving $(length(selected_results)) selected tracks to HDF5...")
        flush(stdout)

        for (idx, sel) in enumerate(selected_results)
            if idx % 10 == 0 || idx == length(selected_results)
                print("\r    Saving track $idx/$(length(selected_results)) to HDF5...")
                flush(stdout)
            end

            r = sel.result
            track = r.track

            track_group_name = "batch_1/track_$(idx)"
            g = create_group(fid, track_group_name)

            # Save event metadata
            g["event_id"] = r.event_id
            g["track_length"] = r.track_length
            g["confidence"] = r.confidence

            # Save MC truth extremes (if available)
            if !isnothing(r.mc_pos1) && !isnothing(r.mc_pos2)
                g["mc_pos1"] = collect(r.mc_pos1)
                g["mc_pos2"] = collect(r.mc_pos2)
            end

            # Save voxels data
            voxels_data = Matrix(track.voxels)
            g["voxels"] = voxels_data
            g["voxel_columns"] = String.(names(track.voxels))

            # Save graph edges
            n_edges = ne(track.graph)
            if n_edges > 0
                edge_matrix = zeros(Int, n_edges, 2)
                for (i, edge) in enumerate(edges(track.graph))
                    edge_matrix[i, 1] = src(edge)
                    edge_matrix[i, 2] = dst(edge)
                end
                g["graph_edges"] = edge_matrix
            else
                g["graph_edges"] = zeros(Int, 0, 2)
            end
            g["n_vertices"] = nv(track.graph)

            # Save components
            if !isempty(track.components)
                max_comp_len = maximum(length.(track.components))
                comp_matrix = zeros(Int, length(track.components), max_comp_len)
                for (i, comp) in enumerate(track.components)
                    for (j, v) in enumerate(comp)
                        comp_matrix[i, j] = v
                    end
                end
                g["components"] = comp_matrix
            else
                g["components"] = zeros(Int, 0, 0)
            end

            # Save central path
            save_central_path_to_hdf5(r.central_path, g)
        end
        println()  # New line after progress
    end
end

"""
    select_events_from_chain(input_dir, pattern, output_file; rb, eb2cut)

Select events from multiple files matching a pattern.

# Arguments
- `input_dir`: Directory containing input files
- `pattern`: Glob pattern to match files (e.g., "electrons_*.h5")
- `output_file`: Name for output file (saved in input_dir)
- `rb`: Blob sphere radius in mm
- `eb2cut`: Minimum Eb2 energy in keV for selection
"""
function select_events_from_chain(input_dir::String, pattern::String, output_file::String;
                                   rb::Float64=12.0,
                                   eb2cut::Float64=200.0)

    println("="^70)
    println("EVENT SELECTION FROM MULTIPLE FILES")
    println("="^70)
    println("Input directory: $input_dir")
    println("Pattern: $pattern")
    println("Output file: $output_file")
    println("Blob radius: $rb mm")
    println("Eb2 cut:     $eb2cut keV")

    # Find matching files
    full_pattern = joinpath(input_dir, pattern)
    files = glob(pattern, input_dir)
    sort!(files)

    if isempty(files)
        error("No files found matching pattern: $full_pattern")
    end

    println("\nFound $(length(files)) files:")
    for f in files
        println("  - $(basename(f))")
    end

    # Chain all results
    println("\nLoading reconstructed tracks from all files...")
    all_results, all_metadata = Petit.chain_reco_results(input_dir, pattern)
    println("  Loaded $(length(all_results)) total tracks")

    # Merge metadata from first file as base (they should have similar structure)
    metadata = if !isempty(all_metadata)
        copy(first(all_metadata))
    else
        Dict{String, Any}()
    end

    if isempty(all_results)
        error("No tracks found in input files")
    end

    # Select events based on Eb2 cut
    println("\nSelecting events with Eb2 >= $eb2cut keV...")
    selected_results = []
    n_failed = 0

    for r in all_results
        try
            blobs = Petit.find_blob_energies(r.track, r.central_path; radius=rb)
            # Convert to keV
            Eb1_keV = blobs.Eb1 * 1e3
            Eb2_keV = blobs.Eb2 * 1e3

            if Eb2_keV >= eb2cut
                push!(selected_results, (
                    result=r,
                    Eb1=Eb1_keV,
                    Eb2=Eb2_keV,
                    asymmetry=blobs.asymmetry
                ))
            end
        catch e
            n_failed += 1
            continue
        end
    end

    n_total = length(all_results)
    n_selected = length(selected_results)
    n_rejected = n_total - n_selected - n_failed

    println("\nSelection summary:")
    println("  Total tracks:    $n_total")
    println("  Failed analysis: $n_failed")
    println("  Rejected:        $n_rejected")
    println("  Selected:        $n_selected")
    println("  Efficiency:      $(round(n_selected/n_total*100, digits=1))%")

    if n_selected == 0
        println("\nWarning: No events passed the selection!")
        return
    end

    # Update metadata with selection info
    metadata["rb_mm"] = rb
    metadata["eb2cut_keV"] = eb2cut
    metadata["n_input_tracks"] = n_total
    metadata["n_selected_tracks"] = n_selected
    metadata["selection_efficiency"] = n_selected / n_total
    metadata["n_input_files"] = length(files)
    metadata["input_pattern"] = pattern

    # Save selected events
    output_path = joinpath(input_dir, output_file)
    println("\nSaving selected events to: $output_path")
    save_selected_results_to_hdf5(selected_results, output_path, metadata)

    println("Done! Saved $n_selected events to $output_path")

    return (n_total=n_total, n_selected=n_selected,
            efficiency=n_selected/n_total, output_path=output_path)
end

"""
    select_events(input_file, output_file; rb, eb2cut)

Select events based on blob energy cut from a single file.

# Arguments
- `input_file`: Path to input HDF5 file with reconstructed tracks
- `output_file`: Name for output file (saved in same directory as input)
- `rb`: Blob sphere radius in mm
- `eb2cut`: Minimum Eb2 energy in keV for selection
"""
function select_events(input_file::String, output_file::String;
                       rb::Float64=12.0,
                       eb2cut::Float64=200.0)

    println("="^70)
    println("EVENT SELECTION BASED ON BLOB ENERGY")
    println("="^70)
    println("Input file:  $input_file")
    println("Output file: $output_file")
    println("Blob radius: $rb mm")
    println("Eb2 cut:     $eb2cut keV")

    # Check input file exists
    if !isfile(input_file)
        error("Input file not found: $input_file")
    end

    # Determine output path (same directory as input)
    input_dir = dirname(input_file)
    output_path = joinpath(input_dir, output_file)

    # Read reconstructed tracks
    println("\nLoading reconstructed tracks...")
    results, metadata = Petit.read_reco_results_from_hdf5(input_file)
    println("  Loaded $(length(results)) tracks")

    if isempty(results)
        error("No tracks found in input file")
    end

    # Select events based on Eb2 cut
    println("\nSelecting events with Eb2 >= $eb2cut keV...")
    selected_results = []
    n_failed = 0

    for r in results
        try
            blobs = Petit.find_blob_energies(r.track, r.central_path; radius=rb)
            # Convert to keV
            Eb1_keV = blobs.Eb1 * 1e3
            Eb2_keV = blobs.Eb2 * 1e3

            if Eb2_keV >= eb2cut
                push!(selected_results, (
                    result=r,
                    Eb1=Eb1_keV,
                    Eb2=Eb2_keV,
                    asymmetry=blobs.asymmetry
                ))
            end
        catch e
            n_failed += 1
            continue
        end
    end

    n_total = length(results)
    n_selected = length(selected_results)
    n_rejected = n_total - n_selected - n_failed

    println("\nSelection summary:")
    println("  Total tracks:    $n_total")
    println("  Failed analysis: $n_failed")
    println("  Rejected:        $n_rejected")
    println("  Selected:        $n_selected")
    println("  Efficiency:      $(round(n_selected/n_total*100, digits=1))%")

    if n_selected == 0
        println("\nWarning: No events passed the selection!")
        return
    end

    # Update metadata with selection info
    metadata["rb_mm"] = rb
    metadata["eb2cut_keV"] = eb2cut
    metadata["n_input_tracks"] = n_total
    metadata["n_selected_tracks"] = n_selected
    metadata["selection_efficiency"] = n_selected / n_total

    # Save selected events to new HDF5 file
    println("\nSaving selected events to: $output_path")
    save_selected_results_to_hdf5(selected_results, output_path, metadata)

    println("Done! Saved $n_selected events to $output_path")

    return (n_total=n_total, n_selected=n_selected,
            efficiency=n_selected/n_total, output_path=output_path)
end

"""
    main()

Main function for command-line usage.
"""
function main()
    # Default parameters
    input_path = ""
    pattern = ""
    output_file = "selected_events.h5"
    rb = 12.0
    eb2cut = 200.0

    # Parse command line arguments
    positional_args = String[]
    for arg in ARGS
        if startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) != 2
                println("Warning: Ignoring malformed argument: $arg")
                continue
            end
            key, value = parts

            try
                if key == "rb"
                    rb = parse(Float64, value)
                elseif key == "eb2cut"
                    eb2cut = parse(Float64, value)
                elseif key == "outfile"
                    output_file = String(value)
                else
                    println("Warning: Unknown argument: --$key")
                end
            catch e
                println("Error parsing argument --$key=$value: $e")
                exit(1)
            end
        else
            push!(positional_args, arg)
        end
    end

    # Determine mode: single file or pattern
    if length(positional_args) == 0
        println("Usage:")
        println("  Single file:")
        println("    julia select_events_rb.jl <input_file> [options]")
        println("")
        println("  Multiple files (pattern):")
        println("    julia select_events_rb.jl <input_dir> <pattern> [options]")
        println("")
        println("Required arguments:")
        println("  input_file       Path to HDF5 file with reconstructed tracks")
        println("  OR")
        println("  input_dir        Directory containing input files")
        println("  pattern          Glob pattern to match files (e.g., \"electrons_*.h5\")")
        println("")
        println("Optional arguments:")
        println("  --rb=X           Blob sphere radius in mm (default: 12.0)")
        println("  --eb2cut=X       Minimum Eb2 energy in keV (default: 200.0)")
        println("  --outfile=X      Output filename (default: selected_events.h5)")
        println("")
        println("Examples:")
        println("  # Single file:")
        println("  julia select_events_rb.jl data/reco.h5 --rb=12 --eb2cut=250 --outfile=selected.h5")
        println("")
        println("  # Multiple files:")
        println("  julia select_events_rb.jl /path/to/data \"electrons_*_reco_*.h5\" --rb=5 --eb2cut=250 --outfile=electrons_selected.h5")
        exit(0)
    elseif length(positional_args) == 1
        # Single file mode
        input_path = positional_args[1]
        if !isfile(input_path)
            error("Input file not found: $input_path")
        end
        select_events(input_path, output_file; rb=rb, eb2cut=eb2cut)
    else
        # Pattern mode: directory + pattern
        input_dir = positional_args[1]
        pattern = positional_args[2]

        if !isdir(input_dir)
            error("Input directory not found: $input_dir")
        end

        select_events_from_chain(input_dir, pattern, output_file; rb=rb, eb2cut=eb2cut)
    end
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
