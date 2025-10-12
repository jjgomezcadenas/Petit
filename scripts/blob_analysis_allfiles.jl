#!/usr/bin/env julia

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

# Now load packages
using ArgParse
using CSV
using DataFrames
using Graphs
using HDF5
using JSON
using Statistics
using Printf

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
using .Petit

"""
    parse_commandline()

Parse command line arguments for multi-file blob analysis script.
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--file-type", "-t"
            help = "Type of files to process: bi214, tl208, 0nubb, or electron_0nubb"
            arg_type = String
            required = true
        "--data-dir", "-d"
            help = "Root data directory containing HDF5 files"
            arg_type = String
            default = "/Users/jjgomezcadenas/Data/HD5t/"
        "--events-per-file", "-n"
            help = "Number of events to process per file (0 for all)"
            arg_type = Int
            default = 0
        "--max-files"
            help = "Maximum number of files to process (0 for all)"
            arg_type = Int
            default = 0
        "--voxel-size"
            help = "Voxel size in mm"
            arg_type = Float64
            default = 3.0
        "--max-distance"
            help = "Maximum distance for track building in mm"
            arg_type = Float64
            default = 10.0
        "--energy-threshold"
            help = "Energy threshold in keV"
            arg_type = Float64
            default = 1.0
        "--xycut"
            help = "Fiducial cut for x,y coordinates"
            arg_type = Float64
            default = 1800.0
        "--zcut"
            help = "Fiducial cut for z coordinate"
            arg_type = Float64
            default = 100.0
        "--blob-radius"
            help = "Radius for blob energy calculation in mm"
            arg_type = Float64
            default = 15.0
        "--eblob-cut"
            help = "Energy cut for blob 2 in keV (tracks pass if blob2 >= cut)"
            arg_type = Float64
            default = 600.0
        "--output-json"
            help = "Output JSON file for tracks"
            arg_type = String
            default = ""  # Will be set based on file-type if empty
        "--output-csv"
            help = "Output CSV file for blob analysis"
            arg_type = String
            default = ""  # Will be set based on file-type if empty
        "--verbose", "-v"
            help = "Verbose output"
            action = :store_true
    end

    return parse_args(s)
end

"""
    find_matching_files(data_dir::String, file_type::String; verbose=false)

Find all HDF5 files in the data directory that match the given file type pattern.
"""
function find_matching_files(data_dir::String, file_type::String; verbose=false)
    # Validate file type
    valid_types = ["bi214", "tl208", "0nubb", "electron_0nubb"]
    if !(file_type in valid_types)
        error("Invalid file type: $file_type. Must be one of: $(join(valid_types, ", "))")
    end

    # Check if directory exists
    if !isdir(data_dir)
        error("Data directory does not exist: $data_dir")
    end

    # Use readdir to get all files, then filter for .h5 files
    all_files = readdir(data_dir, join=true)

    # Filter for .h5 files that contain the file type string
    matching_files = filter(all_files) do f
        isfile(f) &&
        endswith(lowercase(f), ".h5") &&
        occursin(file_type, lowercase(basename(f)))
    end

    if verbose
        println("Found $(length(matching_files)) files matching pattern '$file_type':")
        for (i, f) in enumerate(matching_files)
            println("  $i. $(basename(f))")
        end
    end

    emin = 2400.0
    emax = 2700.0
    if file_type == "bi214"
        emin = 2400.0
        emax = 2500.0
    end

    return matching_files, emin, emax
end

"""
    process_single_file(filepath::String, args::Dict; verbose=false)

Process a single HDF5 file and return the tracks found.
"""
function process_single_file(filepath::String, args::Dict; verbose=false)
    if verbose
        println("\n  Processing: $(basename(filepath))")
    end

    # Determine events to run
    events_to_run = args["events-per-file"]
    if events_to_run == 0
        # If 0, process all events (set to a large number)
        events_to_run = 1000000  # Will process all available events
    end

    try
        # Call event_loop_single_track for this file
        tracks = Petit.event_loop_single_track(
            args["data-dir"];
            input_file = filepath,
            events_to_run = events_to_run,
            voxel_size_mm = args["voxel-size"],
            max_distance_mm = args["max-distance"],
            energy_threshold_kev = args["energy-threshold"],
            xyc = args["xycut"],
            zc = args["zcut"]
        )

        if verbose
            println("    Found $(length(tracks)) single-track events")
        end

        return tracks
    catch e
        println("  WARNING: Failed to process $(basename(filepath)): $e")
        return []
    end
end

"""
    serialize_track_for_json(track::Tracks, walk_result, blob_result)

Convert track data to a serializable format for JSON output.
"""
function serialize_track_for_json(track::Tracks, walk_result, blob_result)
    # Extract voxel information
    voxels = Dict(
        "x" => track.voxels.x,
        "y" => track.voxels.y,
        "z" => track.voxels.z,
        "energy" => track.voxels.energy,
        "event_id" => track.voxels.event_id
    )

    # Extract graph information
    edge_list = []
    for edge in edges(track.graph)
        push!(edge_list, Dict("from" => src(edge), "to" => dst(edge)))
    end

    # Extract walk result
    extremes = if !isnothing(walk_result.extremes[1])
        Dict(
            "start" => Dict(
                "x" => walk_result.extremes[1].x,
                "y" => walk_result.extremes[1].y,
                "z" => walk_result.extremes[1].z,
                "energy" => walk_result.extremes[1].energy
            ),
            "end" => Dict(
                "x" => walk_result.extremes[2].x,
                "y" => walk_result.extremes[2].y,
                "z" => walk_result.extremes[2].z,
                "energy" => walk_result.extremes[2].energy
            )
        )
    else
        nothing
    end

    track_info = Dict{String,Any}(
        "voxels" => voxels,
        "graph" => Dict(
            "num_vertices" => nv(track.graph),
            "edges" => edge_list
        ),
        "walk_result" => Dict(
            "extremes" => extremes,
            "path_indices" => walk_result.path_indices,
            "total_length" => walk_result.total_length,
            "confidence" => walk_result.confidence
        ),
        "blob_result" => Dict(
            "blob1_energy" => blob_result.blob1_energy,
            "blob2_energy" => blob_result.blob2_energy,
            "blob1_voxel_count" => blob_result.blob1_voxel_count,
            "blob2_voxel_count" => blob_result.blob2_voxel_count,
            "blob1_center" => blob_result.blob1_center,
            "blob2_center" => blob_result.blob2_center
        )
    )

    return track_info
end

"""
    main()

Main function for multi-file blob analysis.
"""
function main()
    # Parse command line arguments
    args = parse_commandline()

    # Set default output filenames based on file type if not provided
    if isempty(args["output-json"])
        args["output-json"] = "tracks_$(args["file-type"]).json"
    end
    if isempty(args["output-csv"])
        args["output-csv"] = "blob_analysis_$(args["file-type"]).csv"
    end

    println("="^70)
    println("Multi-File Blob Analysis Script")
    println("="^70)
    println("File type: $(args["file-type"])")
    println("Data directory: $(args["data-dir"])")
    println("Events per file: $(args["events-per-file"] == 0 ? "all" : args["events-per-file"])")
    println("Max files: $(args["max-files"] == 0 ? "all" : args["max-files"])")
    println("Voxel size: $(args["voxel-size"]) mm")
    println("Max distance: $(args["max-distance"]) mm")
    println("Energy threshold: $(args["energy-threshold"]) keV")
    println("Fiducial cuts: xy < $(args["xycut"]) mm, z < $(args["zcut"]) mm")
    println("Blob radius: $(args["blob-radius"]) mm")
    println("Blob 2 energy cut: $(args["eblob-cut"]) keV")
    println("Output JSON: $(args["output-json"])")
    println("Output CSV: $(args["output-csv"])")
    println("="^70)

    # Step 1: Find matching files
    println("\n[1/4] Searching for files...")
    matching_files, emin, emax = find_matching_files(args["data-dir"], args["file-type"];
                                        verbose=args["verbose"])

    if isempty(matching_files)
        error("No files found matching pattern '$(args["file-type"])' in $(args["data-dir"])")
    end

    println("Found $(length(matching_files)) matching files")

    # Limit number of files if requested
    if args["max-files"] > 0 && length(matching_files) > args["max-files"]
        matching_files, emin, emax = matching_files[1:args["max-files"]]
        println("Processing only first $(args["max-files"]) files as requested")
    end

    # Step 2: Process each file and collect tracks
    println("\n[2/4] Processing files...")
    all_tracks = []
    files_processed = 0
    tracks_per_file = Dict{String, Int}()

    for (idx, filepath) in enumerate(matching_files)
        print("\r[$idx/$(length(matching_files))] Processing: $(basename(filepath))...")

        file_tracks = process_single_file(filepath, args; verbose=false)

        if !isempty(file_tracks)
            append!(all_tracks, file_tracks)
            files_processed += 1
            tracks_per_file[basename(filepath)] = length(file_tracks)
        end
    end

    println("\n\nProcessed $files_processed files successfully")
    println("Total single-track events collected: $(length(all_tracks))")

    if args["verbose"] && !isempty(tracks_per_file)
        println("\nTracks per file:")
        for (fname, ntrks) in sort(collect(tracks_per_file))
            println("  $fname: $ntrks tracks")
        end
    end

    if isempty(all_tracks)
        println("\nNo tracks found in any file. Exiting.")
        return
    end

    # Step 3: Analyze each track
    println("\n[3/4] Analyzing $(length(all_tracks)) tracks...")

    # Initialize arrays for CSV output
    source_files = String[]
    event_numbers = Int[]
    eblob1_values = Float64[]
    eblob2_values = Float64[]
    nblob1_values = Int[]
    nblob2_values = Int[]
    confidence_values = Float64[]
    etrk_values = Float64[]
    trkl_values = Float64[]
    select_values = Bool[]

    # Initialize array for JSON output
    tracks_json = []

    # Get the energy cut value
    eblob_cut = args["eblob-cut"]

    # Track statistics for filtering
    n_short_tracks = 0
    n_high_blob1 = 0
    n_high_blob2 = 0
    n_ecut = 0

    for (i, track) in enumerate(all_tracks)
        # Get event number from the first voxel
        event_id = track.voxels.event_id[1]

        # Walk track from extremes
        walk_result = Petit.walk_track_from_extremes(track)

        # Apply track length filter
        if walk_result.total_length < 50.0
            n_short_tracks += 1
            continue
        end

        # Calculate track energy in keV
        track_energy_kev = 1e3 * sum(track.voxels.energy)

        if track_energy_kev < emin || track_energy_kev > emax
            n_ecut += 1
            continue
        end

        # Calculate energy in spheres around extremes
        blob_result = Petit.energy_in_spheres_around_extremes(
            track,
            walk_result,
            args["blob-radius"]
        )

        # Convert blob energies to keV
        blob1_energy_kev = 1e3 * blob_result.blob1_energy
        blob2_energy_kev = 1e3 * blob_result.blob2_energy

        # Apply blob energy filters
        if blob2_energy_kev > 1000.0
            n_high_blob2 += 1
            continue
        end

        if blob1_energy_kev > 1500.0
            n_high_blob1 += 1
            continue
        end

        # Determine if track passes the energy cut
        select = blob2_energy_kev >= eblob_cut

        # Store values for CSV
        push!(source_files, args["file-type"])  # Use file type as source identifier
        push!(event_numbers, event_id)
        push!(eblob1_values, blob1_energy_kev)
        push!(eblob2_values, blob2_energy_kev)
        push!(nblob1_values, blob_result.blob1_voxel_count)
        push!(nblob2_values, blob_result.blob2_voxel_count)
        push!(confidence_values, walk_result.confidence)
        push!(etrk_values, track_energy_kev)
        push!(trkl_values, walk_result.total_length)
        push!(select_values, select)

        # Store track for JSON
        track_json = serialize_track_for_json(track, walk_result, blob_result)
        track_json["event_id"] = event_id
        track_json["source"] = args["file-type"]
        push!(tracks_json, track_json)

        # Progress report
        if i % 100 == 0 || i == length(all_tracks)
            print("\rProcessed: $i/$(length(all_tracks)) tracks")
        end
    end
    println()

    # Print filtering statistics
    if args["verbose"]
        println("\nFiltering statistics:")
         println("  Removed due to ecut between $emin and $emax: $n_ecut")
        println("  Removed due to track length < 50 mm: $n_short_tracks")
        println("  Removed due to blob1 energy > 1500 keV: $n_high_blob1")
        println("  Removed due to blob2 energy > 1000 keV: $n_high_blob2")
        println("  Tracks passing all filters: $(length(etrk_values))")
    end

    # Step 4: Write output files
    println("\n[4/4] Writing output files...")

    # Write JSON file
    open(args["output-json"], "w") do io
        JSON.print(io, tracks_json, 2)
    end
    println("Written track data to: $(args["output-json"])")

    # Create and write CSV DataFrame
    results_df = DataFrame(
        source = source_files,
        event = event_numbers,
        eblob1 = round.(eblob1_values, digits=2),
        eblob2 = round.(eblob2_values, digits=2),
        nblob1 = nblob1_values,
        nblob2 = nblob2_values,
        confidence = round.(confidence_values, digits=4),
        etrk = round.(etrk_values, digits=2),
        trkl = round.(trkl_values, digits=2),
        select = select_values
    )

    CSV.write(args["output-csv"], results_df)
    println("Written blob analysis to: $(args["output-csv"])")

    # Print summary statistics
    if !isempty(etrk_values)
        println("\n" * "="^70)
        println("Summary Statistics:")
        println("="^70)
        println("Files processed: $files_processed")
        println("Total tracks collected: $(length(all_tracks))")
        println("Tracks analyzed (after filters): $(length(etrk_values))")
        println("Tracks passing cut (blob2 >= $(args["eblob-cut"]) keV): $(sum(select_values)) ($(round(100*sum(select_values)/length(select_values), digits=1))%)")
        println("\nTrack properties:")
        println("  Average track energy: $(round(mean(etrk_values), digits=1)) ± $(round(std(etrk_values), digits=1)) keV")
        println("  Average track length: $(round(mean(trkl_values), digits=1)) ± $(round(std(trkl_values), digits=1)) mm")
        println("  Average blob1 energy: $(round(mean(eblob1_values), digits=1)) ± $(round(std(eblob1_values), digits=1)) keV")
        println("  Average blob2 energy: $(round(mean(eblob2_values), digits=1)) ± $(round(std(eblob2_values), digits=1)) keV")
        println("  Average confidence: $(round(mean(confidence_values), digits=3)) ± $(round(std(confidence_values), digits=3))")
        println("="^70)
    else
        println("\nNo tracks survived the filtering criteria.")
    end
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end