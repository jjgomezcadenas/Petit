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

Parse command line arguments for blob analysis script.
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input", "-i"
            help = "Input HDF5 file path"
            arg_type = String
            default = "/Users/jjgomezcadenas/Data/HD5t/0nubb.h5"
        "--events", "-n"
            help = "Number of events to process"
            arg_type = Int
            default = 10
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
            default = "tracks.json"
        "--output-csv"
            help = "Output CSV file for blob analysis"
            arg_type = String
            default = "blob_analysis.csv"
        "--data-dir"
            help = "Data directory path"
            arg_type = String
            default = joinpath(ENV["DATA"], "HD5t")
    end

    return parse_args(s)
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

Main function for blob analysis.
"""
function main()
    # Parse command line arguments
    args = parse_commandline()

    println("="^60)
    println("Blob Analysis Script")
    println("="^60)
    println("Input file: $(args["input"])")
    println("Events to process: $(args["events"])")
    println("Voxel size: $(args["voxel-size"]) mm")
    println("Max distance: $(args["max-distance"]) mm")
    println("Energy threshold: $(args["energy-threshold"]) keV")
    println("Fiducial cuts: xy < $(args["xycut"]) mm, z < $(args["zcut"]) mm")
    println("Blob radius: $(args["blob-radius"]) mm")
    println("Blob 2 energy cut: $(args["eblob-cut"]) keV (select = true if blob2 >= cut)")
    println("="^60)

    # Step 1: Run event_loop_single_track to get tracks
    println("\n[1/3] Processing tracks...")
    tracks = Petit.event_loop_single_track(
        args["data-dir"];
        input_file = args["input"],
        events_to_run = args["events"],
        voxel_size_mm = args["voxel-size"],
        max_distance_mm = args["max-distance"],
        energy_threshold_kev = args["energy-threshold"],
        xyc = args["xycut"],
        zc = args["zcut"]
    )

    println("Found $(length(tracks)) single-track events")

    # Step 2: Analyze each track
    println("\n[2/3] Analyzing tracks...")

    # Initialize arrays for CSV output
    event_numbers = Int[]
    eblob1_values = Float64[]
    eblob2_values = Float64[]
    nblob1_values = Int[]
    nblob2_values = Int[]
    confidence_values = Float64[]
    etrk_values = Float64[]
    trkl_values = Float64[]
    select_values = Bool[]  # Track whether each event passes the cut

    # Initialize array for JSON output
    tracks_json = []

    # Get the energy cut value
    eblob_cut = args["eblob-cut"]

    for (i, track) in enumerate(tracks)
        # Get event number from the first voxel
        event_id = track.voxels.event_id[1]

        # Walk track from extremes
        walk_result = Petit.walk_track_from_extremes(track)

        # Calculate energy in spheres around extremes
        blob_result = Petit.energy_in_spheres_around_extremes(
            track,
            walk_result,
            args["blob-radius"]
        )

        # Calculate track energy in keV
        track_energy_kev = 1e3 * sum(track.voxels.energy)

        if walk_result.total_length < 50.0
            continue
        end
        
        

        # Convert blob energies to keV
        blob1_energy_kev = 1e3 * blob_result.blob1_energy
        blob2_energy_kev = 1e3 * blob_result.blob2_energy

        if blob2_energy_kev > 1000.0
            continue
        end

        if blob1_energy_kev > 1500.0
            continue
        end

        # Determine if track passes the energy cut
        select = blob2_energy_kev >= eblob_cut

        # Store values for CSV
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
        push!(tracks_json, track_json)

        # Progress report
        if i % 10 == 0 || i == length(tracks)
            print("\rProcessed: $i/$(length(tracks)) tracks")
        end
    end
    println()

    # Step 3: Write output files
    println("\n[3/3] Writing output files...")

    # Write JSON file
    open(args["output-json"], "w") do io
        JSON.print(io, tracks_json, 2)
    end
    println("Written track data to: $(args["output-json"])")

    # Create and write CSV DataFrame
    results_df = DataFrame(
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
    println("\n" * "="^60)
    println("Summary Statistics:")
    println("="^60)
    println("Total tracks analyzed: $(length(etrk_values))")
    println("Tracks passing cut (blob2 >= $(args["eblob-cut"]) keV): $(sum(select_values)) ($(round(100*sum(select_values)/length(select_values), digits=1))%)")
    println("Average track energy: $(round(mean(etrk_values), digits=1)) keV")
    println("Average track length: $(round(mean(trkl_values), digits=1)) mm")
    println("Average blob1 energy: $(round(mean(eblob1_values), digits=1)) keV")
    println("Average blob2 energy: $(round(mean(eblob2_values), digits=1)) keV")
    println("Average confidence: $(round(mean(confidence_values), digits=3))")
    println("="^60)
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end