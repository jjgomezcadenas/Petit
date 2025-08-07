#!/usr/bin/env julia

# Batch processing script for Petit HDF5 analysis
# Usage: julia run_hd5t.jl <input_dir> [options]

using Pkg
# Activate the project environment
Pkg.activate(dirname(@__DIR__))

using Petit
using Printf

# Parse command line arguments
function parse_args(args)
    if length(args) < 1
        println("Usage: julia run_hd5t.jl <input_dir> [options]")
        println("Options:")
        println("  --input-file=<file>      Input HDF5 file (default: 0nubb.next.h5)")
        println("  --events=<n>            Number of events to process (default: 100)")
        println("  --voxel-size=<mm>       Voxel size in mm (default: 5)")
        println("  --max-distance=<mm>     Max distance in mm (default: 10)")
        println("  --energy-threshold=<keV> Energy threshold in keV (default: 10)")
        println("  --output-dir=<dir>      Output directory (default: znubb)")
        exit(1)
    end
    
    # Default values
    params = Dict(
        :input_dir => args[1],
        :input_file => "0nubb.next.h5",
        :events_to_run => 100,
        :voxel_size_mm => 5.0,
        :max_distance_mm => 10.0,
        :energy_threshold_kev => 10.0,
        :output_dir => "znubb"
    )
    
    # Parse optional arguments
    for i in 2:length(args)
        arg = args[i]
        if startswith(arg, "--input-file=")
            params[:input_file] = split(arg, "=")[2]
        elseif startswith(arg, "--events=")
            params[:events_to_run] = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--voxel-size=")
            params[:voxel_size_mm] = parse(Float64, split(arg, "=")[2])
        elseif startswith(arg, "--max-distance=")
            params[:max_distance_mm] = parse(Float64, split(arg, "=")[2])
        elseif startswith(arg, "--energy-threshold=")
            params[:energy_threshold_kev] = parse(Float64, split(arg, "=")[2])
        elseif startswith(arg, "--output-dir=")
            params[:output_dir] = split(arg, "=")[2]
        else
            println("Warning: Unknown argument $arg")
        end
    end
    
    return params
end

function main()
    # Parse command line arguments
    params = parse_args(ARGS)
    
    # Create output directory if it doesn't exist
    if !isdir(params[:output_dir])
        mkpath(params[:output_dir])
        println("Created output directory: $(params[:output_dir])")
    end
    
    # Print run parameters
    println("\n=== Petit HDF5 Analysis Batch Job ===")
    println("Input directory: $(params[:input_dir])")
    println("Input file: $(params[:input_file])")
    println("Events to process: $(params[:events_to_run])")
    println("Voxel size: $(params[:voxel_size_mm]) mm")
    println("Max distance: $(params[:max_distance_mm]) mm")
    println("Energy threshold: $(params[:energy_threshold_kev]) keV")
    println("Output directory: $(params[:output_dir])")
    println("=====================================\n")
    
    # Run the event loop
    println("Starting analysis...")
    start_time = time()
    
    try
        cmdir = joinpath(ENV["DATA"], "HD5t")
        xfile = joinpath(cmdir, params[:input_dir])
        rsgn = event_loop(xfile; 
                         input_file=params[:input_file],
                         events_to_run=params[:events_to_run], 
                         voxel_size_mm=params[:voxel_size_mm],
                         max_distance_mm=params[:max_distance_mm], 
                         energy_threshold_kev=params[:energy_threshold_kev])
        
        println("\nAnalysis complete!")
        println("Events processed: $(rsgn.n_events_processed)")
        println("Single track events: $(rsgn.n_single_track)")
        println("Two track events: $(rsgn.n_two_track)")
        println("Three+ track events: $(rsgn.n_three_plus_track)")
        println("Failed events: $(rsgn.n_failed)")
        
        # Generate and save histograms for single track events
        println("\nGenerating histograms for single track events...")
        hx, hy, hz, he = histogram_results(rsgn.single_track)
        HSt1 = Dict("hx" => hx, "hy" => hy, "hz" => hz, "he" => he)
        output_file = joinpath(params[:output_dir], "HSt1.txt")
        save_histos(HSt1, output_file)
        println("Saved: $output_file")
        
        # Generate and save histograms for two track primary
        println("Generating histograms for two track primary...")
        hx, hy, hz, he = histogram_results(rsgn.two_track_primary)
        HSt2p = Dict("hx" => hx, "hy" => hy, "hz" => hz, "he" => he)
        output_file = joinpath(params[:output_dir], "HSt2p.txt")
        save_histos(HSt2p, output_file)
        println("Saved: $output_file")
        
        # Generate and save histograms for two track secondary
        println("Generating histograms for two track secondary...")
        hx, hy, hz, he = histogram_results(rsgn.two_track_secondary)
        HSt2s = Dict("hx" => hx, "hy" => hy, "hz" => hz, "he" => he)
        output_file = joinpath(params[:output_dir], "HSt2s.txt")
        save_histos(HSt2s, output_file)
        println("Saved: $output_file")
        
        # Generate and save histograms for three track primary
        println("Generating histograms for three track primary...")
        hx, hy, hz, he = histogram_results(rsgn.three_track_primary)
        HSt3p = Dict("hx" => hx, "hy" => hy, "hz" => hz, "he" => he)
        output_file = joinpath(params[:output_dir], "HSt3p.txt")
        save_histos(HSt3p, output_file)
        println("Saved: $output_file")
        
        # Generate and save histograms for three track secondary
        println("Generating histograms for three track secondary...")
        hx, hy, hz, he = histogram_results(rsgn.three_track_secondary)
        HSt3s = Dict("hx" => hx, "hy" => hy, "hz" => hz, "he" => he)
        output_file = joinpath(params[:output_dir], "HSt3s.txt")
        save_histos(HSt3s, output_file)
        println("Saved: $output_file")
        
        # Save summary statistics
        println("\nSaving summary statistics...")
        summary_file = joinpath(params[:output_dir], "analysis_summary.txt")
        open(summary_file, "w") do io
            println(io, "Petit HDF5 Analysis Summary")
            println(io, "==========================")
            println(io, "Input directory: $(params[:input_dir])")
            println(io, "Input file: $(params[:input_file])")
            println(io, "Events to process: $(params[:events_to_run])")
            println(io, "Voxel size: $(params[:voxel_size_mm]) mm")
            println(io, "Max distance: $(params[:max_distance_mm]) mm")
            println(io, "Energy threshold: $(params[:energy_threshold_kev]) keV")
            println(io, "")
            println(io, "Results:")
            println(io, "--------")
            println(io, "Events processed: $(rsgn.n_events_processed)")
            println(io, "Single track events: $(rsgn.n_single_track)")
            println(io, "Two track events: $(rsgn.n_two_track)")
            println(io, "Three+ track events: $(rsgn.n_three_plus_track)")
            println(io, "Failed events: $(rsgn.n_failed)")
            println(io, "")
            println(io, "Processing time: $(round(time() - start_time, digits=2)) seconds")
        end
        println("Saved: $summary_file")
        
        elapsed_time = time() - start_time
        println("\n=== Analysis Complete ===")
        println("Total processing time: $(round(elapsed_time, digits=2)) seconds")
        println("========================\n")
        
    catch e
        println("\nError during analysis:")
        println(e)
        if isa(e, SystemError) || isa(e, ArgumentError)
            println("\nPlease check that:")
            println("1. The input directory exists: $(params[:input_dir])")
            println("2. The input file exists: $(joinpath(params[:input_dir], params[:input_file]))")
        end
        rethrow(e)
    end
end

# Run the main function if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end