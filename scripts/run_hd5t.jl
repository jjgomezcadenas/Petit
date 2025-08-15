#!/usr/bin/env julia

# Batch processing script for Petit HDF5 analysis - Version 2
# This version saves the AnalysisResults structure directly to disk
# Usage: julia run_hd5t2.jl <input_dir> [options]

using Pkg
# Activate the project environment
Pkg.activate(dirname(@__DIR__))

using Petit
using DataFrames
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
        println("  --output-dir=<dir>      Output base directory (default: AnalysisResults)")
        exit(1)
    end
    
    # Default values
    params = Dict(
        :input_dir => args[1],
        :input_file => "0nubb.next.h5",
        :events_to_run => 100,
        :voxel_size_mm => 2.0,
        :max_distance_mm => 10.0,
        :energy_threshold_kev => 10.0,
        :output_dir => "AnalysisResults",
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
    
    # Create base output directory if it doesn't exist
    if !isdir(params[:output_dir])
        mkpath(params[:output_dir])
        println("Created base output directory: $(params[:output_dir])")
    end
    
    # Generate specific output directory name based on input filename
    input_basename = splitext(params[:input_file])[1]  # Remove .h5 extension
    specific_output_dir = joinpath(params[:output_dir], input_basename)
    
    if !isdir(specific_output_dir)
        mkpath(specific_output_dir)
        println("Created specific output directory: $(specific_output_dir)")
    end
    
    # Print run parameters
    println("\n=== Petit HDF5 Analysis Batch Job (Version 2) ===")
    println("Input directory: $(params[:input_dir])")
    println("Input file: $(params[:input_file])")
    println("Events to process: $(params[:events_to_run])")
    println("Voxel size: $(params[:voxel_size_mm]) mm")
    println("Max distance: $(params[:max_distance_mm]) mm")
    println("Energy threshold: $(params[:energy_threshold_kev]) keV")
    println("Base output directory: $(params[:output_dir])")
    println("Specific output directory: $(specific_output_dir)")
    
    # Show output files that will be generated
    summary_filename = input_basename * "_summary.txt"
    println("Results will be saved to: $(specific_output_dir)/ (HDF5 + JSON files)")
    println("Summary file: $(summary_filename)")
    println("==================================================\n")
    
    # Run the event loop
    println("Starting analysis...")
    start_time = time()
    
    try
        # Correct way to call event_loop - it expects cmdir as the base directory
        cmdir = joinpath(ENV["DATA"], "HD5t") 
        
        # The event_loop function expects cmdir and then constructs the full path internally
        # So we need to pass the directory name within cmdir, not the full path
        rsgn = event_loop(cmdir; 
                         input_file=joinpath(params[:input_dir], params[:input_file]),
                         events_to_run=params[:events_to_run], 
                         voxel_size_mm=params[:voxel_size_mm],
                         max_distance_mm=params[:max_distance_mm], 
                         energy_threshold_kev=params[:energy_threshold_kev])
        
        elapsed_time = time() - start_time
        
        println("\nAnalysis complete!")
        println("Events processed: $(rsgn.n_events_processed)")
        println("Single track events: $(rsgn.n_single_track)")
        println("Two track events: $(rsgn.n_two_track)")
        println("Three+ track events: $(rsgn.n_three_plus_track)")
        println("Failed events: $(rsgn.n_failed)")
        println("Processing time: $(round(elapsed_time, digits=2)) seconds")
        
        # Save the complete results structure to the specific directory
        # (specific_output_dir was already created at the beginning)
        save_analysis_results(rsgn, specific_output_dir)
        
        # Save a human-readable summary
        println("\nSaving summary statistics...")
        summary_file = joinpath(specific_output_dir, summary_filename)
        open(summary_file, "w") do io
            println(io, "Petit HDF5 Analysis Summary (Version 2)")
            println(io, "======================================")
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
            println(io, "Data Summary:")
            println(io, "-------------")
            println(io, "Single track entries: $(nrow(rsgn.single_track)) rows")
            println(io, "Two track primary entries: $(nrow(rsgn.two_track_primary)) rows")
            println(io, "Two track secondary entries: $(nrow(rsgn.two_track_secondary)) rows")
            println(io, "Three+ track primary entries: $(nrow(rsgn.three_track_primary)) rows")
            println(io, "Three+ track secondary entries: $(nrow(rsgn.three_track_secondary)) rows")
            println(io, "")
            println(io, "Files Generated:")
            println(io, "---------------")
            println(io, "Results directory: $(specific_output_dir)/ (contains HDF5 + JSON files)")
            println(io, "Summary: $(summary_filename)")
            println(io, "")
            println(io, "Processing time: $(round(elapsed_time, digits=2)) seconds")
            println(io, "")
            println(io, "To load the results in Julia:")
            println(io, "using Petit")
            println(io, "results = load_analysis_results(\"$(specific_output_dir)\")")
        end
        println("Saved: $summary_file")
        
        println("\n=== Analysis Complete ===")
        println("Results saved to: $specific_output_dir")
        println("Summary saved to: $summary_file")
        println("Total processing time: $(round(elapsed_time, digits=2)) seconds")
        println("=========================\n")
        
    catch e
        println("\nError during analysis:")
        println(e)
        if isa(e, SystemError) || isa(e, ArgumentError)
            println("\nPlease check that:")
            println("1. The input directory exists in \$(ENV[\"DATA\"])/HD5t/: $(params[:input_dir])")
            println("2. The input file exists: $(params[:input_file])")
            println("3. Environment variable DATA is set correctly")
        end
        rethrow(e)
    end
end

# Run the main function if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end