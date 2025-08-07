#!/usr/bin/env julia

"""
Script to read H5 data files and count events
Reuses code from hd5t.jl and the Petit module
"""

using Pkg
project_dir = dirname(@__DIR__)
Pkg.activate(project_dir)

using DataFrames
using HDF5
using Printf
using Petit

function count_events_in_file(filename::String)
    """
    Count the number of events in an H5 file.
    
    Args:
        filename: Path to the H5 file
        
    Returns:
        Tuple of (n_events_particles, n_events_hits, filename)
    """
    try
        # Read the datasets using the Petit function
        dfs = Petit.get_dataset_dfs(filename)
        
        # Extract the dataframes
        hits_df = dfs["hits"]
        particles_df = dfs["particles"]
        
        # Count events in each dataset
        n_events_particles = Petit.number_of_events(particles_df)
        n_events_hits = Petit.number_of_events(hits_df)
        
        return (n_events_particles, n_events_hits, basename(filename))
        
    catch e
        println("Error processing file $filename: $e")
        return (0, 0, basename(filename))
    end
end

function main()
    # Data directory
    data_dir = "/Users/jjgomezcadenas/Data/HD5t"
    
    if !isdir(data_dir)
        println("Error: Data directory $data_dir not found")
        return
    end
    
    # Get all .h5 files in the directory
    h5_files = filter(f -> endswith(f, ".h5"), readdir(data_dir, join=true))
    
    if isempty(h5_files)
        println("No .h5 files found in $data_dir")
        return
    end
    
    println("Found $(length(h5_files)) H5 files in $data_dir")
    println("=" ^ 70)
    
    total_events_particles = 0
    total_events_hits = 0
    
    # Process each file
    for h5_file in h5_files
        n_particles, n_hits, filename = count_events_in_file(h5_file)
        
        @printf("%-30s | Particles: %6d | Hits: %6d\n", 
                filename, n_particles, n_hits)
        
        total_events_particles += n_particles
        total_events_hits += n_hits
    end
    
    println("=" ^ 70)
    @printf("%-30s | Particles: %6d | Hits: %6d\n", 
            "TOTAL", total_events_particles, total_events_hits)
    
    if total_events_particles != total_events_hits
        println("⚠️  Warning: Event counts differ between particles and hits datasets")
    end
end

# Run the script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end