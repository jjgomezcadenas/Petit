#!/usr/bin/env julia

"""
Script to export comparison_df from hd5t_sa.jl analysis results
Creates analysis_results.txt and comparison_df.csv files
"""

using DataFrames
using CSV
using Printf
using JSON
using OrderedCollections
using Dates

"""
Find all .js files in the AnalysisSummary directory
"""
function find_summary_files(summary_dir)
    summary_files = filter(f -> endswith(f, ".js"), readdir(summary_dir))
    return summary_files
end

"""
Safely read a JSON analysis summary file, handling missing keys gracefully.
"""
function safe_read_json_summary(file_path::String)
    try
        # Read the JSON file
        json_string = read(file_path, String)
        data = JSON.parse(json_string, dicttype=OrderedDict)
        
        # Define expected keys with default values
        default_values = OrderedDict(
            "data_type" => "unknown",
            "energy_resolution_keV" => 0.0,
            "cut_track_length_nhits" => 0,
            "cut_radial_mm" => 0.0,
            "cut_z_inner_left_mm" => 0.0,
            "cut_z_inner_right_mm" => 0.0,
            "cut_z_outer_left_mm" => 0.0,
            "cut_z_outer_right_mm" => 0.0,
            "cut_ROI_left_keV" => 0.0,
            "cut_ROI_right_keV" => 0.0,
            "eff_contained" => 0.0,
            "eff_1trk" => 0.0,
            "eff_radial" => 0.0,
            "eff_trkl" => 0.0,
            "eff_zfid" => 0.0,
            "eff_roi" => 0.0,
            "eff_2e" => 0.0,
            "eff_total" => 0.0
        )
        
        # Fill in missing keys with defaults
        result = OrderedDict()
        for (key, default_val) in default_values
            if haskey(data, key)
                result[key] = data[key]
            else
                result[key] = default_val
                @warn "Key '$key' not found in $file_path, using default value: $default_val"
            end
        end
        
        return result
        
    catch e
        @error "Failed to read JSON file: $file_path" exception=(e, catch_backtrace())
        # Return empty result with defaults
        return OrderedDict(
            "data_type" => "error",
            "energy_resolution_keV" => 0.0,
            "cut_track_length_nhits" => 0,
            "cut_radial_mm" => 0.0,
            "cut_z_inner_left_mm" => 0.0,
            "cut_z_inner_right_mm" => 0.0,
            "cut_z_outer_left_mm" => 0.0,
            "cut_z_outer_right_mm" => 0.0,
            "cut_ROI_left_keV" => 0.0,
            "cut_ROI_right_keV" => 0.0,
            "eff_contained" => 0.0,
            "eff_1trk" => 0.0,
            "eff_radial" => 0.0,
            "eff_trkl" => 0.0,
            "eff_zfid" => 0.0,
            "eff_roi" => 0.0,
            "eff_2e" => 0.0,
            "eff_total" => 0.0
        )
    end
end

"""
Load all summary files from the directory.
"""
function load_all_summaries(summary_dir::String, summary_files::Vector{String})
    all_summaries = OrderedDict{String, OrderedDict}()
    
    for filename in summary_files
        local full_path = joinpath(summary_dir, filename)
        local data = safe_read_json_summary(full_path)
        # Use filename without extension as key
        local key = replace(filename, ".js" => "")
        all_summaries[key] = data
    end
    
    return all_summaries
end

"""
Create a comparison DataFrame from all summary data.
"""
function create_comparison_table(all_summaries::OrderedDict)
    # Create DataFrame directly with explicit columns
    data_type = String[]
    contained = Float64[]
    single_track = Float64[]
    radial = Float64[]
    track_length = Float64[]
    z_fiducial = Float64[]
    roi = Float64[]
    two_electron = Float64[]
    total = Float64[]
    
    for (name, data) in all_summaries
        push!(data_type, name)
        push!(contained, round(data["eff_contained"], digits=5))
        push!(single_track, round(data["eff_1trk"], digits=5))
        push!(radial, round(data["eff_radial"], digits=3))
        push!(track_length, round(data["eff_trkl"], digits=3))
        push!(z_fiducial, round(data["eff_zfid"], digits=3))
        push!(roi, round(data["eff_roi"], digits=3))
        push!(two_electron, round(data["eff_2e"], digits=3))
        push!(total, round(data["eff_total"], digits=8))
    end
    
    return DataFrame(
        DataType = data_type,
        Contained = contained,
        SingleTrack = single_track,
        Radial = radial,
        TrackLength = track_length,
        ZFiducial = z_fiducial,
        ROI = roi,
        TwoElectron = two_electron,
        Total = total
    )
end

function main()
    # Define paths
    summary_dir = joinpath(@__DIR__, "AnalysisSummary")
    
    println("Loading summary files from: $summary_dir")
    
    # Find summary files
    summary_files = find_summary_files(summary_dir)
    
    if length(summary_files) == 0
        println("❌ No summary files found!")
        return
    end
    
    println("Found $(length(summary_files)) summary files:")
    for file in summary_files
        println("  - $file")
    end
    
    # Load all summaries
    all_summaries = load_all_summaries(summary_dir, summary_files)
    
    # Create comparison DataFrame
    comparison_df = create_comparison_table(all_summaries)
    
    println("\n📊 Comparison DataFrame created with $(nrow(comparison_df)) rows and $(ncol(comparison_df)) columns")
    
    # Save as CSV
    csv_filename = joinpath(@__DIR__, "comparison_df.csv")
    CSV.write(csv_filename, comparison_df)
    println("✅ DataFrame saved to: $csv_filename")
    
    # Save as formatted text
    txt_filename = joinpath(@__DIR__, "analysis_results.txt")
    open(txt_filename, "w") do io
        println(io, "# HD5t Analysis Results - Efficiency Comparison")
        println(io, "# Generated on: $(now())")
        println(io, "# Number of datasets: $(nrow(comparison_df))")
        println(io, "#" * "="^70)
        println(io, "")
        
        # Write column headers
        col_names = names(comparison_df)
        header_line = @sprintf("%-25s", col_names[1])
        for i in 2:length(col_names)
            header_line *= @sprintf(" %12s", col_names[i])
        end
        println(io, header_line)
        println(io, "-"^(25 + 13*(length(col_names)-1)))
        
        # Write data rows
        for row in eachrow(comparison_df)
            data_line = @sprintf("%-25s", row[1])  # DataType column
            
            # Format numeric columns
            for i in 2:length(col_names)
                data_line *= @sprintf(" %12.6f", row[i])
            end
            
            println(io, data_line)
        end
        
        println(io, "")
        println(io, "#" * "="^70)
        println(io, "# Column Descriptions:")
        println(io, "# DataType    : Analysis dataset name")
        println(io, "# Contained   : Contained events efficiency")
        println(io, "# SingleTrack : Single track selection efficiency") 
        println(io, "# Radial      : Radial cut efficiency")
        println(io, "# TrackLength : Track length cut efficiency")
        println(io, "# ZFiducial   : Z fiducial cut efficiency")
        println(io, "# ROI         : ROI selection efficiency")
        println(io, "# TwoElectron : Two electron identification efficiency")
        println(io, "# Total       : Overall total efficiency")
        println(io, "#" * "="^70)
    end
    
    println("✅ Analysis results saved to: $txt_filename")
    
    # Display summary
    println("\n📋 Summary:")
    println("Dataset with highest total efficiency: $(comparison_df[argmax(comparison_df.Total), :DataType]) ($(maximum(comparison_df.Total)))")
    println("Dataset with lowest total efficiency:  $(comparison_df[argmin(comparison_df.Total), :DataType]) ($(minimum(comparison_df.Total)))")
    
    # Show the DataFrame
    println("\n📊 Comparison DataFrame:")
    show(comparison_df, allrows=true, allcols=true)
    println("")
end

# Run the script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end