using JSON
using OrderedCollections
using DataFrames

"""
# List all .js files in the AnalysisSummary directory
"""
function find_summary_files(summary_dir)
	
	summary_files = filter(f -> endswith(f, ".js"), readdir(summary_dir))
	
	if length(summary_files) == 0
		md"""
		⚠️ No summary files found in $(summary_dir)
		
		Please run the hd5t_pa.jl notebook first to generate analysis summaries.
		"""
	else
		return summary_files
	end
end

"""
    safe_read_json_summary(file_path::String)

Safely read a JSON analysis summary file, handling missing keys gracefully.
Returns an OrderedDict with default values for missing keys.
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
    load_all_summaries(summary_dir::String, summary_files::Vector{String})

Load all summary files from the directory, handling scope issues properly.
Returns an OrderedDict with filename (without extension) as keys.
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
    create_comparison_table(all_summaries::OrderedDict)

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

"""
    create_cuts_table(summary_data::OrderedDict)

Create a table of analysis cuts/parameters.
"""
function create_cuts_table(summary_data::OrderedDict)
    return DataFrame(
        Parameter = [
            "Data Type",
            "Energy Resolution (keV)",
            "Track Length Cut (nhits)",
            "Radial Cut (mm)",
            "Z Inner Left (mm)",
            "Z Inner Right (mm)",
            "Z Outer Left (mm)",
            "Z Outer Right (mm)",
            "ROI Lower Bound (keV)",
            "ROI Upper Bound (keV)"
        ],
        Value = [
            summary_data["data_type"],
            summary_data["energy_resolution_keV"],
            summary_data["cut_track_length_nhits"],
            summary_data["cut_radial_mm"],
            summary_data["cut_z_inner_left_mm"],
            summary_data["cut_z_inner_right_mm"],
            summary_data["cut_z_outer_left_mm"],
            summary_data["cut_z_outer_right_mm"],
            summary_data["cut_ROI_left_keV"],
            summary_data["cut_ROI_right_keV"]
        ]
    )
end

"""
    create_efficiency_table(summary_data::OrderedDict)

Create a table of efficiency components.
"""
function create_efficiency_table(summary_data::OrderedDict)
    return DataFrame(
        "Efficiency Component" => [
            "Contained Events",
            "Single Track",
            "Radial Cut",
            "Track Length Cut",
            "Z Fiducial Cut",
            "ROI Selection",
            "Two Electron ID",
            "**Total**"
        ],
        "Value" => [
            round(summary_data["eff_contained"], digits=5),
            round(summary_data["eff_1trk"], digits=5),
            round(summary_data["eff_radial"], digits=4),
            round(summary_data["eff_trkl"], digits=4),
            round(summary_data["eff_zfid"], digits=4),
            round(summary_data["eff_roi"], digits=4),
            round(summary_data["eff_2e"], digits=4),
            round(summary_data["eff_total"], digits=8)
        ],
        "Percentage" => [
            "$(round(100*summary_data["eff_contained"], digits=2))%",
            "$(round(100*summary_data["eff_1trk"], digits=2))%",
            "$(round(100*summary_data["eff_radial"], digits=2))%",
            "$(round(100*summary_data["eff_trkl"], digits=2))%",
            "$(round(100*summary_data["eff_zfid"], digits=2))%",
            "$(round(100*summary_data["eff_roi"], digits=2))%",
            "$(round(100*summary_data["eff_2e"], digits=2))%",
            "**$(round(100*summary_data["eff_total"], digits=6))%**"
        ]
    )
end

"""
    create_signal_background_table(all_summaries::OrderedDict)

Create signal vs background comparison table.
Returns tuple (signal_data, sb_table) where signal_data is the signal efficiency
and sb_table is a DataFrame with S/B analysis.
"""
function create_signal_background_table(all_summaries::OrderedDict)
    # Find signal and backgrounds
    signal_data = nothing
    backgrounds = OrderedDict()
    
    for (key, data) in all_summaries
        if key == "znubb"
            signal_data = data
        else
            backgrounds[key] = data
        end
    end
    
    if signal_data === nothing || length(backgrounds) == 0
        return (nothing, DataFrame())
    end
    
    # Calculate signal to background ratios
    sb_data = []
    
    for (name, bkg_data) in backgrounds
        local signal_eff = signal_data["eff_total"]
        local bkg_eff = bkg_data["eff_total"]
        
        # Avoid division by zero
        local sb_ratio, rejection
        if bkg_eff > 0
            sb_ratio = signal_eff / bkg_eff
            rejection = 1 / bkg_eff
        else
            sb_ratio = Inf
            rejection = Inf
        end
        
        push!(sb_data, [
            name,
            round(signal_eff, digits=6),
            round(bkg_eff, sigdigits=3),
            round(sb_ratio, sigdigits=3),
            round(rejection, sigdigits=3)
        ])
    end
    
    sb_df = DataFrame(
        sb_data,
        [:Background, :SignalEff, :BackgroundEff, :SBRatio, :Rejection]
    )
    
    return (signal_data, sb_df)
end

"""
    create_export_table(all_summaries::OrderedDict)

Create a combined export table with all summary data.
"""
function create_export_table(all_summaries::OrderedDict)
    export_data = []
    
    for (key, data) in all_summaries
        local filename = key * ".js"
        push!(export_data, OrderedDict(
            "File" => filename,
            "DataType" => data["data_type"],
            "EnergyRes_keV" => data["energy_resolution_keV"],
            "RadialCut_mm" => data["cut_radial_mm"],
            "TrackLengthCut" => data["cut_track_length_nhits"],
            "ROI_low_keV" => data["cut_ROI_left_keV"],
            "ROI_high_keV" => data["cut_ROI_right_keV"],
            "Eff_Contained" => round(data["eff_contained"], digits=4),
            "Eff_1Track" => round(data["eff_1trk"], digits=4),
            "Eff_Radial" => round(data["eff_radial"], digits=4),
            "Eff_TrackLength" => round(data["eff_trkl"], digits=4),
            "Eff_ZFiducial" => round(data["eff_zfid"], digits=4),
            "Eff_ROI" => round(data["eff_roi"], digits=4),
            "Eff_2e" => round(data["eff_2e"], digits=4),
            "Eff_Total" => round(data["eff_total"], digits=6)
        ))
    end
    
    return DataFrame(export_data)
end

"""
    inspect_json_file(file_path::String)

Debug function to inspect the structure of a JSON file.
Returns a string description of the file contents.
"""
function inspect_json_file(file_path::String)
    try
        json_string = read(file_path, String)
        data = JSON.parse(json_string)
        
        inspection = "File: $file_path\n"
        inspection *= "Keys found: $(join(keys(data), ", "))\n"
        inspection *= "Data types: $(Dict(k => typeof(v) for (k,v) in data))\n"
        
        return inspection
    catch e
        return "Error reading file $file_path: $e"
    end
end