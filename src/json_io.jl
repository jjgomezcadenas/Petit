using JSON
using OrderedCollections

"""
    json_analysis_summary(path::String, data_type::String,
                                  trklm::Integer, xr::Float64,
                                  zil::Float64, zir::Float64, 
                                  zol::Float64, zor::Float64,
                                  erex::Float64, roi_xi::Float64, roi_xu::Float64,
                                  eff_contained::Float64,
                                  eff_1trk::Float64,
                                  eff_radial::Float64,
                                  eff_trkl::Float64,
                                  eff_zfid::Float64,
                                  eff_roi::Float64,
                                  eff_2e::Float64,  
                                  eff_total::Float64)

Write analysis summary to JSON file with preserved field order.

Uses OrderedDict to maintain the exact order of fields as they appear in the function.

# Arguments
- `path::String`: Output file path
- `data_type::String`: Type of data (e.g., "znubb", "bi214_copper")
- `trklm::Integer`: Track length cut in number of hits
- `xr::Float64`: Radial cut in mm
- `zil::Float64`: Z inner left boundary in mm
- `zir::Float64`: Z inner right boundary in mm
- `zol::Float64`: Z outer left boundary in mm
- `zor::Float64`: Z outer right boundary in mm
- `erex::Float64`: Energy resolution in keV
- `roi_xi::Float64`: ROI lower bound in keV
- `roi_xu::Float64`: ROI upper bound in keV
- `eff_contained::Float64`: Containment efficiency
- `eff_1trk::Float64`: Single track efficiency
- `eff_radial::Float64`: Radial cut efficiency
- `eff_trkl::Float64`: Track length cut efficiency
- `eff_zfid::Float64`: Z fiducial cut efficiency
- `eff_roi::Float64`: ROI efficiency
- `eff_2e::Float64`: Two electron ID efficiency
- `eff_total::Float64`: Total efficiency

# Example
```julia
json_analysis_summary("analysis.json", "znubb", 50, 1300.0,
                     -100.0, 100.0, -2000.0, 2000.0,
                     12.5, 2440.0, 2500.0,
                     0.95, 0.8, 0.9, 0.85, 0.92, 0.88, 0.8, 0.45)
```
"""
function json_analysis_summary(path::String, data_type::String,
                                       trklm::Integer, xr::Float64,
                                       zil::Float64, zir::Float64, 
                                       zol::Float64, zor::Float64,
                                       erex::Float64, roi_xi::Float64, roi_xu::Float64,
                                       eff_contained::Float64,
                                       eff_1trk::Float64,
                                       eff_radial::Float64,
                                       eff_trkl::Float64,
                                       eff_zfid::Float64,
                                       eff_roi::Float64,
                                       eff_2e::Float64,  
                                       eff_total::Float64)
    
    # Use OrderedDict to preserve field order
    data = OrderedDict(
        "data_type" => data_type,
        "energy_resolution_keV" => erex,
        "cut_track_length_nhits" => trklm,
        "cut_radial_mm" => xr,
        "cut_z_inner_left_mm" => zil,
        "cut_z_inner_right_mm" => zir,
        "cut_z_outer_left_mm" => zol,
        "cut_z_outer_right_mm" => zor,
        "cut_ROI_left_keV" => roi_xi,
        "cut_ROI_right_keV" => roi_xu,
        "eff_contained" => eff_contained,
        "eff_1trk" => eff_1trk,
        "eff_radial" => eff_radial,
        "eff_trkl" => eff_trkl,
        "eff_zfid" => eff_zfid,
        "eff_roi" => eff_roi,
        "eff_2e" => eff_2e,
        "eff_total" => eff_total
    )
    
    open(path, "w") do io
        JSON.print(io, data, 4)  # 4 spaces indentation for readability
    end
end

"""
    read_json_analysis_summary(path::String)

Read analysis summary from JSON file, preserving field order.

# Arguments
- `path::String`: Input file path

# Returns
- `OrderedDict`: Ordered dictionary with all analysis parameters

# Example
```julia
summary = read_json_analysis_summary("analysis.json")
println(summary["data_type"])
println(summary["eff_total"])
```
"""
function read_json_analysis_summary(path::String)
    json_string = read(path, String)
    # Parse JSON into OrderedDict to preserve order
    return JSON.parse(json_string, dicttype=OrderedDict)
end