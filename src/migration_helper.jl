# Migration helper for converting old TracksSummary-based results to DataFrame-based results

using DataFrames

"""
    extract_unique_energies(df::DataFrame)

Extract unique energies per event from a DataFrame.
Since energy is repeated for each voxel in an event, this returns one energy value per event.
"""
function extract_unique_energies(df::DataFrame)
    if nrow(df) == 0
        return Float64[]
    end
    unique_df = unique(select(df, [:event_id, :energy]))
    return unique_df.energy
end

"""
    get_track_energies(results::AnalysisResults, track_type::Symbol)

Get unique energy values for a specific track type.

# Arguments
- `results`: AnalysisResults structure
- `track_type`: One of :single_track, :two_track_primary, :two_track_secondary, 
                :three_track_primary, :three_track_secondary

# Returns
Vector of unique energy values (one per event)
"""
function get_track_energies(results::AnalysisResults, track_type::Symbol)
    df = getfield(results, track_type)
    return extract_unique_energies(df)
end

"""
    get_track_positions(results::AnalysisResults, track_type::Symbol)

Get position arrays (x, y, z) for a specific track type.

# Returns
Named tuple (xs=..., ys=..., zs=...)
"""
function get_track_positions(results::AnalysisResults, track_type::Symbol)
    df = getfield(results, track_type)
    return (xs=df.x, ys=df.y, zs=df.z)
end

# Export helper functions
export extract_unique_energies, get_track_energies, get_track_positions