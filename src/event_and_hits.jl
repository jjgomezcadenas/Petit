"""
    select_events(hitsdf::DataFrame, nevent::Int; voxel_size_mm=2.0, max_distance_mm=5.0, energy_threshold_kev=10.0)

Process a single event to find tracks by voxelizing hits and applying clustering.

# Arguments
- `hitsdf::DataFrame`: DataFrame containing hit data with columns (event_id, x, y, z, energy)
- `nevent::Int`: Event ID to process
- `voxel_size_mm::Float64=2.0`: Voxel size in mm for spatial discretization
- `max_distance_mm::Float64=5.0`: Maximum distance in mm for connecting voxels into tracks
- `energy_threshold_kev::Float64=10.0`: Minimum energy threshold in keV for including voxels

# Returns
- `Vector{Tracks}`: Array of track objects found in the event, sorted by total energy deposition (highest energy first)
"""
function select_events(hitsdf::DataFrame, nevent::Int; 
                       voxel_size_mm::Float64=2.0, 
                       max_distance_mm::Float64=5.0, 
                       energy_threshold_kev::Float64=10.0)
    # Convert keV threshold to energy units (assuming data is in MeV)
    energy_threshold = energy_threshold_kev * 1e-3
    
    # First get the specific event to avoid processing unnecessary data
    event_data = get_event(hitsdf, nevent)
    
    # Create a temporary DataFrame with just this event for voxelization
    temp_df = DataFrame(
        event_id = [nevent for _ in 1:nrow(event_data)],
        x = event_data.x,
        y = event_data.y, 
        z = event_data.z,
        energy = event_data.energy
    )
    
    # Voxelize just this event
    vdf = voxelize_hits(temp_df, voxel_size_mm)
    
    # Build tracks from voxelized data
    tracks = build_tracks(vdf, nevent; 
                         max_distance=max_distance_mm, 
                         energy_threshold=energy_threshold)
    
    # Sort tracks by total energy deposition (sum of all voxel energies) in descending order
    if length(tracks) > 0
        track_energies = [sum(track.voxels.energy) for track in tracks]
        sorted_indices = sortperm(track_energies, rev=true)
        tracks = tracks[sorted_indices]
    end
    
    return tracks
end

number_of_events(partdf) = length(Set(partdf.event_id))


function get_event(hitsdf::DataFrame, event_id::Int)

    ghdf = groupby(hitsdf, :event_id)
    
    event_data = nothing
    for group in ghdf
        if first(group.event_id) == event_id
            event_data = group
            break
        end
    end

    # Check if event exists
    if event_data === nothing
        error("Event ID $event_id not found in grouped dataframe")
    end

	DataFrame(event_data)
end


function energy_primary(partdf::DataFrame)
    # Group by event_id and process each group
    ghdf = groupby(partdf, :event_id)
    energies = Vector{Float64}()
    
    for group in ghdf
        # Filter for primary particles (mother_id == 0) within this event group
        primary_particles = filter(row -> row.mother_id == 0, group)
        egamma = sum(primary_particles.kin_energy)
        push!(energies, egamma)
    end
    
    return energies
end


function energy_primary(partdf::DataFrame, event_id::Int)
    # Get specific event
    event_particles = filter(:event_id => ==(event_id), partdf)
    
    # Check if event exists
    if isempty(event_particles)
        error("Event ID $event_id not found in dataframe")
    end
    
    # Get primary particles (mother_id == 0)
    primary_particles = filter(:mother_id => ==(0), event_particles)
    
    # Sum their kinetic energy
    return sum(primary_particles.kin_energy)
end


function energy_deposited(hitsdf::DataFrame, event_id::Int)
    # Get specific event
    event_hits = filter(:event_id => ==(event_id), hitsdf)
    
    # Check if event exists
    if isempty(event_hits)
        error("Event ID $event_id not found in dataframe")
    end
    
    return sum(event_hits.energy)
end


function energy_deposited(hitsdf::DataFrame)
    # Group by event_id and process each group
    ghdf = groupby(hitsdf, :event_id)
    energies = Vector{Float64}()
    
    for group in ghdf
        egamma = sum(group.energy)
        push!(energies, egamma)
    end
    
    return energies
end

function hits_per_event(hitsdf::DataFrame, event_id::Int)
    ghdf = groupby(hitsdf, :event_id)
    
    # Find the specific event
    for group in ghdf
        if first(group.event_id) == event_id
            return nrow(group)
        end
    end
    
    # Event not found
    error("Event ID $event_id not found in dataframe")
end


function hits_per_all_events(hitsdf::DataFrame)
    ghdf = groupby(hitsdf, :event_id)
    counts = Vector{Int}(undef, length(ghdf))
    
    for (i, subdf) in enumerate(ghdf)
        counts[i] = nrow(subdf)
    end
    
    return counts
end
