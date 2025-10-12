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
                       energy_threshold_kev::Float64=10.0,
                       emin::Float64=-Inf,
                       emax::Float64=Inf)
    # Convert keV threshold to energy units (assuming data is in MeV)
    energy_threshold = energy_threshold_kev * 1e-3

    # First get the specific event to avoid processing unnecessary data
    event_data = get_event(hitsdf, nevent)

    # Calculate total event energy in keV
    energy = 1e+3 * sum(event_data.energy)
    #println("energy = $energy")

    # Check if event energy is within range
    if energy < emin || energy > emax
        #println("Event $nevent rejected: energy=$energy keV not in range [$emin, $emax] keV")
        return Tracks[]  # Return empty tracks array
    end

    #println("+++event passes energy cut+++")

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


function find_events_with_alphas(partdf::DataFrame)
    # Group by event_id and process each group
    ghdf = groupby(partdf, :event_id)
    energies = Vector{Float64}()
    xs = Vector{Float64}()
    ys = Vector{Float64}()
    zs = Vector{Float64}()
    ids = Vector{Int64}()
    fvs = Vector{String}()  # Assuming final_volume is a String, adjust if needed
    
    for group in ghdf
        # Filter for alphas within this event group
        alphas = filter(row -> row.particle_name == "alpha", group)
        for alpha in eachrow(alphas)
            push!(ids, alpha.event_id)
            push!(energies, alpha.kin_energy)
            push!(xs, alpha.initial_x)
            push!(ys, alpha.initial_y)
            push!(zs, alpha.initial_z)
            push!(fvs, alpha.final_volume)
        end
    end
    
    return DataFrame(event_id=ids, x=xs, y=ys, z=zs, energy=energies, finalv=fvs)
end




function filter_fiducial_events(hitsdf::DataFrame, xyc::Float64, zc::Float64)
    """
    Filter events based on fiducial volume cuts.
    
    An event passes the fiducial cuts if ALL hits in that event satisfy:
    - |x| < xyc (x within ±xyc)
    - |y| < xyc (y within ±xyc)
    - z > zc OR z < -zc (z outside the central band [-zc, zc])
    
    Parameters:
    - hitsdf: DataFrame with columns event_id, x, y, z
    - xyc: Cut value for x and y coordinates
    - zc: Cut value for z coordinate
    
    Returns:
    - DataFrame with events that pass all fiducial cuts
    """
    
    # Group by event_id and apply fiducial cuts
    grouped_df = groupby(hitsdf, :event_id)
    
    result = combine(grouped_df, 
        [:x, :y, :z] => ((x, y, z) -> 
            all(abs.(x) .< xyc) &&                    # All |x| < xyc
            all(abs.(y) .< xyc) &&                    # All |y| < xyc  
            all((z .> zc) .| (z .< -zc))              # All z outside [-zc, zc]
        ) => :passes_fiducial_cuts
    )
    
    # Get event IDs that pass all cuts
    passing_event_ids = result[result.passes_fiducial_cuts, :event_id]
    
    # Filter original DataFrame to keep only passing events
    filtered_df = filter(row -> row.event_id in passing_event_ids, hitsdf)
    
    return filtered_df
end


"""
    filter_radial(hitsdf::DataFrame, xc::Float64, yc::Float64)

Filter events based on radial containment within a circle.

An event passes the radial cut if ALL hits in that event satisfy:
- sqrt(x^2 + y^2) < r, where r = sqrt(xc^2 + yc^2)

This ensures all hits of an event are contained within a circle of radius r
centered at the origin.

# Arguments
- `hitsdf::DataFrame`: DataFrame with columns event_id, x, y, z
- `xc::Float64`: X-coordinate used to define radius (r = sqrt(xc^2 + yc^2))
- `yc::Float64`: Y-coordinate used to define radius (r = sqrt(xc^2 + yc^2))

# Returns
- `DataFrame`: Events where all hits are within the radial cut

# Example
```julia
# Filter events to keep only those with all hits within r=1800mm
filtered_df = filter_radial(hitsdf, 1800.0, 0.0)

# Filter with r = sqrt(1000^2 + 1000^2) ≈ 1414mm
filtered_df = filter_radial(hitsdf, 1000.0, 1000.0)
```
"""
function filter_radial(hitsdf::DataFrame, xc::Float64, yc::Float64)
    # Calculate the radius from xc and yc
    r = sqrt(xc^2 + yc^2)
    
    # Group by event_id and apply radial cuts
    grouped_df = groupby(hitsdf, :event_id)
    
    result = combine(grouped_df, 
        [:x, :y] => ((x, y) -> 
            all(sqrt.(x.^2 .+ y.^2) .< r)  # All hits within radius r
        ) => :passes_radial_cut
    )
    
    # Get event IDs that pass the radial cut
    passing_event_ids = result[result.passes_radial_cut, :event_id]
    
    # Filter original DataFrame to keep only passing events
    filtered_df = filter(row -> row.event_id in passing_event_ids, hitsdf)
    
    return filtered_df
end

"""
    filter_z(hitsdf::DataFrame, zil::Float64, zir::Float64, zol::Float64, zor::Float64)

Filter events based on z-coordinate ranges.

An event passes the z-coordinate filter if ALL hits in that event satisfy:
- (z > zol && z < zil) OR (z > zir && z < zor)

This creates two allowed z-ranges: (zol, zil) and (zir, zor). All hits of an event 
must be within one or both of these ranges.

# Arguments
- `hitsdf::DataFrame`: DataFrame with columns event_id, x, y, z
- `zil::Float64`: Inner left z boundary
- `zir::Float64`: Inner right z boundary  
- `zol::Float64`: Outer left z boundary
- `zor::Float64`: Outer right z boundary

# Returns
- `DataFrame`: Events where all hits satisfy the z-coordinate conditions

# Example
```julia
# Filter events to keep only those with hits in ranges (-200, -100) or (100, 200)
filtered_df = filter_z(hitsdf, -100.0, 100.0, -200.0, 200.0)
```
"""
function filter_z(hitsdf::DataFrame, zil::Float64, zir::Float64, zol::Float64, zor::Float64)
    # Group by event_id and apply z cuts
    grouped_df = groupby(hitsdf, :event_id)
    
    result = combine(grouped_df, 
        [:z] => (z -> 
            all((z .> zol .&& z .< zil) .| (z .> zir .&& z .< zor))  # All hits in allowed z ranges
        ) => :passes_z_cut
    )
    
    # Get event IDs that pass the z cut
    passing_event_ids = result[result.passes_z_cut, :event_id]
    
    # Filter original DataFrame to keep only passing events
    filtered_df = filter(row -> row.event_id in passing_event_ids, hitsdf)
    
    return filtered_df
end

"""
    filter_short_tracks(hitsdf::DataFrame, trkl::Int)

Filter events based on track length (number of hits per event).

An event is kept only if its track length (number of hits) is >= trkl.
Events with fewer hits than trkl are filtered out.

# Arguments
- `hitsdf::DataFrame`: DataFrame with columns event_id, x, y, z, energy
- `trkl::Int`: Minimum track length (number of hits) required to keep an event

# Returns
- `DataFrame`: Events where the number of hits >= trkl

# Example
```julia
# Keep only events with at least 10 hits
filtered_df = filter_short_tracks(hitsdf, 10)

# Keep only events with at least 5 hits
filtered_df = filter_short_tracks(hitsdf, 5)
```
"""
function filter_short_tracks(hitsdf::DataFrame, trkl::Int)
    # Group by event_id and count the number of rows (hits) in each group
    grouped_df = groupby(hitsdf, :event_id)
    
    # Calculate track length for each event
    result = combine(grouped_df, nrow => :track_length)
    
    # Filter to keep only events with track_length >= trkl
    passing_event_ids = result[result.track_length .>= trkl, :event_id]
    
    # Filter original DataFrame to keep only passing events
    filtered_df = filter(row -> row.event_id in passing_event_ids, hitsdf)
    
    return filtered_df
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

"""
    nof_events(xfile::String)

Read the number of events from the HDF5 file configuration.

# Arguments
- `xfile::String`: Full path to the HDF5 file

# Returns
- `Int`: Number of events in the file from configuration
"""
function nof_events(xfile::String)
    h5open(xfile, "r") do fid
        config_array = read(fid["MC/configuration"])
        param_value = config_array[2].param_value
        return parse(Int, param_value)
    end
end

"""
    count_events(cmdir::String, input_file::String)

Count the total number of events with hits in the input HDF5 file.

# Arguments
- `cmdir::String`: Directory containing the input file
- `input_file::String`: Name of the HDF5 input file

# Returns
- `Int`: Total number of unique events with hits in the file
"""
function count_events(cmdir::String, input_file::String)
    filepath = joinpath(cmdir, input_file)
    dfs = get_dataset_dfs(filepath)

    # Assuming the hits data contains event_id column
    if haskey(dfs, "hits")
        return length(unique(dfs["hits"].event_id))
    else
        error("No 'hits' dataset found in HDF5 file")
    end
end
