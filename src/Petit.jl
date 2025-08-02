module Petit

using DataFrames
using Graphs
using StatsBase
using Plots

include("histos.jl")
using .histos



struct Tracks
    voxels::DataFrame
    graph::SimpleGraph{Int}
    components::Vector{Vector{Int}}
end

struct AnalysisResults
    single_track_energies::Vector{Float64}
    two_track_primary::Vector{Float64}
    two_track_secondary::Vector{Float64}
    three_track_primary::Vector{Float64}
    three_track_secondary::Vector{Float64}
    n_events_processed::Int
    n_single_track::Int
    n_two_track::Int
    n_three_plus_track::Int
    n_failed::Int
end



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
- `Vector{Tracks}`: Array of track objects found in the event
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
    
    return tracks
end

"""
    analysis_loop(hitsdf::DataFrame; events_to_run=1:100, voxel_size_mm=2.0, max_distance_mm=5.0, energy_threshold_kev=1.0)

Analyze multiple events to extract track energy distributions.

# Arguments
- `hitsdf::DataFrame`: DataFrame containing hit data
- `events_to_run`: Range or collection of event IDs to process (default: 1:100)
- `voxel_size_mm::Float64=2.0`: Voxel size in mm
- `max_distance_mm::Float64=5.0`: Maximum distance for track clustering
- `energy_threshold_kev::Float64=1.0`: Energy threshold in keV

# Returns
- `AnalysisResults`: Struct containing energy arrays and statistics
"""
function analysis_loop(hitsdf::DataFrame; 
                      events_to_run=1:100, 
                      voxel_size_mm::Float64=2.0,
                      max_distance_mm::Float64=5.0, 
                      energy_threshold_kev::Float64=1.0)
    
    # Initialize energy arrays with proper types
    single_track_energies = Float64[]
    two_track_primary = Float64[]
    two_track_secondary = Float64[]
    three_track_primary = Float64[]
    three_track_secondary = Float64[]
    
    # Initialize counters
    n_single_track = 0
    n_two_track = 0
    n_three_plus_track = 0
    n_failed = 0
    n_events_processed = 0
    
    for nevent in events_to_run
        n_events_processed += 1
        
        try
            tracks = select_events(hitsdf, nevent; 
                                 voxel_size_mm=voxel_size_mm, 
                                 max_distance_mm=max_distance_mm, 
                                 energy_threshold_kev=energy_threshold_kev)

            if length(tracks) == 1
                # Single track event
                energy_kev = 1e+3 * sum(tracks[1].voxels.energy)
                push!(single_track_energies, energy_kev)
                n_single_track += 1
                
            elseif length(tracks) == 2 
                # Two track event
                primary_energy = 1e+3 * sum(tracks[1].voxels.energy)
                secondary_energy = 1e+3 * sum(tracks[2].voxels.energy)
                push!(two_track_primary, primary_energy)
                push!(two_track_secondary, secondary_energy)
                n_two_track += 1
                
            elseif length(tracks) >= 3 
                # Three or more track event
                primary_energy = 1e+3 * sum(tracks[1].voxels.energy)
                secondary_energy = 1e+3 * sum(tracks[2].voxels.energy)
                push!(three_track_primary, primary_energy)
                push!(three_track_secondary, secondary_energy)
                n_three_plus_track += 1
            end
            
        catch e
            println("Warning: Error processing event $nevent: $e")
            n_failed += 1
            continue
        end
    end
    
    return AnalysisResults(
        single_track_energies,
        two_track_primary,
        two_track_secondary,
        three_track_primary,
        three_track_secondary,
        n_events_processed,
        n_single_track,
        n_two_track,
        n_three_plus_track,
        n_failed
    )
end


function euclidean_distance(x1, y1, z1, x2, y2, z2)
    return sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
end


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


function voxel_distances(hitsdf::DataFrame, event_id::Int; max_distance=Inf)
    evtdf = get_event(hitsdf, event_id)
    voxel_distances(evtdf; max_distance=max_distance)
end
    
    
function voxel_distances(df::DataFrame; max_distance=Inf)
    all_distances = Float64[]

    n_voxels = nrow(df)
    
    # Return empty array for events with insufficient voxels
    if n_voxels < 2
        return Float64[]
    end
    
    # Compute all pairwise distances within this event
    for i in 1:(n_voxels-1)
        for j in (i+1):n_voxels
            dist = euclidean_distance(
                df.x[i], df.y[i], df.z[i],
                df.x[j], df.y[j], df.z[j]
            )
            
            # Only include distances within the specified maximum
            if dist <= max_distance
                push!(all_distances, dist)
            end
        end
    end

    return all_distances
end


function voxel_closest_distance(hitsdf::DataFrame, event_id::Int; max_distance=Inf)
    evtdf = get_event(hitsdf, event_id)
    voxel_closest_distance(evtdf; max_distance=max_distance)
end


function voxel_closest_distance(df::DataFrame; max_distance=Inf)

    closest_distances = Float64[]
    
    n_voxels = nrow(df)
    
    # Return empty array for events with insufficient voxels
    if n_voxels < 2
        return Float64[]
    end
        
    # For each voxel, find the closest neighbor
    for i in 1:n_voxels
        min_dist = Inf
        
        # Check distance to all other voxels in the same event
        for j in 1:n_voxels
            if i != j
                dist = euclidean_distance(
                    df.x[i], df.y[i], df.z[i],
                    df.x[j], df.y[j], df.z[j]
                )
                
                if dist < min_dist
                    min_dist = dist
                end
            end
        end
        
        # Only include distances within the specified maximum
        if min_dist <= max_distance && min_dist != Inf
            push!(closest_distances, min_dist)
        end
    end
    
    return closest_distances
end


function voxel_energy(hitsdf::DataFrame, event_id::Int)
    evtdf = get_event(hitsdf, event_id)
    voxel_energy(evtdf)
end


function voxel_energy(df::DataFrame)
    # Simply return the energy column as Float64 array
    return Float64.(df.energy)
end


function voxelize_hits(hitsdf::DataFrame, voxel_size::Float64)
    ghdf = groupby(hitsdf, :event_id)
    voxel_rows = Vector{NamedTuple{(:event_id, :x, :y, :z, :energy), Tuple{Int64, Float64, Float64, Float64, Float64}}}()

    for group in ghdf
        eid = first(group.event_id)
        
        x = Float64.(group.x)
        y = Float64.(group.y)
        z = Float64.(group.z)
        e = Float64.(group.energy)

        # Compute voxel indices (integers), then voxel center positions
        ix = floor.(Int, x ./ voxel_size)
        iy = floor.(Int, y ./ voxel_size)
        iz = floor.(Int, z ./ voxel_size)

        # Combine into a tuple key
        voxel_keys = zip(ix, iy, iz)

        # Aggregate energy per voxel
        voxel_energy = Dict{Tuple{Int, Int, Int}, Float64}()

        for (key, en) in zip(voxel_keys, e)
            voxel_energy[key] = get(voxel_energy, key, 0.0) + en
        end

        # Push one row per voxel
        for ((i, j, k), total_e) in voxel_energy
            push!(voxel_rows, (
                event_id = eid,
                x = (i + 0.5) * voxel_size,
                y = (j + 0.5) * voxel_size,
                z = (k + 0.5) * voxel_size,
                energy = total_e
            ))
        end
    end

    return DataFrame(voxel_rows)
end


function build_tracks(hitsdf::DataFrame, event_id::Int; 
                     max_distance::Float64=1.5, 
                     energy_threshold::Float64=0.0)
    
    event_data = get_event(hitsdf, event_id)
    n_voxels = nrow(event_data)
    if n_voxels == 0
        return Tracks[]
    end
    
    # Filter voxels by energy threshold
    valid_indices = findall(e -> e >= 1e-3*energy_threshold, event_data.energy)
    if isempty(valid_indices)
        return Tracks[]
    end
    
    filtered_data = event_data[valid_indices, :]
    n_filtered = length(valid_indices)
    
    # Create graph with filtered voxels
    g = SimpleGraph(n_filtered)
    
    # Add edges based on distance
    for i in 1:(n_filtered-1)
        for j in (i+1):n_filtered
            dist = euclidean_distance(
                filtered_data.x[i], filtered_data.y[i], filtered_data.z[i],
                filtered_data.x[j], filtered_data.y[j], filtered_data.z[j]
            )
            if dist <= max_distance
                add_edge!(g, i, j)
            end
        end
    end
    
    # Find connected components using Graphs.jl optimized algorithm
    components = connected_components(g)
    
    # Create VGraph objects for each componentI
    graphs = Tracks[]
    
    for component in components
        if !isempty(component)
            component_data = filtered_data[component, :]
            
            # Create subgraph for this component
            subgraph_vertices = length(component)
            subgraph = SimpleGraph(subgraph_vertices)
            
            # Remap edges to new indices
            vertex_map = Dict(old_v => new_v for (new_v, old_v) in enumerate(component))
            
            for edge in edges(g)
                src_old, dst_old = src(edge), dst(edge)
                if src_old in component && dst_old in component
                    src_new = vertex_map[src_old]
                    dst_new = vertex_map[dst_old]
                    add_edge!(subgraph, src_new, dst_new)
                end
            end
            
            push!(graphs, Tracks(component_data, subgraph, [component]))
        end
    end
    
    return graphs
end



function plot_hits_trk(trk::Tracks; nbins::Int=100)
	plot_hits(trk.voxels; nbins)
end


function plot_hits_evt(hitsdf::DataFrame, index::Int; nbins=100)
    eventdf = get_event(hitsdf, index)
	plot_hits(eventdf; nbins)
end


function plot_hits(df::DataFrame; nbins::Int=100)
    
    x = Float64.(df.x)
    y = Float64.(df.y)
    z = Float64.(df.z)
    e = Float64.(df.energy)

    # Compute padded limits (1.3x range)
    xmid, xrange = mean((minimum(x), maximum(x))), maximum(x) - minimum(x)
    ymid, yrange = mean((minimum(y), maximum(y))), maximum(y) - minimum(y)
    zmid, zrange = mean((minimum(z), maximum(z))), maximum(z) - minimum(z)

    xlim = (xmid - 0.65 * xrange, xmid + 0.65 * xrange)
    ylim = (ymid - 0.65 * yrange, ymid + 0.65 * yrange)
    zlim = (zmid - 0.65 * zrange, zmid + 0.65 * zrange)

    cmap = cgrad(:viridis, alpha=1.0)

    # === XY ===
    h_xy = fit(Histogram, (x, y), nbins=nbins)
    wxy = h_xy.weights
    wxy_masked = map(v -> v == 0.0 ? NaN : v, wxy)
    xcenters_xy = diff(h_xy.edges[1]) ./ 2 .+ h_xy.edges[1][1:end-1]
    ycenters_xy = diff(h_xy.edges[2]) ./ 2 .+ h_xy.edges[2][1:end-1]
    p1 = heatmap(xcenters_xy, ycenters_xy, wxy_masked';
        xlabel="x", ylabel="y", title="XY Heatmap",
        xlims=xlim, ylims=ylim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === XZ ===
    h_xz = fit(Histogram, (x, z), nbins=nbins)
    wxz_masked = map(v -> v == 0.0 ? NaN : v, h_xz.weights)
    xcenters_xz = diff(h_xz.edges[1]) ./ 2 .+ h_xz.edges[1][1:end-1]
    zcenters_xz = diff(h_xz.edges[2]) ./ 2 .+ h_xz.edges[2][1:end-1]
    p2 = heatmap(xcenters_xz, zcenters_xz, wxz_masked';
        xlabel="x", ylabel="z", title="XZ Heatmap",
        xlims=xlim, ylims=zlim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === YZ ===
    h_yz = fit(Histogram, (y, z), nbins=nbins)
    wyz_masked = map(v -> v == 0.0 ? NaN : v, h_yz.weights)
    ycenters_yz = diff(h_yz.edges[1]) ./ 2 .+ h_yz.edges[1][1:end-1]
    zcenters_yz = diff(h_yz.edges[2]) ./ 2 .+ h_yz.edges[2][1:end-1]
    p3 = heatmap(ycenters_yz, zcenters_yz, wyz_masked';
        xlabel="y", ylabel="z", title="YZ Heatmap",
        xlims=ylim, ylims=zlim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === 3D scatter ===
    p4 = scatter(x, y, z, marker_z=e, ms=2,
        xlabel="x", ylabel="y", zlabel="z", title="3D Scatter",
        xlims=xlim, ylims=ylim, zlims=zlim,
        colorbar_title="Energy", legend=false, cgrad=cmap)

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
end

export Tracks, AnalysisResults
export get_event, voxelize_hits, euclidean_distance, build_tracks, select_events, analysis_loop
export voxel_distances, voxel_closest_distance, voxel_energy
export hits_per_event, hits_per_all_events
export plot_hits_evt, plot_hits_trk, plot_hits

# Re-export histos functions
export hist1d, hist2d, p1df, step_hist, get_histo1d, Histo1d

end
