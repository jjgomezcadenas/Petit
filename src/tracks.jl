struct Tracks
    voxels::DataFrame
    graph::SimpleGraph{Int}
    components::Vector{Vector{Int}}
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


function make_tracks(vdf, nevent; max_dist=10.0, energy_thr=1.0)
    """
    Create tracks for a specific event using build_tracks function.
    
    Parameters:
    - vdf: DataFrame with voxelized hits
    - nevent: Event ID to process
    - max_dist: Maximum distance for connecting voxels (mm)
    - energy_thr: Energy threshold for voxels (keV)
    
    Returns:
    - Vector of Tracks objects sorted by energy (highest first)
    """
    tracks = build_tracks(vdf, nevent; max_distance=max_dist, 
                         energy_threshold=energy_thr)

    if length(tracks) > 0
        track_energies = [sum(track.voxels.energy) for track in tracks]
        sorted_indices = sortperm(track_energies, rev=true)
        tracks = tracks[sorted_indices]
    end
	
    return tracks
end