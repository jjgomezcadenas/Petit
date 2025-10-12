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

function find_track_extremes(track::Tracks)
    """
    Find the extreme points (endpoints) of a track by analyzing its graph structure.

    This function now uses an improved algorithm that handles sharp turns better
    by combining topology analysis and curvature calculations.

    Parameters:
    - track: A Tracks object

    Returns:
    - Tuple of (extreme1_idx, extreme2_idx, path, confidence) where:
      - extreme1_idx: Index of first extreme vertex
      - extreme2_idx: Index of second extreme vertex
      - path: Vector of vertex indices showing the path between extremes
      - confidence: Float64 confidence score (0.0-1.0) indicating result reliability
    """
    return find_track_extremes_improved(track; method=:combined)
end

function find_track_extremes_legacy(track::Tracks)
    """
    Legacy version of find_track_extremes using the original distance-based algorithm.

    This function is kept for backward compatibility and comparison purposes.
    For new code, use find_track_extremes() which uses the improved algorithm.

    Parameters:
    - track: A Tracks object

    Returns:
    - Tuple of (extreme1_idx, extreme2_idx, path) where:
      - extreme1_idx: Index of first extreme vertex
      - extreme2_idx: Index of second extreme vertex
      - path: Vector of vertex indices showing the path between extremes
    """
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[])
    elseif n_vertices == 1
        return (1, 1, [1])
    end

    # Find vertices with degree 1 (endpoints)
    endpoints = Int[]
    for v in vertices(g)
        if degree(g, v) == 1
            push!(endpoints, v)
        end
    end

    # Case 1: Linear track with clear endpoints
    if length(endpoints) >= 2
        # Find the pair of endpoints with maximum distance
        max_dist = 0.0
        best_pair = (endpoints[1], endpoints[2])

        for i in 1:length(endpoints)-1
            for j in i+1:length(endpoints)
                v1, v2 = endpoints[i], endpoints[j]
                dist = euclidean_distance(
                    track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
                    track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
                )
                if dist > max_dist
                    max_dist = dist
                    best_pair = (v1, v2)
                end
            end
        end

        # Find path between extremes using BFS
        path = find_path_bfs(g, best_pair[1], best_pair[2])
        return (best_pair[1], best_pair[2], path)

    # Case 2: Circular track or single endpoint - find furthest points
    else
        # Use all vertices if no clear endpoints
        vertices_to_check = length(endpoints) == 1 ?
                            [endpoints[1]; setdiff(vertices(g), endpoints)] :
                            collect(vertices(g))

        # Find pair of vertices with maximum distance
        max_dist = 0.0
        best_pair = (vertices_to_check[1], vertices_to_check[min(2, end)])

        for i in 1:length(vertices_to_check)-1
            for j in i+1:length(vertices_to_check)
                v1, v2 = vertices_to_check[i], vertices_to_check[j]
                dist = euclidean_distance(
                    track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
                    track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
                )
                if dist > max_dist
                    max_dist = dist
                    best_pair = (v1, v2)
                end
            end
        end

        # Find shortest path between extremes
        path = find_path_bfs(g, best_pair[1], best_pair[2])
        return (best_pair[1], best_pair[2], path)
    end
end

function find_path_bfs(g::SimpleGraph, start_vertex::Int, end_vertex::Int)
    """
    Find shortest path between two vertices using breadth-first search.

    Parameters:
    - g: SimpleGraph
    - start_vertex: Starting vertex index
    - end_vertex: Target vertex index

    Returns:
    - Vector of vertex indices representing the path
    """
    if start_vertex == end_vertex
        return [start_vertex]
    end

    # BFS to find shortest path
    queue = [(start_vertex, [start_vertex])]
    visited = Set{Int}()

    while !isempty(queue)
        current, path = popfirst!(queue)

        if current == end_vertex
            return path
        end

        if current in visited
            continue
        end

        push!(visited, current)

        for neighbor in neighbors(g, current)
            if neighbor ∉ visited
                push!(queue, (neighbor, vcat(path, neighbor)))
            end
        end
    end

    # No path found
    return Int[]
end

function walk_track_from_extremes(track::Tracks)
    """
    Walk through a track from one extreme to the other.

    Parameters:
    - track: A Tracks object

    Returns:
    - NamedTuple with:
      - extremes: Tuple of (start_voxel_data, end_voxel_data)
      - path_indices: Vector of vertex indices along the path
      - path_voxels: DataFrame of voxels along the path in order
      - total_length: Total path length in mm
      - confidence: Float64 confidence score for extreme identification
    """
    extreme1, extreme2, path, confidence = find_track_extremes(track)

    if isnothing(extreme1)
        return (extremes = (nothing, nothing),
                path_indices = Int[],
                path_voxels = DataFrame(),
                total_length = 0.0,
                confidence = 0.0)
    end

    # Get voxel data for extremes
    start_voxel = track.voxels[extreme1, :]
    end_voxel = track.voxels[extreme2, :]

    # Get voxels along the path
    path_voxels = track.voxels[path, :]

    # Calculate total path length
    total_length = 0.0
    for i in 1:length(path)-1
        v1, v2 = path[i], path[i+1]
        total_length += euclidean_distance(
            track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
            track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
        )
    end

    return (extremes = (start_voxel, end_voxel),
            path_indices = path,
            path_voxels = path_voxels,
            total_length = total_length,
            confidence = confidence)
end

function energy_in_spheres_around_extremes(track::Tracks, walk_result, radius::Float64)
    """
    Calculate the total energy contained within spheres of given radius around
    the start and end voxels of a track.

    Parameters:
    - track: A Tracks object
    - walk_result: Result from walk_track_from_extremes function
    - radius: Radius of the spheres in mm

    Returns:
    - NamedTuple with:
      - start_sphere_energy: Total energy within radius of start voxel
      - end_sphere_energy: Total energy within radius of end voxel
      - start_voxel_count: Number of voxels within start sphere
      - end_voxel_count: Number of voxels within end sphere
      - start_center: Coordinates of start sphere center (x, y, z)
      - end_center: Coordinates of end sphere center (x, y, z)
    """

    # Check if we have valid extremes
    if isnothing(walk_result.extremes[1])
        return (start_sphere_energy = 0.0,
                end_sphere_energy = 0.0,
                start_voxel_count = 0,
                end_voxel_count = 0,
                start_center = (NaN, NaN, NaN),
                end_center = (NaN, NaN, NaN))
    end

    start_voxel, end_voxel = walk_result.extremes

    # Get coordinates of extreme voxels
    start_center = (start_voxel.x, start_voxel.y, start_voxel.z)
    end_center = (end_voxel.x, end_voxel.y, end_voxel.z)

    # Initialize counters
    start_sphere_energy = 0.0
    end_sphere_energy = 0.0
    start_voxel_count = 0
    end_voxel_count = 0

    # Check each voxel in the track
    for i in 1:nrow(track.voxels)
        voxel_x = track.voxels.x[i]
        voxel_y = track.voxels.y[i]
        voxel_z = track.voxels.z[i]
        voxel_energy = track.voxels.energy[i]

        # Calculate distance to start voxel
        dist_to_start = euclidean_distance(voxel_x, voxel_y, voxel_z,
                                          start_center[1], start_center[2], start_center[3])

        # Calculate distance to end voxel
        dist_to_end = euclidean_distance(voxel_x, voxel_y, voxel_z,
                                        end_center[1], end_center[2], end_center[3])

        # Check if voxel is within start sphere
        if dist_to_start <= radius
            start_sphere_energy += voxel_energy
            start_voxel_count += 1
        end

        # Check if voxel is within end sphere
        if dist_to_end <= radius
            end_sphere_energy += voxel_energy
            end_voxel_count += 1
        end
    end

    # Determine which sphere has more energy
    # blob1 is the sphere with larger energy
    if start_sphere_energy >= end_sphere_energy
        blob1_energy = start_sphere_energy
        blob1_voxel_count = start_voxel_count
        blob1_center = start_center
        blob2_energy = end_sphere_energy
        blob2_voxel_count = end_voxel_count
        blob2_center = end_center
    else
        blob1_energy = end_sphere_energy
        blob1_voxel_count = end_voxel_count
        blob1_center = end_center
        blob2_energy = start_sphere_energy
        blob2_voxel_count = start_voxel_count
        blob2_center = start_center
    end

    # Return both old format (for compatibility) and new format
    return (start_sphere_energy = start_sphere_energy,
            end_sphere_energy = end_sphere_energy,
            start_voxel_count = start_voxel_count,
            end_voxel_count = end_voxel_count,
            start_center = start_center,
            end_center = end_center,
            # New format: blob1 has higher energy than blob2
            blob1_energy = blob1_energy,
            blob1_voxel_count = blob1_voxel_count,
            blob1_center = blob1_center,
            blob2_energy = blob2_energy,
            blob2_voxel_count = blob2_voxel_count,
            blob2_center = blob2_center)
end

function energy_in_spheres_around_extremes(track::Tracks, radius::Float64)
    """
    Convenience function that calls walk_track_from_extremes automatically and
    then calculates sphere energies.

    Parameters:
    - track: A Tracks object
    - radius: Radius of the spheres in mm

    Returns:
    - Same NamedTuple as the main function
    """
    walk_result = walk_track_from_extremes(track)
    return energy_in_spheres_around_extremes(track, walk_result, radius)
end

# =============================================================================
# IMPROVED TRACK EXTREMES ALGORITHMS
# =============================================================================

using LinearAlgebra  # For dot product in curvature calculations

function find_track_extremes_improved(track::Tracks; method::Symbol=:combined)
    """
    Improved algorithm to find track extremes that handles sharp turns better.

    Parameters:
    - track: A Tracks object
    - method: Algorithm to use (:topology, :curvature, :combined)

    Returns:
    - Tuple of (extreme1_idx, extreme2_idx, path, confidence_score)
    """
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0)
    elseif n_vertices == 1
        return (1, 1, [1], 1.0)
    end

    if method == :topology
        return find_extremes_topology_based(track)
    elseif method == :curvature
        return find_extremes_curvature_based(track)
    else  # :combined
        return find_extremes_combined(track)
    end
end

function find_extremes_topology_based(track::Tracks)
    """
    Find extremes based on graph topology - prioritizes degree-1 vertices
    and uses path-length-weighted distance for better handling of curved tracks.
    """
    g = track.graph
    n_vertices = nv(g)

    # Find vertices with degree 1 (true endpoints)
    endpoints = Int[]
    for v in vertices(g)
        if degree(g, v) == 1
            push!(endpoints, v)
        end
    end

    # If we have exactly 2 endpoints, use them (high confidence)
    if length(endpoints) == 2
        path = find_path_bfs(g, endpoints[1], endpoints[2])
        return (endpoints[1], endpoints[2], path, 0.9)
    end

    # If we have more than 2 endpoints, find the pair with maximum path length
    if length(endpoints) > 2
        max_path_length = 0.0
        best_pair = (endpoints[1], endpoints[2])
        best_path = Int[]

        for i in 1:length(endpoints)-1
            for j in i+1:length(endpoints)
                path = find_path_bfs(g, endpoints[i], endpoints[j])
                path_length = calculate_path_length(track, path)

                if path_length > max_path_length
                    max_path_length = path_length
                    best_pair = (endpoints[i], endpoints[j])
                    best_path = path
                end
            end
        end

        return (best_pair[1], best_pair[2], best_path, 0.8)
    end

    # Fallback: no clear endpoints, use distance-based method
    return find_extremes_distance_fallback(track, 0.6)
end

function find_extremes_curvature_based(track::Tracks)
    """
    Find extremes by analyzing local curvature - vertices with lowest curvature
    are likely to be at track ends.
    """
    g = track.graph
    n_vertices = nv(g)

    if n_vertices < 3
        return find_extremes_distance_fallback(track, 0.7)
    end

    # Calculate curvature at each vertex
    curvatures = calculate_vertex_curvatures(track)

    # Find vertices with lowest curvature (most "straight")
    # These are likely to be at track ends
    sorted_indices = sortperm(curvatures)

    # Try pairs of low-curvature vertices
    max_path_length = 0.0
    best_pair = (sorted_indices[1], sorted_indices[2])
    best_path = Int[]
    confidence = 0.0

    # Check first several low-curvature vertices
    n_candidates = min(6, n_vertices)

    for i in 1:n_candidates-1
        for j in i+1:n_candidates
            v1, v2 = sorted_indices[i], sorted_indices[j]

            # Skip if vertices are too close (likely same region)
            if are_vertices_too_close(track, v1, v2)
                continue
            end

            path = find_path_bfs(g, v1, v2)
            if !isempty(path)
                path_length = calculate_path_length(track, path)

                # Bonus for low curvature at endpoints
                curvature_bonus = 2.0 / (curvatures[v1] + curvatures[v2] + 1e-6)
                score = path_length + curvature_bonus

                if score > max_path_length
                    max_path_length = score
                    best_pair = (v1, v2)
                    best_path = path
                    confidence = 0.7 + 0.2 * curvature_bonus / 10.0  # Scale confidence
                end
            end
        end
    end

    return (best_pair[1], best_pair[2], best_path, min(confidence, 0.9))
end

function find_extremes_combined(track::Tracks)
    """
    Combined approach: try multiple methods and pick the best result.
    """
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0)
    end

    # Check if this is a dense track (high connectivity)
    avg_degree = 2 * ne(g) / n_vertices
    is_dense_track = avg_degree > 6.0  # Threshold for dense tracks

    if is_dense_track
        # For dense tracks, use spatial-based method
        spatial_result = find_extremes_spatial_based(track)

        # Also try topology and curvature for comparison
        topo_result = find_extremes_topology_based(track)
        curv_result = find_extremes_curvature_based(track)

        # Compare path lengths and choose the best
        spatial_length = calculate_path_length(track, spatial_result[3])
        topo_length = calculate_path_length(track, topo_result[3])
        curv_length = calculate_path_length(track, curv_result[3])

        # Choose the method that gives the longest path
        if spatial_length >= topo_length && spatial_length >= curv_length
            return spatial_result
        elseif topo_length >= curv_length
            return topo_result
        else
            return curv_result
        end
    else
        # For sparse tracks, use original logic
        # Try topology-based first
        topo_result = find_extremes_topology_based(track)

        # Try curvature-based
        curv_result = find_extremes_curvature_based(track)

        # If topology gives high confidence, use it
        if topo_result[4] >= 0.8
            return topo_result
        end

        # Otherwise, compare path lengths and choose longer path
        topo_length = calculate_path_length(track, topo_result[3])
        curv_length = calculate_path_length(track, curv_result[3])

        if curv_length > topo_length * 1.1  # 10% improvement threshold
            return curv_result
        else
            return topo_result
        end
    end
end

function find_extremes_spatial_based(track::Tracks)
    """
    Find extremes for dense tracks using spatial analysis.
    Works well for tracks with high connectivity where topology-based methods fail.
    """
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0)
    elseif n_vertices == 1
        return (1, 1, [1], 1.0)
    end

    # Get spatial extreme candidates
    spatial_candidates = [
        argmin(track.voxels.x),  # X min
        argmax(track.voxels.x),  # X max
        argmin(track.voxels.y),  # Y min
        argmax(track.voxels.y),  # Y max
        argmin(track.voxels.z),  # Z min
        argmax(track.voxels.z)   # Z max
    ]

    # Add minimum degree vertices (potential extremes even in dense tracks)
    min_degree = minimum([degree(g, v) for v in vertices(g)])
    min_degree_vertices = [v for v in vertices(g) if degree(g, v) == min_degree]

    # Combine and remove duplicates
    all_candidates = unique([spatial_candidates; min_degree_vertices])

    # Find pair with maximum path length
    max_path_length = 0.0
    best_pair = (all_candidates[1], all_candidates[min(2, end)])
    best_path = Int[]

    for i in 1:length(all_candidates)-1
        for j in i+1:length(all_candidates)
            v1, v2 = all_candidates[i], all_candidates[j]
            path = find_path_bfs(g, v1, v2)
            if !isempty(path)
                path_length = calculate_path_length(track, path)
                if path_length > max_path_length
                    max_path_length = path_length
                    best_pair = (v1, v2)
                    best_path = path
                end
            end
        end
    end

    # Calculate confidence based on path efficiency and characteristics
    if !isempty(best_path)
        pos1 = (track.voxels.x[best_pair[1]], track.voxels.y[best_pair[1]], track.voxels.z[best_pair[1]])
        pos2 = (track.voxels.x[best_pair[2]], track.voxels.y[best_pair[2]], track.voxels.z[best_pair[2]])
        straight_dist = euclidean_distance(pos1[1], pos1[2], pos1[3], pos2[1], pos2[2], pos2[3])

        if max_path_length > 0
            efficiency = straight_dist / max_path_length

            # Base confidence for spatial method
            confidence = 0.75

            # Bonus for good efficiency (straight tracks)
            if efficiency > 0.8
                confidence += 0.1
            end

            # Bonus for long paths
            confidence += min(0.1, max_path_length / 500.0)

            # Penalty if endpoints have very high degree (not true extremes)
            deg1 = degree(g, best_pair[1])
            deg2 = degree(g, best_pair[2])
            avg_endpoint_degree = (deg1 + deg2) / 2.0
            avg_track_degree = 2 * ne(g) / nv(g)

            if avg_endpoint_degree < avg_track_degree * 0.7  # Endpoints have lower connectivity
                confidence += 0.05
            end

            confidence = min(confidence, 0.95)  # Cap at 0.95
        else
            confidence = 0.5
        end
    else
        confidence = 0.3
    end

    return (best_pair[1], best_pair[2], best_path, confidence)
end

function calculate_vertex_curvatures(track::Tracks)
    """
    Calculate local curvature at each vertex based on neighboring vertices.
    Lower curvature indicates straighter regions (likely track ends).
    """
    g = track.graph
    n_vertices = nv(g)
    curvatures = zeros(n_vertices)

    for v in vertices(g)
        neighbors = collect(Graphs.neighbors(g, v))

        if length(neighbors) < 2
            curvatures[v] = 0.0  # Endpoint - zero curvature
        elseif length(neighbors) == 2
            # Calculate angle between neighbors
            v1, v2 = neighbors[1], neighbors[2]

            # Vectors from current vertex to neighbors
            vec1 = [track.voxels.x[v1] - track.voxels.x[v],
                   track.voxels.y[v1] - track.voxels.y[v],
                   track.voxels.z[v1] - track.voxels.z[v]]

            vec2 = [track.voxels.x[v2] - track.voxels.x[v],
                   track.voxels.y[v2] - track.voxels.y[v],
                   track.voxels.z[v2] - track.voxels.z[v]]

            # Normalize vectors
            norm1 = sqrt(sum(vec1.^2))
            norm2 = sqrt(sum(vec2.^2))

            if norm1 > 1e-6 && norm2 > 1e-6
                vec1 ./= norm1
                vec2 ./= norm2

                # Dot product gives cos(angle)
                cos_angle = dot(vec1, vec2)
                cos_angle = clamp(cos_angle, -1.0, 1.0)

                # Curvature = 1 - |cos(angle)| (0 for straight, 1 for sharp turn)
                curvatures[v] = 1.0 - abs(cos_angle)
            else
                curvatures[v] = 0.0
            end
        else
            # Multiple neighbors - use average curvature
            total_curvature = 0.0
            count = 0

            for i in 1:length(neighbors)-1
                for j in i+1:length(neighbors)
                    v1, v2 = neighbors[i], neighbors[j]

                    vec1 = [track.voxels.x[v1] - track.voxels.x[v],
                           track.voxels.y[v1] - track.voxels.y[v],
                           track.voxels.z[v1] - track.voxels.z[v]]

                    vec2 = [track.voxels.x[v2] - track.voxels.x[v],
                           track.voxels.y[v2] - track.voxels.y[v],
                           track.voxels.z[v2] - track.voxels.z[v]]

                    norm1 = sqrt(sum(vec1.^2))
                    norm2 = sqrt(sum(vec2.^2))

                    if norm1 > 1e-6 && norm2 > 1e-6
                        vec1 ./= norm1
                        vec2 ./= norm2
                        cos_angle = clamp(dot(vec1, vec2), -1.0, 1.0)
                        total_curvature += 1.0 - abs(cos_angle)
                        count += 1
                    end
                end
            end

            curvatures[v] = count > 0 ? total_curvature / count : 1.0  # High curvature for junctions
        end
    end

    return curvatures
end

function calculate_path_length(track::Tracks, path::Vector{Int})
    """
    Calculate the total length along a path through the track.
    """
    if length(path) < 2
        return 0.0
    end

    total_length = 0.0
    for i in 1:length(path)-1
        v1, v2 = path[i], path[i+1]
        dist = euclidean_distance(
            track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
            track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
        )
        total_length += dist
    end

    return total_length
end

function are_vertices_too_close(track::Tracks, v1::Int, v2::Int, threshold::Float64=2.0)
    """
    Check if two vertices are too close to be considered extremes.
    """
    dist = euclidean_distance(
        track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
        track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
    )
    return dist < threshold
end

function find_extremes_distance_fallback(track::Tracks, confidence::Float64)
    """
    Fallback to original distance-based method with specified confidence.
    """
    g = track.graph
    vertices_to_check = collect(vertices(g))

    max_dist = 0.0
    best_pair = (vertices_to_check[1], vertices_to_check[min(2, end)])

    for i in 1:length(vertices_to_check)-1
        for j in i+1:length(vertices_to_check)
            v1, v2 = vertices_to_check[i], vertices_to_check[j]
            dist = euclidean_distance(
                track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
                track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
            )
            if dist > max_dist
                max_dist = dist
                best_pair = (v1, v2)
            end
        end
    end

    path = find_path_bfs(g, best_pair[1], best_pair[2])
    return (best_pair[1], best_pair[2], path, confidence)
end