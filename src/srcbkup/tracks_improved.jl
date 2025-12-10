"""
Improved track extremes finding algorithms that handle sharp turns better.
"""

using LinearAlgebra  # For dot product

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