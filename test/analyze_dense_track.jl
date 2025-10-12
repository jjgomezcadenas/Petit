"""
Analyze the dense track with sharp bends and develop improved algorithm.
"""

using Petit
using Graphs
using DataFrames
using Statistics

# Load the track
track = load_track_csv("pluto/track_bends_sharp")

println("=== ANALYZING DENSE TRACK WITH SHARP BENDS ===")
println("Track has $(nv(track.graph)) vertices and $(ne(track.graph)) edges")
println("Average degree: $(round(2*ne(track.graph)/nv(track.graph), digits=2))")

# Find spatial extreme candidates
spatial_candidates = [
    argmin(track.voxels.x),  # X min
    argmax(track.voxels.x),  # X max
    argmin(track.voxels.y),  # Y min
    argmax(track.voxels.y),  # Y max
    argmin(track.voxels.z),  # Z min
    argmax(track.voxels.z)   # Z max
]

# Add minimum degree vertices
min_degree = minimum([degree(track.graph, v) for v in vertices(track.graph)])
min_degree_vertices = [v for v in vertices(track.graph) if degree(track.graph, v) == min_degree]

spatial_candidates = unique([spatial_candidates; min_degree_vertices])
println("Spatial extreme candidates: $spatial_candidates")

# Test all pairs to find the one with maximum path length
max_path_length = 0.0
best_spatial_pair = (0, 0)
best_path = Int[]

println("\nTesting all pairs:")
for i in 1:length(spatial_candidates)-1
    for j in i+1:length(spatial_candidates)
        v1, v2 = spatial_candidates[i], spatial_candidates[j]
        path = find_path_bfs(track.graph, v1, v2)
        if !isempty(path)
            path_length = calculate_path_length(track, path)
            println("  Pair $v1 ↔ $v2: path length $(round(path_length, digits=1)) mm, vertices: $(length(path))")
            if path_length > max_path_length
                max_path_length = path_length
                best_spatial_pair = (v1, v2)
                best_path = path
            end
        end
    end
end

println("\n=== SPATIAL ALGORITHM RESULTS ===")
println("Best spatial pair: $(best_spatial_pair[1]) ↔ $(best_spatial_pair[2])")
println("Path length: $(round(max_path_length, digits=2)) mm")
println("Path vertices: $(length(best_path))")

# Get positions
v1, v2 = best_spatial_pair
pos1 = (track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1])
pos2 = (track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2])
straight_dist = euclidean_distance(pos1[1], pos1[2], pos1[3], pos2[1], pos2[2], pos2[3])
efficiency = straight_dist / max_path_length

println("Positions: $pos1 ↔ $pos2")
println("Straight distance: $(round(straight_dist, digits=2)) mm")
println("Path efficiency: $(round(efficiency, digits=3))")

# Compare with current algorithm result
current_result = find_track_extremes(track)
current_path_length = calculate_path_length(track, current_result[3])

println("\n=== CURRENT ALGORITHM RESULTS ===")
println("Current algorithm: $(current_result[1]) ↔ $(current_result[2])")
println("Current path length: $(round(current_path_length, digits=2)) mm")
println("Current path vertices: $(length(current_result[3]))")
println("Current confidence: $(round(current_result[4], digits=3))")

println("\n=== COMPARISON ===")
println("Improvement in path length: $(round(max_path_length - current_path_length, digits=2)) mm")
println("Improvement in path vertices: $(length(best_path) - length(current_result[3])) more vertices")
println("Spatial algorithm is $(round(100 * (max_path_length - current_path_length) / current_path_length, digits=1))% longer")

function find_extremes_spatial_based(track::Tracks)
    """
    Find extremes for dense tracks using spatial analysis.
    """
    # Get spatial extreme candidates
    spatial_candidates = [
        argmin(track.voxels.x),  # X min
        argmax(track.voxels.x),  # X max
        argmin(track.voxels.y),  # Y min
        argmax(track.voxels.y),  # Y max
        argmin(track.voxels.z),  # Z min
        argmax(track.voxels.z)   # Z max
    ]

    # Add minimum degree vertices (if any exist)
    if nv(track.graph) > 0
        min_degree = minimum([degree(track.graph, v) for v in vertices(track.graph)])
        min_degree_vertices = [v for v in vertices(track.graph) if degree(track.graph, v) == min_degree]
        spatial_candidates = unique([spatial_candidates; min_degree_vertices])
    end

    # Find pair with maximum path length
    max_path_length = 0.0
    best_pair = (spatial_candidates[1], spatial_candidates[min(2, end)])
    best_path = Int[]

    for i in 1:length(spatial_candidates)-1
        for j in i+1:length(spatial_candidates)
            v1, v2 = spatial_candidates[i], spatial_candidates[j]
            path = find_path_bfs(track.graph, v1, v2)
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

    # Calculate confidence based on how much better this is than random pairs
    # Higher confidence for longer paths and better efficiency
    pos1 = (track.voxels.x[best_pair[1]], track.voxels.y[best_pair[1]], track.voxels.z[best_pair[1]])
    pos2 = (track.voxels.x[best_pair[2]], track.voxels.y[best_pair[2]], track.voxels.z[best_pair[2]])
    straight_dist = euclidean_distance(pos1[1], pos1[2], pos1[3], pos2[1], pos2[2], pos2[3])
    efficiency = straight_dist / max_path_length

    # Confidence based on path efficiency and length
    confidence = 0.7 + 0.2 * efficiency + min(0.1, max_path_length / 500.0)
    confidence = min(confidence, 0.95)  # Cap at 0.95

    return (best_pair[1], best_pair[2], best_path, confidence)
end

# Test the new spatial-based algorithm
println("\n=== TESTING NEW SPATIAL ALGORITHM ===")
spatial_result = find_extremes_spatial_based(track)
println("Spatial algorithm: $(spatial_result[1]) ↔ $(spatial_result[2])")
println("Confidence: $(round(spatial_result[4], digits=3))")
println("Path length: $(round(calculate_path_length(track, spatial_result[3]), digits=2)) mm")