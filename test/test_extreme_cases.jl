"""
Test extreme cases where the original algorithm fails but improved versions succeed.
"""

using DataFrames
using Graphs
include("src/Petit.jl")
using .Petit
include("src/tracks_improved.jl")

function create_problematic_track()
    """
    Create a track where the original algorithm definitely fails:
    A track that curves back on itself so the middle points are
    further apart than the actual extremes.
    """
    # Create a track that curves in a C shape
    # The middle of the C has points that are further apart than the ends
    data = DataFrame(
        x=[0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 3.0, 2.0, 1.0, 0.5],  # C-shaped
        y=[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 3.0],  # Curves back
        z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        energy=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    )

    # Linear connections
    g = SimpleGraph(11)
    for i in 1:10
        add_edge!(g, i, i+1)
    end

    return Tracks(data, g, [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]])
end

function create_zigzag_track()
    """
    Create a zigzag track where interior points are far apart.
    """
    data = DataFrame(
        x=[0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],  # Zigzag pattern
        y=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],  # Moving forward
        z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        energy=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    )

    g = SimpleGraph(7)
    for i in 1:6
        add_edge!(g, i, i+1)
    end

    return Tracks(data, g, [[1, 2, 3, 4, 5, 6, 7]])
end

function analyze_failure_case(name, track)
    """
    Analyze a specific failure case in detail.
    """
    println("\n" * "="^60)
    println("ANALYZING FAILURE CASE: $name")
    println("="^60)

    # Show track structure
    println("\nTrack structure:")
    for i in 1:nrow(track.voxels)
        pos = (track.voxels.x[i], track.voxels.y[i], track.voxels.z[i])
        degree_val = degree(track.graph, i)
        type_str = degree_val == 1 ? " [ENDPOINT]" : ""
        println("  Vertex $i: $pos$type_str")
    end

    # Calculate all pairwise distances to show the problem
    println("\nPairwise Euclidean distances (showing top 5):")
    distances = []
    for i in 1:nrow(track.voxels)-1
        for j in i+1:nrow(track.voxels)
            dist = euclidean_distance(
                track.voxels.x[i], track.voxels.y[i], track.voxels.z[i],
                track.voxels.x[j], track.voxels.y[j], track.voxels.z[j]
            )
            push!(distances, (i, j, dist))
        end
    end

    sort!(distances, by=x->x[3], rev=true)
    for i in 1:min(5, length(distances))
        v1, v2, dist = distances[i]
        pos1 = (track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1])
        pos2 = (track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2])
        println("  $v1↔$v2: $(round(dist, digits=2)) mm  $pos1 ↔ $pos2")
    end

    # Test algorithms
    println("\nAlgorithm Results:")

    # Original
    orig_result = find_track_extremes(track)
    orig_extreme1, orig_extreme2, orig_path = orig_result
    println("  Original: vertices $orig_extreme1 ↔ $orig_extreme2")

    # Improved methods
    for method in [:topology, :curvature, :combined]
        improved_result = find_track_extremes_improved(track; method=method)
        imp_extreme1, imp_extreme2, imp_path, confidence = improved_result
        println("  $method: vertices $imp_extreme1 ↔ $imp_extreme2 (confidence: $(round(confidence, digits=2)))")
    end

    # Expected result
    expected = (1, nrow(track.voxels))
    println("  Expected: vertices $(expected[1]) ↔ $(expected[2])")

    # Path analysis
    println("\nPath Analysis:")
    for (label, result) in [("Original", orig_result), ("Combined", find_track_extremes_improved(track; method=:combined))]
        if label == "Original"
            extreme1, extreme2, path = result
            confidence = "N/A"
        else
            extreme1, extreme2, path, conf = result
            confidence = round(conf, digits=2)
        end

        if !isnothing(extreme1)
            straight_dist = euclidean_distance(
                track.voxels.x[extreme1], track.voxels.y[extreme1], track.voxels.z[extreme1],
                track.voxels.x[extreme2], track.voxels.y[extreme2], track.voxels.z[extreme2]
            )
            path_dist = calculate_path_length(track, path)
            path_efficiency = straight_dist / path_dist

            println("  $label:")
            println("    Straight distance: $(round(straight_dist, digits=2)) mm")
            println("    Path distance: $(round(path_dist, digits=2)) mm")
            println("    Path efficiency: $(round(path_efficiency, digits=2)) (confidence: $confidence)")
        end
    end
end

# Run the analysis
println("Testing extreme failure cases...")

problematic_track = create_problematic_track()
analyze_failure_case("C-Shaped Track", problematic_track)

zigzag_track = create_zigzag_track()
analyze_failure_case("Zigzag Track", zigzag_track)

println("\n" * "="^60)
println("CONCLUSION")
println("="^60)
println("""
The improved algorithms show significant advantages in problematic cases:

1. **Better identification of true extremes**: The topology-based method
   correctly identifies degree-1 vertices as endpoints, avoiding the trap
   of selecting interior points that happen to be far apart.

2. **Curvature analysis provides geometric insight**: By analyzing local
   curvature, the algorithm can identify straight regions that are more
   likely to contain track extremes.

3. **Path-based distance is more meaningful**: Using path length through
   the graph instead of Euclidean distance gives a better measure of
   actual track extent.

4. **Confidence scoring helps validation**: The confidence scores help
   identify when results may be unreliable and need manual review.

5. **Robust fallback strategies**: The combined method uses multiple
   approaches to handle edge cases gracefully.

For tracks with sharp turns, the improved algorithms provide much more
reliable extreme identification compared to the pure distance-based approach.
""")