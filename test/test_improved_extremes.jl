"""
Test script to demonstrate improved track extremes finding for tracks with sharp turns.
"""

using DataFrames
using Graphs
include("src/Petit.jl")
using .Petit
include("src/tracks_improved.jl")

function create_sharp_turn_track()
    """
    Create a test track that makes a sharp 90-degree turn.
    Current algorithm often fails on this type of track.
    """
    # L-shaped track: horizontal segment then vertical segment
    data = DataFrame(
        x=[0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0],  # Horizontal then vertical
        y=[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0],  # Sharp turn at (4,0)
        z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # All in same plane
        energy=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    )

    # Create linear graph: 1-2-3-4-5-6-7-8-9
    g = SimpleGraph(9)
    for i in 1:8
        add_edge!(g, i, i+1)
    end

    return Tracks(data, g, [[1, 2, 3, 4, 5, 6, 7, 8, 9]])
end

function create_s_curve_track()
    """
    Create a track with an S-curve - another challenging case.
    """
    data = DataFrame(
        x=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
        y=[0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0],  # S-shaped curve
        z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        energy=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    )

    g = SimpleGraph(9)
    for i in 1:8
        add_edge!(g, i, i+1)
    end

    return Tracks(data, g, [[1, 2, 3, 4, 5, 6, 7, 8, 9]])
end

function create_u_turn_track()
    """
    Create a track with a U-turn.
    """
    data = DataFrame(
        x=[0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0, 0.0],
        y=[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0],  # U-shaped
        z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        energy=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    )

    g = SimpleGraph(10)
    for i in 1:9
        add_edge!(g, i, i+1)
    end

    return Tracks(data, g, [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]])
end

function test_extremes_comparison()
    """
    Compare original and improved algorithms on challenging tracks.
    """
    println("=" ^ 80)
    println("COMPARING TRACK EXTREMES ALGORITHMS")
    println("=" ^ 80)

    test_tracks = [
        ("Sharp 90° Turn (L-shape)", create_sharp_turn_track()),
        ("S-Curve", create_s_curve_track()),
        ("U-Turn", create_u_turn_track())
    ]

    for (name, track) in test_tracks
        println("\n📊 Testing: $name")
        println("-" ^ 40)

        # Original algorithm
        orig_result = find_track_extremes(track)
        orig_extreme1, orig_extreme2, orig_path = orig_result

        println("Original Algorithm:")
        println("  Extremes: vertex $orig_extreme1 ↔ vertex $orig_extreme2")
        if !isnothing(orig_extreme1)
            pos1 = (track.voxels.x[orig_extreme1], track.voxels.y[orig_extreme1], track.voxels.z[orig_extreme1])
            pos2 = (track.voxels.x[orig_extreme2], track.voxels.y[orig_extreme2], track.voxels.z[orig_extreme2])
            println("  Positions: $pos1 ↔ $pos2")
            println("  Path length: $(length(orig_path)) vertices")
            straight_dist = euclidean_distance(pos1[1], pos1[2], pos1[3], pos2[1], pos2[2], pos2[3])
            path_dist = calculate_path_length(track, orig_path)
            println("  Straight distance: $(round(straight_dist, digits=2)) mm")
            println("  Path distance: $(round(path_dist, digits=2)) mm")
        end

        # Improved algorithms
        methods = [:topology, :curvature, :combined]
        for method in methods
            improved_result = find_track_extremes_improved(track; method=method)
            imp_extreme1, imp_extreme2, imp_path, confidence = improved_result

            println("\nImproved Algorithm ($method):")
            println("  Extremes: vertex $imp_extreme1 ↔ vertex $imp_extreme2")
            println("  Confidence: $(round(confidence, digits=2))")

            if !isnothing(imp_extreme1)
                pos1 = (track.voxels.x[imp_extreme1], track.voxels.y[imp_extreme1], track.voxels.z[imp_extreme1])
                pos2 = (track.voxels.x[imp_extreme2], track.voxels.y[imp_extreme2], track.voxels.z[imp_extreme2])
                println("  Positions: $pos1 ↔ $pos2")
                println("  Path length: $(length(imp_path)) vertices")
                straight_dist = euclidean_distance(pos1[1], pos1[2], pos1[3], pos2[1], pos2[2], pos2[3])
                path_dist = calculate_path_length(track, imp_path)
                println("  Straight distance: $(round(straight_dist, digits=2)) mm")
                println("  Path distance: $(round(path_dist, digits=2)) mm")
            end
        end

        # Analysis
        println("\n🔍 Analysis:")
        expected_extremes = (1, nrow(track.voxels))  # Should be first and last vertices
        println("  Expected extremes: vertex $(expected_extremes[1]) ↔ vertex $(expected_extremes[2])")

        orig_correct = (orig_extreme1, orig_extreme2) == expected_extremes || (orig_extreme2, orig_extreme1) == expected_extremes
        println("  Original algorithm correct: $(orig_correct ? "✅" : "❌")")

        for method in methods
            improved_result = find_track_extremes_improved(track; method=method)
            imp_extreme1, imp_extreme2, _, confidence = improved_result
            imp_correct = (imp_extreme1, imp_extreme2) == expected_extremes || (imp_extreme2, imp_extreme1) == expected_extremes
            println("  $method algorithm correct: $(imp_correct ? "✅" : "❌") (confidence: $(round(confidence, digits=2)))")
        end

        println()
    end
end

function demonstrate_curvature_analysis()
    """
    Show how curvature analysis helps identify track extremes.
    """
    println("\n" * "=" ^ 80)
    println("CURVATURE ANALYSIS DEMONSTRATION")
    println("=" ^ 80)

    track = create_sharp_turn_track()

    println("\nAnalyzing L-shaped track with sharp turn:")
    println("Vertices and their positions:")

    curvatures = calculate_vertex_curvatures(track)

    for i in 1:nrow(track.voxels)
        pos = (track.voxels.x[i], track.voxels.y[i], track.voxels.z[i])
        curv = curvatures[i]
        degree_val = degree(track.graph, i)

        vertex_type = if degree_val == 1
            "ENDPOINT"
        elseif curv > 0.5
            "SHARP TURN"
        elseif curv > 0.2
            "MILD CURVE"
        else
            "STRAIGHT"
        end

        println("  Vertex $i: $pos, curvature=$(round(curv, digits=3)), degree=$degree_val [$vertex_type]")
    end

    # Show sorted by curvature
    sorted_indices = sortperm(curvatures)
    println("\nVertices sorted by curvature (lowest first - likely extremes):")
    for i in 1:min(5, length(sorted_indices))
        idx = sorted_indices[i]
        pos = (track.voxels.x[idx], track.voxels.y[idx], track.voxels.z[idx])
        println("  Vertex $idx: curvature=$(round(curvatures[idx], digits=3)), position=$pos")
    end
end

# Run the tests
println("Testing improved track extremes algorithms...")
test_extremes_comparison()
demonstrate_curvature_analysis()

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)
println("""
The improved algorithms offer several advantages over the original:

1. **Topology-based method**:
   - Prioritizes degree-1 vertices (true endpoints)
   - Uses path length instead of Euclidean distance
   - More reliable for tracks with clear endpoints

2. **Curvature-based method**:
   - Analyzes local geometry to find straight regions
   - Identifies vertices with low curvature as likely extremes
   - Effective for tracks without clear degree-1 endpoints

3. **Combined method**:
   - Uses multiple strategies for robustness
   - Provides confidence scores for result quality
   - Falls back gracefully when primary methods fail

4. **Key improvements for sharp turns**:
   - Path length vs. straight-line distance
   - Local curvature analysis
   - Better handling of junction vertices
   - Confidence scoring for result validation
""")