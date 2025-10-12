"""
Example script showing how to save tracks to disk for sharing.
"""

using Petit
using DataFrames
using Graphs

println("Creating example tracks...")

# Create a linear track
println("\n1. Creating linear track...")
linear_data = DataFrame(
    x=[0.0, 1.0, 2.0, 3.0, 4.0],
    y=[0.0, 0.0, 0.0, 0.0, 0.0],
    z=[0.0, 0.0, 0.0, 0.0, 0.0],
    energy=[0.1, 0.2, 0.3, 0.4, 0.5]
)

linear_graph = SimpleGraph(5)
for i in 1:4
    add_edge!(linear_graph, i, i+1)
end
linear_track = Tracks(linear_data, linear_graph, [[1, 2, 3, 4, 5]])

print_track_summary(linear_track)
save_track_with_analysis(linear_track, "linear_track")

# Create a sharp turn track
println("\n2. Creating sharp turn track...")
sharp_turn_data = DataFrame(
    x=[0.0, 1.0, 2.0, 3.0, 3.0, 3.0],  # L-shaped
    y=[0.0, 0.0, 0.0, 0.0, 1.0, 2.0],  # Sharp 90° turn
    z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    energy=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
)

sharp_turn_graph = SimpleGraph(6)
for i in 1:5
    add_edge!(sharp_turn_graph, i, i+1)
end
sharp_turn_track = Tracks(sharp_turn_data, sharp_turn_graph, [[1, 2, 3, 4, 5, 6]])

print_track_summary(sharp_turn_track)
save_track_with_analysis(sharp_turn_track, "sharp_turn_track")

# Create a circular track
println("\n3. Creating circular track...")
circle_data = DataFrame(
    x=[0.0, 1.0, 1.0, 0.0],
    y=[0.0, 0.0, 1.0, 1.0],
    z=[0.0, 0.0, 0.0, 0.0],
    energy=[0.2, 0.2, 0.2, 0.2]
)

circle_graph = SimpleGraph(4)
add_edge!(circle_graph, 1, 2)
add_edge!(circle_graph, 2, 3)
add_edge!(circle_graph, 3, 4)
add_edge!(circle_graph, 4, 1)  # Close the circle

circle_track = Tracks(circle_data, circle_graph, [[1, 2, 3, 4]])

print_track_summary(circle_track)
save_track_with_analysis(circle_track, "circular_track")

# Also save as binary for faster loading
println("\n4. Saving binary versions...")
save_track_binary(linear_track, "linear_track.jls")
save_track_binary(sharp_turn_track, "sharp_turn_track.jls")
save_track_binary(circle_track, "circular_track.jls")

println("\n" * "="^60)
println("SAVED TRACKS:")
println("="^60)
println("Files created:")
println("  linear_track_voxels.csv + linear_track_graph.json + linear_track_analysis.json")
println("  sharp_turn_track_voxels.csv + sharp_turn_track_graph.json + sharp_turn_track_analysis.json")
println("  circular_track_voxels.csv + circular_track_graph.json + circular_track_analysis.json")
println("  linear_track.jls")
println("  sharp_turn_track.jls")
println("  circular_track.jls")
println("")
println("To load a track:")
println("  track = load_track_csv(\"linear_track\")")
println("  # or")
println("  track = load_track_binary(\"linear_track.jls\")")
println("")
println("To inspect analysis:")
println("  analysis = JSON.parsefile(\"linear_track_analysis.json\")")

# Test loading to make sure it works
println("\n5. Testing load functionality...")
loaded_track = load_track_csv("linear_track")
println("Loaded track has $(nrow(loaded_track.voxels)) voxels and $(nv(loaded_track.graph)) vertices")

# Compare with original
extremes_original = find_track_extremes(linear_track)
extremes_loaded = find_track_extremes(loaded_track)
println("Original extremes: $(extremes_original[1]) ↔ $(extremes_original[2]) (confidence: $(round(extremes_original[4], digits=3)))")
println("Loaded extremes: $(extremes_loaded[1]) ↔ $(extremes_loaded[2]) (confidence: $(round(extremes_loaded[4], digits=3)))")

println("\nTracks saved successfully! You can now share these files.")