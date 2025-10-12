using Petit
using DataFrames
using Plots

# Create test data with a simple linear track
test_data = DataFrame(
    event_id = fill(1, 10),
    x = collect(0.0:1.0:9.0),
    y = collect(0.0:0.5:4.5),
    z = fill(0.0, 10),
    energy = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.4, 0.3, 0.2, 0.1]
)

println("Creating tracks...")
tracks = make_tracks(test_data, 1; max_dist=1.5, energy_thr=0.05)

if length(tracks) > 0
    println("Found $(length(tracks)) track(s)")
    track = tracks[1]

    println("Walking through track to find extremes...")
    walk_result = walk_track_from_extremes(track)

    println("Track has $(nrow(track.voxels)) voxels")
    println("Path length: $(round(walk_result.total_length, digits=2)) mm")

    # Plot without extremes
    println("Creating plot without extremes...")
    p1 = plot_track_with_extremes(track)

    # Plot with extremes and path
    println("Creating plot with extremes and path...")
    p2 = plot_track_with_extremes(track, walk_result)

    # Save plots
    savefig(p1, "track_plot_basic.png")
    savefig(p2, "track_plot_with_extremes.png")

    println("Plots saved as track_plot_basic.png and track_plot_with_extremes.png")
else
    println("No tracks found!")
end