#!/usr/bin/env julia

"""
Analyze extreme-finding errors in selected/reconstructed events.

This script reads events from a reconstructed/selected HDF5 file and
compares the reco extremes with MC truth (stored as mc_pos1, mc_pos2)
to understand where the algorithm fails.

Usage:
  julia analyze_reco_extremes.jl <reco_file.h5> [--outlier_threshold=10] [--outdir=DIR]
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using Printf
using DataFrames
using Statistics
using Plots
using Graphs

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit


function distance_3d(p1, p2)
    return sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2 + (p1[3]-p2[3])^2)
end

function get_extremes_from_central_path(central_path::DataFrame)
    """Get extremes positions from a central path DataFrame."""
    if nrow(central_path) < 2
        return nothing, nothing
    end
    # First and last rows of central path are the extremes
    ext1 = (central_path.x[1], central_path.y[1], central_path.z[1])
    ext2 = (central_path.x[end], central_path.y[end], central_path.z[end])
    return ext1, ext2
end

function match_extremes(mc_pos1, mc_pos2, reco_pos1, reco_pos2)
    """Match reco extremes to MC extremes, trying both pairings."""
    d11 = distance_3d(mc_pos1, reco_pos1)
    d22 = distance_3d(mc_pos2, reco_pos2)
    pairing1_total = d11 + d22

    d12 = distance_3d(mc_pos1, reco_pos2)
    d21 = distance_3d(mc_pos2, reco_pos1)
    pairing2_total = d12 + d21

    if pairing1_total <= pairing2_total
        return (d1=d11, d2=d22, pairing=1)
    else
        return (d1=d12, d2=d21, pairing=2)
    end
end

function analyze_track_topology(track::Petit.Tracks)
    """Analyze track topology for diagnostic purposes."""
    g = track.graph
    degrees = [Graphs.degree(g, v) for v in Graphs.vertices(g)]

    deg1_vertices = findall(d -> d == 1, degrees)
    deg3plus_vertices = findall(d -> d >= 3, degrees)

    # Find extent
    vx = track.voxels
    extent = sqrt((maximum(vx.x) - minimum(vx.x))^2 +
                  (maximum(vx.y) - minimum(vx.y))^2 +
                  (maximum(vx.z) - minimum(vx.z))^2)

    # Energy distribution
    total_energy = sum(vx.energy)

    return (
        n_voxels = nrow(vx),
        n_deg1 = length(deg1_vertices),
        n_deg3plus = length(deg3plus_vertices),
        avg_degree = mean(degrees),
        max_degree = maximum(degrees),
        extent = extent,
        total_energy = total_energy,
        deg1_vertices = deg1_vertices
    )
end

function get_vertex_position(track::Petit.Tracks, vertex_idx::Int)
    """Get (x,y,z) position of a vertex in the track."""
    vx = track.voxels
    return (vx.x[vertex_idx], vx.y[vertex_idx], vx.z[vertex_idx])
end

function analyze_event(reco_result; outlier_threshold::Float64=10.0)
    """Analyze a single reconstructed event, comparing reco extremes with MC truth."""
    event_id = reco_result.event_id
    track = reco_result.track
    central_path = reco_result.central_path

    # Get reco extremes from central path
    reco_pos1, reco_pos2 = get_extremes_from_central_path(central_path)

    if isnothing(reco_pos1)
        return (event_type="reco_no_central_path", event_id=event_id)
    end

    # Get MC truth extremes from stored values
    mc_pos1 = reco_result.mc_pos1
    mc_pos2 = reco_result.mc_pos2

    if isnothing(mc_pos1) || isnothing(mc_pos2)
        return (event_type="no_mc_positions", event_id=event_id)
    end

    # Match extremes
    match = match_extremes(mc_pos1, mc_pos2, reco_pos1, reco_pos2)

    # Analyze track topology
    topo = analyze_track_topology(track)

    # Check if reco extremes are at degree-1 vertices
    reco_ext1_deg1 = false
    reco_ext2_deg1 = false

    vx = track.voxels
    for vi in topo.deg1_vertices
        pos = get_vertex_position(track, vi)
        if distance_3d(pos, reco_pos1) < 0.1
            reco_ext1_deg1 = true
        end
        if distance_3d(pos, reco_pos2) < 0.1
            reco_ext2_deg1 = true
        end
    end

    # Find closest reco voxel to each MC extreme
    mc1_closest_dist = Inf
    mc1_closest_idx = 0
    mc2_closest_dist = Inf
    mc2_closest_idx = 0

    for i in 1:nrow(vx)
        pos = (vx.x[i], vx.y[i], vx.z[i])
        d1 = distance_3d(pos, mc_pos1)
        d2 = distance_3d(pos, mc_pos2)
        if d1 < mc1_closest_dist
            mc1_closest_dist = d1
            mc1_closest_idx = i
        end
        if d2 < mc2_closest_dist
            mc2_closest_dist = d2
            mc2_closest_idx = i
        end
    end

    mc1_at_deg1 = mc1_closest_idx in topo.deg1_vertices
    mc2_at_deg1 = mc2_closest_idx in topo.deg1_vertices

    mc1_degree = Graphs.degree(track.graph, mc1_closest_idx)
    mc2_degree = Graphs.degree(track.graph, mc2_closest_idx)

    # Classify event type based on whether MC extremes are at degree-1 vertices
    event_type = "both_mc_at_endpoint"

    if match.d1 > outlier_threshold || match.d2 > outlier_threshold
        if !mc1_at_deg1 && !mc2_at_deg1
            event_type = "both_mc_not_at_endpoint"
        elseif !mc1_at_deg1
            event_type = "mc1_not_at_endpoint"
        elseif !mc2_at_deg1
            event_type = "mc2_not_at_endpoint"
        elseif topo.n_deg1 > 2
            event_type = "spurious_endpoints"
        else
            event_type = "wrong_endpoint_selected"
        end
    end

    return (
        event_type = event_type,
        event_id = event_id,
        d1 = match.d1,
        d2 = match.d2,
        pairing = match.pairing,
        mc_pos1 = mc_pos1,
        mc_pos2 = mc_pos2,
        reco_pos1 = reco_pos1,
        reco_pos2 = reco_pos2,
        n_voxels = topo.n_voxels,
        n_deg1 = topo.n_deg1,
        n_deg3plus = topo.n_deg3plus,
        avg_degree = topo.avg_degree,
        max_degree = topo.max_degree,
        extent = topo.extent,
        tortuosity = reco_result.track_length / topo.extent,
        total_energy = topo.total_energy,
        reco_ext1_deg1 = reco_ext1_deg1,
        reco_ext2_deg1 = reco_ext2_deg1,
        mc1_at_deg1 = mc1_at_deg1,
        mc2_at_deg1 = mc2_at_deg1,
        mc1_degree = mc1_degree,
        mc2_degree = mc2_degree,
        mc1_closest_dist = mc1_closest_dist,
        mc2_closest_dist = mc2_closest_dist,
        confidence = reco_result.confidence,
        track_length = reco_result.track_length
    )
end

function main(reco_file::String;
              outlier_threshold::Float64=10.0,
              output_dir::String="")

    println("="^70)
    println("RECO EXTREMES ANALYSIS")
    println("="^70)
    println("Reco file: $reco_file")
    println("Outlier threshold: $outlier_threshold mm")

    # Load reconstructed events
    println("\nLoading reconstructed events...")
    reco_results, reco_metadata = Petit.read_reco_results_from_hdf5(reco_file)
    println("  Loaded $(length(reco_results)) events")

    # Check that MC positions are available
    if isempty(reco_results)
        error("No events in reco file")
    end

    if !hasproperty(first(reco_results), :mc_pos1) || isnothing(first(reco_results).mc_pos1)
        error("Reco file does not contain MC positions (mc_pos1, mc_pos2)")
    end

    # Analyze each event
    println("\nAnalyzing events...")
    results = []

    for (i, reco_r) in enumerate(reco_results)
        if i % 10 == 0
            print("\r  Processing event $i/$(length(reco_results))...")
            flush(stdout)
        end

        try
            result = analyze_event(reco_r; outlier_threshold=outlier_threshold)
            push!(results, result)
        catch e
            push!(results, (event_type="analysis_error: $e", event_id=reco_r.event_id))
        end
    end
    println("\r  Processed $(length(reco_results)) events.          ")

    # Filter out incomplete events and count them
    n_no_central_path = count(r -> r.event_type == "reco_no_central_path", results)
    n_no_mc_positions = count(r -> r.event_type == "no_mc_positions", results)
    n_analysis_error = count(r -> startswith(r.event_type, "analysis_error"), results)

    complete = filter(r -> r.event_type != "reco_no_central_path" &&
                           r.event_type != "no_mc_positions" &&
                           !startswith(r.event_type, "analysis_error"), results)

    d1_values = [r.d1 for r in complete]
    d2_values = [r.d2 for r in complete]

    println("\n" * "="^70)
    println("RESULTS SUMMARY")
    println("="^70)
    println("\nTotal events:        $(length(results))")
    println("Complete:            $(length(complete))")
    println("No central path:     $n_no_central_path")
    println("No MC positions:     $n_no_mc_positions")
    println("Analysis errors:     $n_analysis_error")

    if isempty(complete)
        println("\nNo complete analyses to report.")
        return
    end

    println("\n--- EXTREME FINDING ACCURACY ---")
    println("d1 (extreme 1): $(round(mean(d1_values), digits=2)) ± $(round(std(d1_values), digits=2)) mm")
    println("d2 (extreme 2): $(round(mean(d2_values), digits=2)) ± $(round(std(d2_values), digits=2)) mm")
    println("Max d1: $(round(maximum(d1_values), digits=1)) mm")
    println("Max d2: $(round(maximum(d2_values), digits=1)) mm")

    # Good vs outlier events
    good_events = filter(r -> r.d1 <= outlier_threshold && r.d2 <= outlier_threshold, complete)
    outlier_events = filter(r -> r.d1 > outlier_threshold || r.d2 > outlier_threshold, complete)

    println("\n--- CLASSIFICATION ---")
    println("Good events (d < $outlier_threshold mm): $(length(good_events)) ($(round(100*length(good_events)/length(complete), digits=1))%)")
    println("Outlier events:                         $(length(outlier_events)) ($(round(100*length(outlier_events)/length(complete), digits=1))%)")

    # Analyze topology differences
    println("\n--- TOPOLOGY ANALYSIS ---")

    if !isempty(good_events)
        println("\nGood events (n=$(length(good_events))):")
        println("  n_voxels:     $(round(mean(r.n_voxels for r in good_events), digits=0)) ± $(round(std([r.n_voxels for r in good_events]), digits=0))")
        println("  n_deg1:       $(round(mean(r.n_deg1 for r in good_events), digits=1)) ± $(round(std([r.n_deg1 for r in good_events]), digits=1))")
        println("  n_deg3+:      $(round(mean(r.n_deg3plus for r in good_events), digits=1)) ± $(round(std([r.n_deg3plus for r in good_events]), digits=1))")
        println("  avg_degree:   $(round(mean(r.avg_degree for r in good_events), digits=2)) ± $(round(std([r.avg_degree for r in good_events]), digits=2))")
        println("  extent:       $(round(mean(r.extent for r in good_events), digits=1)) ± $(round(std([r.extent for r in good_events]), digits=1)) mm")
        println("  track_length: $(round(mean(r.track_length for r in good_events), digits=1)) ± $(round(std([r.track_length for r in good_events]), digits=1)) mm")
        println("  tortuosity:   $(round(mean(r.tortuosity for r in good_events), digits=2)) ± $(round(std([r.tortuosity for r in good_events]), digits=2))")
        println("  mc1_at_deg1:  $(count(r.mc1_at_deg1 for r in good_events))/$(length(good_events))")
        println("  mc2_at_deg1:  $(count(r.mc2_at_deg1 for r in good_events))/$(length(good_events))")
    end

    if !isempty(outlier_events)
        println("\nOutlier events (n=$(length(outlier_events))):")
        println("  n_voxels:     $(round(mean(r.n_voxels for r in outlier_events), digits=0)) ± $(round(std([r.n_voxels for r in outlier_events]), digits=0))")
        println("  n_deg1:       $(round(mean(r.n_deg1 for r in outlier_events), digits=1)) ± $(round(std([r.n_deg1 for r in outlier_events]), digits=1))")
        println("  n_deg3+:      $(round(mean(r.n_deg3plus for r in outlier_events), digits=1)) ± $(round(std([r.n_deg3plus for r in outlier_events]), digits=1))")
        println("  avg_degree:   $(round(mean(r.avg_degree for r in outlier_events), digits=2)) ± $(round(std([r.avg_degree for r in outlier_events]), digits=2))")
        println("  extent:       $(round(mean(r.extent for r in outlier_events), digits=1)) ± $(round(std([r.extent for r in outlier_events]), digits=1)) mm")
        println("  track_length: $(round(mean(r.track_length for r in outlier_events), digits=1)) ± $(round(std([r.track_length for r in outlier_events]), digits=1)) mm")
        println("  tortuosity:   $(round(mean(r.tortuosity for r in outlier_events), digits=2)) ± $(round(std([r.tortuosity for r in outlier_events]), digits=2))")
        println("  mc1_at_deg1:  $(count(r.mc1_at_deg1 for r in outlier_events))/$(length(outlier_events))")
        println("  mc2_at_deg1:  $(count(r.mc2_at_deg1 for r in outlier_events))/$(length(outlier_events))")

        # Error type breakdown
        println("\n  Error type breakdown:")
        event_types = [r.event_type for r in outlier_events]
        for et in unique(event_types)
            n = count(==(et), event_types)
            println("    $et: $n ($(round(100*n/length(outlier_events), digits=1))%)")
        end
    end

    # Detailed analysis of worst outliers
    if !isempty(outlier_events)
        println("\n" * "-"^70)
        println("WORST OUTLIERS (top 10):")
        println("-"^70)

        sorted_outliers = sort(outlier_events, by=r->max(r.d1, r.d2), rev=true)

        for (i, r) in enumerate(sorted_outliers[1:min(10, length(sorted_outliers))])
            println("\n[$i] Event $(r.event_id): d1=$(round(r.d1, digits=1))mm, d2=$(round(r.d2, digits=1))mm")
            println("    MC pos1:   $(round.(r.mc_pos1, digits=1))")
            println("    MC pos2:   $(round.(r.mc_pos2, digits=1))")
            println("    Reco pos1: $(round.(r.reco_pos1, digits=1))")
            println("    Reco pos2: $(round.(r.reco_pos2, digits=1))")
            println("    Topology: $(r.n_voxels) voxels, $(r.n_deg1) deg-1, $(r.n_deg3plus) deg-3+, avg_deg=$(round(r.avg_degree, digits=2))")
            println("    MC1 at deg-1: $(r.mc1_at_deg1) (deg=$(r.mc1_degree), dist_to_reco=$(round(r.mc1_closest_dist, digits=1))mm)")
            println("    MC2 at deg-1: $(r.mc2_at_deg1) (deg=$(r.mc2_degree), dist_to_reco=$(round(r.mc2_closest_dist, digits=1))mm)")
            println("    Event type: $(r.event_type)")
            println("    Extent: $(round(r.extent, digits=1))mm, Track length: $(round(r.track_length, digits=1))mm, Tortuosity: $(round(r.tortuosity, digits=2))")
            println("    Confidence: $(round(r.confidence, digits=2))")
        end
    end

    # Generate plots
    println("\n" * "="^70)
    println("GENERATING PLOTS")
    println("="^70)

    # Determine output directory
    if isempty(output_dir)
        output_dir = dirname(reco_file)
    end
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Plot 1: d1 and d2 histograms
    p1 = histogram(d1_values, bins=range(0, min(60, maximum(d1_values)*1.1), length=30),
                   xlabel="d1 (mm)", ylabel="Counts",
                   title="Extreme 1 Error\nmean=$(round(mean(d1_values), digits=1))mm",
                   label="", fillalpha=0.7, color=:blue)
    vline!(p1, [outlier_threshold], label="threshold", color=:red, linestyle=:dash)

    p2 = histogram(d2_values, bins=range(0, min(60, maximum(d2_values)*1.1), length=30),
                   xlabel="d2 (mm)", ylabel="Counts",
                   title="Extreme 2 Error\nmean=$(round(mean(d2_values), digits=1))mm",
                   label="", fillalpha=0.7, color=:green)
    vline!(p2, [outlier_threshold], label="threshold", color=:red, linestyle=:dash)

    # Plot 2: d1 vs d2 scatter
    p3 = scatter(d1_values, d2_values,
                 xlabel="d1 (mm)", ylabel="d2 (mm)",
                 title="d1 vs d2",
                 label="", alpha=0.6, markersize=4)
    plot!(p3, [0, 60], [0, 60], label="d1=d2", linestyle=:dash, color=:gray)
    hline!(p3, [outlier_threshold], label="", color=:red, linestyle=:dash)
    vline!(p3, [outlier_threshold], label="threshold", color=:red, linestyle=:dash)

    # Plot 3: n_deg1 distribution for good vs outliers
    if !isempty(good_events) && !isempty(outlier_events)
        good_deg1 = [r.n_deg1 for r in good_events]
        outlier_deg1 = [r.n_deg1 for r in outlier_events]

        p4 = histogram(good_deg1, bins=2:12, label="Good", alpha=0.6, color=:blue,
                       xlabel="Number of deg-1 vertices", ylabel="Counts",
                       title="Degree-1 Vertex Count")
        histogram!(p4, outlier_deg1, bins=2:12, label="Outlier", alpha=0.6, color=:red)
    else
        p4 = plot(title="Degree-1 distribution\n(insufficient data)")
    end

    # Combined plot
    p = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))

    plot_file = joinpath(output_dir, "reco_extremes_analysis.png")
    savefig(p, plot_file)
    println("Saved plot to: $plot_file")

    # Recommendations
    println("\n" * "="^70)
    println("RECOMMENDATIONS FOR ALGORITHM IMPROVEMENT")
    println("="^70)

    if !isempty(outlier_events)
        # Check if MC extremes are not at deg-1 vertices
        mc_not_at_endpoint = count(r -> !r.mc1_at_deg1 || !r.mc2_at_deg1, outlier_events)
        spurious = count(r -> r.event_type == "spurious_endpoints", outlier_events)
        wrong_selected = count(r -> r.event_type == "wrong_endpoint_selected", outlier_events)

        println("\nError breakdown:")
        println("  MC extreme not at deg-1 vertex: $mc_not_at_endpoint / $(length(outlier_events))")
        println("  Spurious endpoints:             $spurious / $(length(outlier_events))")
        println("  Wrong endpoint selected:        $wrong_selected / $(length(outlier_events))")

        if mc_not_at_endpoint > length(outlier_events) / 2
            println("\n1. Primary issue: MC truth extremes are NOT at degree-1 vertices")
            println("   -> The topology-based approach fundamentally cannot find these")
            println("   -> Consider: energy-based extreme finding (Bragg peak detection)")
            println("   -> Consider: spatial extent analysis (furthest apart voxels)")
        end

        if spurious > length(outlier_events) / 4
            println("\n2. Issue: Spurious endpoints (more than 2 deg-1 vertices)")
            println("   -> Current algorithm may pick wrong endpoints")
            println("   -> Consider: ranking endpoints by energy (Bragg peak has high E)")
            println("   -> Consider: path length criterion (true extremes maximize path)")
        end

        avg_outlier_deg1 = mean(r.n_deg1 for r in outlier_events)
        avg_good_deg1 = isempty(good_events) ? 0.0 : mean(r.n_deg1 for r in good_events)

        if avg_outlier_deg1 > avg_good_deg1 + 1
            println("\n3. Outliers have MORE deg-1 vertices ($avg_outlier_deg1 vs $avg_good_deg1)")
            println("   -> Need better endpoint selection when multiple candidates exist")
        elseif avg_outlier_deg1 < avg_good_deg1 - 0.5
            println("\n3. Outliers have FEWER deg-1 vertices ($avg_outlier_deg1 vs $avg_good_deg1)")
            println("   -> Dense tracks obscure true endpoints")
            println("   -> Consider different voxel/connectivity parameters for dense tracks")
        end
    end

    return (results=results, good=good_events, outliers=outlier_events,
            d1=d1_values, d2=d2_values)
end

# Command line interface
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia analyze_reco_extremes.jl <reco_file.h5> [options]")
        println("")
        println("Arguments:")
        println("  reco_file.h5       Reconstructed events file with MC positions")
        println("")
        println("Options:")
        println("  --outlier_threshold=X  Distance threshold for outliers in mm (default: 10)")
        println("  --outdir=DIR           Output directory for plots (default: same as reco_file)")
        exit(0)
    end

    reco_file = ARGS[1]
    outlier_threshold = 10.0
    output_dir = ""

    for arg in ARGS[2:end]
        if startswith(arg, "--outlier_threshold=")
            outlier_threshold = parse(Float64, split(arg, "=")[2])
        elseif startswith(arg, "--outdir=")
            output_dir = String(split(arg, "=")[2])
        end
    end

    main(reco_file; outlier_threshold=outlier_threshold, output_dir=output_dir)
end
