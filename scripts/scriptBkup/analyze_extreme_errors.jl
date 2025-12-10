#!/usr/bin/env julia

"""
Analyze why extreme-finding algorithm fails on selected events.

Reads reconstructed/selected track files and compares found extremes
with MC truth to identify failure modes and suggest improvements.
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
using LinearAlgebra

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

"""
    get_mc_extremes(mc_hits_df, voxel_size=1.0, max_distance=3.0)

Get MC truth extremes from original MC hits.
"""
function get_mc_extremes(mc_hits_df::DataFrame; voxel_size::Float64=1.0, max_distance::Float64=3.0)
    mc_vx = Petit.voxelize_event(mc_hits_df, voxel_size)
    mc_dfpars = Petit.DiffusionParams()

    mc_tracks = Petit.make_tracks(mc_vx;
                                  max_distance_mm=max_distance,
                                  energy_threshold_kev=0.0,
                                  diffusion=mc_dfpars)

    if length(mc_tracks) != 1
        return nothing
    end

    mc_track = mc_tracks[1]
    mc_walk = Petit.walk_track_from_extremes(mc_track)

    ext1, ext2 = mc_walk.extremes
    if isnothing(ext1)
        return nothing
    end

    return (pos1=(ext1.x, ext1.y, ext1.z),
            pos2=(ext2.x, ext2.y, ext2.z),
            track=mc_track,
            walk=mc_walk)
end

"""
    distance_3d(p1, p2)

Euclidean distance between two 3D points.
"""
function distance_3d(p1, p2)
    return sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2 + (p1[3]-p2[3])^2)
end

"""
    match_extremes(mc_pos1, mc_pos2, reco_pos1, reco_pos2)

Match MC and reco extremes with minimum total distance.
"""
function match_extremes(mc_pos1, mc_pos2, reco_pos1, reco_pos2)
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

"""
    analyze_graph_structure(track)

Analyze the graph structure to identify potential issues.
"""
function analyze_graph_structure(track::Petit.Tracks)
    g = track.graph
    n_vertices = nv(g)
    n_edges = ne(g)

    degrees = [degree(g, v) for v in vertices(g)]
    deg1_vertices = findall(d -> d == 1, degrees)
    deg3plus_vertices = findall(d -> d >= 3, degrees)  # branch points

    # Find all paths between degree-1 vertices
    branch_info = []
    for v in deg3plus_vertices
        neighbors_list = neighbors(g, v)
        push!(branch_info, (vertex=v, degree=degrees[v], neighbors=neighbors_list))
    end

    return (
        n_vertices=n_vertices,
        n_edges=n_edges,
        avg_degree=mean(degrees),
        max_degree=maximum(degrees),
        n_deg1=length(deg1_vertices),
        n_branches=length(deg3plus_vertices),
        deg1_vertices=deg1_vertices,
        branch_vertices=deg3plus_vertices,
        branch_info=branch_info,
        degrees=degrees
    )
end

"""
    find_closest_vertex(track, pos)

Find the vertex in track closest to a given position.
"""
function find_closest_vertex(track::Petit.Tracks, pos::Tuple{Float64, Float64, Float64})
    voxels = track.voxels
    min_dist = Inf
    closest_v = 1

    for v in 1:nrow(voxels)
        d = distance_3d((voxels.x[v], voxels.y[v], voxels.z[v]), pos)
        if d < min_dist
            min_dist = d
            closest_v = v
        end
    end

    return (vertex=closest_v, distance=min_dist)
end

"""
    analyze_path_to_true_extreme(track, found_extreme, true_extreme)

Analyze why the algorithm didn't find the true extreme.
"""
function analyze_path_to_true_extreme(track::Petit.Tracks,
                                       found_pos::Tuple, true_pos::Tuple)
    g = track.graph
    voxels = track.voxels

    # Find vertices closest to found and true positions
    found_v = find_closest_vertex(track, found_pos)
    true_v = find_closest_vertex(track, true_pos)

    if found_v.vertex == true_v.vertex
        return (same_vertex=true, reason="correct")
    end

    # Check if there's a path between them
    # Use BFS to find shortest path
    path = []
    try
        path = Petit.find_path_bfs(g, found_v.vertex, true_v.vertex)
    catch
        path = []
    end

    # Check degrees along the path
    path_degrees = [degree(g, v) for v in path]
    branches_in_path = findall(d -> d >= 3, path_degrees)

    # Check energies along path
    path_energies = [voxels.energy[v] for v in path]

    # Identify the issue
    reason = "unknown"
    details = Dict{String, Any}()

    if isempty(path)
        reason = "disconnected"
    elseif length(branches_in_path) > 0
        reason = "wrong_branch"
        details["n_branches"] = length(branches_in_path)
        details["branch_positions"] = branches_in_path
        # Find where the algorithm went wrong
        for (i, bp) in enumerate(branches_in_path)
            branch_v = path[bp]
            branch_neighbors = neighbors(g, branch_v)
            details["branch_$(i)_vertex"] = branch_v
            details["branch_$(i)_degree"] = path_degrees[bp]
            details["branch_$(i)_neighbors"] = branch_neighbors
        end
    elseif degree(g, true_v.vertex) != 1
        reason = "true_not_endpoint"
        details["true_degree"] = degree(g, true_v.vertex)
    else
        reason = "topology_confusion"
    end

    return (
        same_vertex=false,
        reason=reason,
        found_vertex=found_v.vertex,
        found_degree=degree(g, found_v.vertex),
        true_vertex=true_v.vertex,
        true_degree=degree(g, true_v.vertex),
        path_length=length(path),
        n_branches_in_path=length(branches_in_path),
        path=path,
        path_degrees=path_degrees,
        path_energies=path_energies,
        details=details
    )
end

"""
    analyze_single_event(reco_result, mc_hits_df; verbose=true)

Analyze a single event comparing reco extremes with MC truth.
"""
function analyze_single_event(reco_result, mc_hits_df::DataFrame; verbose::Bool=true)
    track = reco_result.track
    central_path = reco_result.central_path
    event_id = reco_result.event_id

    # Get reco extremes from the walk
    reco_walk = Petit.walk_track_from_extremes(track)
    ext1, ext2 = reco_walk.extremes

    if isnothing(ext1)
        return (success=false, reason="no_reco_extremes", event_id=event_id)
    end

    reco_pos1 = (ext1.x, ext1.y, ext1.z)
    reco_pos2 = (ext2.x, ext2.y, ext2.z)

    # Get MC truth extremes
    mc_result = get_mc_extremes(mc_hits_df)
    if isnothing(mc_result)
        return (success=false, reason="no_mc_extremes", event_id=event_id)
    end

    mc_pos1, mc_pos2 = mc_result.pos1, mc_result.pos2

    # Match extremes
    match = match_extremes(mc_pos1, mc_pos2, reco_pos1, reco_pos2)

    # Analyze graph structure
    graph_info = analyze_graph_structure(track)

    # Determine which extreme(s) are wrong
    is_ext1_wrong = match.d1 > 5.0  # 5mm threshold
    is_ext2_wrong = match.d2 > 5.0

    analysis1 = nothing
    analysis2 = nothing

    if is_ext1_wrong
        if match.pairing == 1
            analysis1 = analyze_path_to_true_extreme(track, reco_pos1, mc_pos1)
        else
            analysis1 = analyze_path_to_true_extreme(track, reco_pos2, mc_pos1)
        end
    end

    if is_ext2_wrong
        if match.pairing == 1
            analysis2 = analyze_path_to_true_extreme(track, reco_pos2, mc_pos2)
        else
            analysis2 = analyze_path_to_true_extreme(track, reco_pos1, mc_pos2)
        end
    end

    if verbose && (is_ext1_wrong || is_ext2_wrong)
        println("\n" * "=" * 70)
        println("Event $event_id - EXTREME ERROR DETECTED")
        println("=" * 70)
        println("Distance errors: d1=$(round(match.d1, digits=1)) mm, d2=$(round(match.d2, digits=1)) mm")
        println("\nGraph structure:")
        println("  Vertices: $(graph_info.n_vertices), Edges: $(graph_info.n_edges)")
        println("  Deg-1 vertices: $(graph_info.n_deg1), Branch points: $(graph_info.n_branches)")
        println("  Avg degree: $(round(graph_info.avg_degree, digits=2)), Max degree: $(graph_info.max_degree)")

        if is_ext1_wrong && !isnothing(analysis1)
            println("\nExtreme 1 error analysis:")
            println("  Reason: $(analysis1.reason)")
            println("  Found vertex: $(analysis1.found_vertex) (degree=$(analysis1.found_degree))")
            println("  True vertex: $(analysis1.true_vertex) (degree=$(analysis1.true_degree))")
            println("  Path length: $(analysis1.path_length) vertices")
            println("  Branches in path: $(analysis1.n_branches_in_path)")
            if haskey(analysis1.details, "n_branches")
                println("  Branch details: $(analysis1.details)")
            end
        end

        if is_ext2_wrong && !isnothing(analysis2)
            println("\nExtreme 2 error analysis:")
            println("  Reason: $(analysis2.reason)")
            println("  Found vertex: $(analysis2.found_vertex) (degree=$(analysis2.found_degree))")
            println("  True vertex: $(analysis2.true_vertex) (degree=$(analysis2.true_degree))")
            println("  Path length: $(analysis2.path_length) vertices")
            println("  Branches in path: $(analysis2.n_branches_in_path)")
        end
    end

    return (
        success=true,
        event_id=event_id,
        d1=match.d1,
        d2=match.d2,
        is_ext1_wrong=is_ext1_wrong,
        is_ext2_wrong=is_ext2_wrong,
        graph_info=graph_info,
        analysis1=analysis1,
        analysis2=analysis2,
        reco_pos1=reco_pos1,
        reco_pos2=reco_pos2,
        mc_pos1=mc_pos1,
        mc_pos2=mc_pos2,
        reco_confidence=reco_walk.confidence
    )
end

"""
    load_mc_hits(mc_file, event_id)

Load MC hits for a specific event.
"""
function load_mc_hits(mc_file::String, event_id::Int)
    dfs = Petit.get_dataset_dfs(mc_file)
    hits = dfs["hits"]
    evt_hits = Petit.get_event(hits, event_id)

    # Transform to expected format
    df2 = select(evt_hits, Not([:time, :label, :particle_id, :hit_id]))
    f = 1e+5/2.5
    df2.electrons = round.(Int, df2.energy .* f)

    return df2
end

"""
    run_analysis(selected_file, mc_file; n_events=-1, verbose=true)

Run analysis on selected events file.
"""
function run_analysis(selected_file::String, mc_file::String;
                      n_events::Int=-1, verbose::Bool=true,
                      error_threshold::Float64=5.0)

    println("="^70)
    println("EXTREME-FINDING ERROR ANALYSIS")
    println("="^70)
    println("Selected file: $selected_file")
    println("MC file: $mc_file")

    # Load selected events
    println("\nLoading selected events...")
    results, metadata = Petit.read_reco_results_from_hdf5(selected_file)
    println("  Loaded $(length(results)) events")

    if n_events > 0
        results = results[1:min(n_events, length(results))]
    end

    # Analyze each event
    all_analyses = []
    error_reasons = Dict{String, Int}()

    n_ext1_wrong = 0
    n_ext2_wrong = 0
    n_both_wrong = 0
    d1_values = Float64[]
    d2_values = Float64[]

    for (i, r) in enumerate(results)
        if i % 10 == 0
            print("\r  Analyzing event $i/$(length(results))...")
            flush(stdout)
        end

        # Load corresponding MC hits
        mc_hits = load_mc_hits(mc_file, r.event_id)

        # Analyze
        analysis = analyze_single_event(r, mc_hits; verbose=verbose)
        push!(all_analyses, analysis)

        if analysis.success
            push!(d1_values, analysis.d1)
            push!(d2_values, analysis.d2)

            if analysis.is_ext1_wrong
                n_ext1_wrong += 1
                if !isnothing(analysis.analysis1)
                    reason = analysis.analysis1.reason
                    error_reasons[reason] = get(error_reasons, reason, 0) + 1
                end
            end
            if analysis.is_ext2_wrong
                n_ext2_wrong += 1
                if !isnothing(analysis.analysis2)
                    reason = analysis.analysis2.reason
                    error_reasons[reason] = get(error_reasons, reason, 0) + 1
                end
            end
            if analysis.is_ext1_wrong && analysis.is_ext2_wrong
                n_both_wrong += 1
            end
        end
    end
    println("\r  Analyzed $(length(results)) events                ")

    # Summary statistics
    n_total = length(results)
    n_good = count(a -> a.success && !a.is_ext1_wrong && !a.is_ext2_wrong, all_analyses)
    n_any_wrong = n_ext1_wrong + n_ext2_wrong - n_both_wrong

    println("\n" * "-"^70)
    println("SUMMARY")
    println("-"^70)
    println("Total events analyzed: $n_total")
    println("Events with correct extremes: $n_good ($(round(100*n_good/n_total, digits=1))%)")
    println("Events with any error: $n_any_wrong ($(round(100*n_any_wrong/n_total, digits=1))%)")
    println("  - Extreme 1 wrong: $n_ext1_wrong")
    println("  - Extreme 2 wrong: $n_ext2_wrong")
    println("  - Both wrong: $n_both_wrong")

    println("\nError distance statistics:")
    println("  d1: $(round(mean(d1_values), digits=2)) ± $(round(std(d1_values), digits=2)) mm")
    println("  d2: $(round(mean(d2_values), digits=2)) ± $(round(std(d2_values), digits=2)) mm")

    println("\nError reasons breakdown:")
    for (reason, count) in sort(collect(error_reasons), by=x->-x[2])
        println("  $reason: $count ($(round(100*count/(n_ext1_wrong+n_ext2_wrong), digits=1))%)")
    end

    # Analyze wrong_branch cases in detail
    wrong_branch_analyses = filter(a -> a.success &&
        ((!isnothing(a.analysis1) && a.analysis1.reason == "wrong_branch") ||
         (!isnothing(a.analysis2) && a.analysis2.reason == "wrong_branch")),
        all_analyses)

    if length(wrong_branch_analyses) > 0
        println("\n" * "-"^70)
        println("WRONG BRANCH ANALYSIS")
        println("-"^70)
        println("Events with wrong branch: $(length(wrong_branch_analyses))")

        # Collect statistics about branches
        n_branches_when_wrong = Int[]
        branch_degrees = Int[]

        for a in wrong_branch_analyses
            if !isnothing(a.analysis1) && a.analysis1.reason == "wrong_branch"
                push!(n_branches_when_wrong, a.analysis1.n_branches_in_path)
            end
            if !isnothing(a.analysis2) && a.analysis2.reason == "wrong_branch"
                push!(n_branches_when_wrong, a.analysis2.n_branches_in_path)
            end
            push!(branch_degrees, a.graph_info.n_branches)
        end

        println("Branches in path to true extreme: $(round(mean(n_branches_when_wrong), digits=1)) ± $(round(std(n_branches_when_wrong), digits=1))")
        println("Total branch points in graph: $(round(mean(branch_degrees), digits=1)) ± $(round(std(branch_degrees), digits=1))")
    end

    # Suggestions for improvement
    println("\n" * "="^70)
    println("SUGGESTIONS FOR IMPROVEMENT")
    println("="^70)

    if get(error_reasons, "wrong_branch", 0) > 0
        println("""
1. WRONG BRANCH PROBLEM ($(get(error_reasons, "wrong_branch", 0)) cases):
   The algorithm takes a wrong turn at branch points.

   Possible solutions:
   a) Energy-guided branch selection: At each branch point, prefer the
      direction with higher energy gradient (towards Bragg peak).

   b) Path straightness heuristic: Prefer branches that continue in a
      similar direction to the incoming path (tracks are roughly straight).

   c) Distance-to-boundary: Prefer branches that lead towards the track's
      spatial boundary (extremes should be at track edges).

   d) Backtracking: If a branch leads to a dead-end that's not a high-energy
      region, backtrack and try alternative branches.
""")
    end

    if get(error_reasons, "true_not_endpoint", 0) > 0
        println("""
2. TRUE EXTREME NOT A DEGREE-1 VERTEX ($(get(error_reasons, "true_not_endpoint", 0)) cases):
   The true extreme position has degree > 1 in the voxelized graph.

   Possible solutions:
   a) Don't rely solely on topology; use spatial extent as primary criterion.

   b) For dense tracks, find the two points that maximize path length.

   c) Use energy weighting: true extremes have high energy (Bragg peaks).
""")
    end

    if get(error_reasons, "topology_confusion", 0) > 0
        println("""
3. TOPOLOGY CONFUSION ($(get(error_reasons, "topology_confusion", 0)) cases):
   Multiple degree-1 vertices, algorithm picks the wrong ones.

   Possible solutions:
   a) Among all degree-1 pairs, choose the pair with maximum graph distance.

   b) Weight by energy: prefer high-energy degree-1 vertices.

   c) Check that selected extremes span most of the track's spatial extent.
""")
    end

    # Create diagnostic plots
    println("\nGenerating diagnostic plots...")

    # Plot 1: d1 vs d2 scatter
    p1 = scatter(d1_values, d2_values,
                 xlabel="d1 (mm)", ylabel="d2 (mm)",
                 title="Extreme Position Errors",
                 label="", alpha=0.6, markersize=4)
    hline!(p1, [error_threshold], linestyle=:dash, color=:red, label="threshold")
    vline!(p1, [error_threshold], linestyle=:dash, color=:red, label="")

    # Plot 2: Error histogram
    all_errors = vcat(d1_values, d2_values)
    p2 = histogram(all_errors,
                   bins=range(0, min(50, maximum(all_errors)*1.1), length=30),
                   xlabel="Error (mm)", ylabel="Counts",
                   title="Error Distribution",
                   label="", fillalpha=0.7)
    vline!(p2, [error_threshold], linestyle=:dash, color=:red, label="threshold")

    # Plot 3: Error reasons pie chart (as bar chart)
    if !isempty(error_reasons)
        reasons = collect(keys(error_reasons))
        counts = [error_reasons[r] for r in reasons]
        p3 = bar(reasons, counts,
                 xlabel="Reason", ylabel="Count",
                 title="Error Reasons",
                 label="", xrotation=45)
    else
        p3 = plot(title="No errors detected")
    end

    p_combined = plot(p1, p2, p3, layout=(1, 3), size=(1500, 400))
    display(p_combined)

    savefig(p_combined, joinpath(pdir, "scripts", "extreme_error_analysis.png"))
    println("Saved: extreme_error_analysis.png")

    return (analyses=all_analyses, error_reasons=error_reasons,
            d1=d1_values, d2=d2_values,
            n_good=n_good, n_any_wrong=n_any_wrong)
end

"""
    main()

Main function for command-line usage.
"""
function main()
    # Default paths
    selected_file = ""
    mc_file = ""
    n_events = -1
    verbose = false

    # Parse command line arguments
    positional_args = String[]
    for arg in ARGS
        if startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) != 2
                continue
            end
            key, value = parts

            if key == "nevents"
                n_events = parse(Int, value)
            elseif key == "verbose"
                verbose = parse(Bool, value)
            end
        else
            push!(positional_args, arg)
        end
    end

    if length(positional_args) >= 2
        selected_file = positional_args[1]
        mc_file = positional_args[2]
    end

    if isempty(selected_file) || isempty(mc_file)
        println("Usage: julia analyze_extreme_errors.jl <selected_file> <mc_file> [options]")
        println("\nRequired arguments:")
        println("  selected_file    Path to HDF5 file with selected/reconstructed tracks")
        println("  mc_file          Path to original MC HDF5 file")
        println("\nOptional arguments:")
        println("  --nevents=N      Number of events to analyze (default: all)")
        println("  --verbose=true   Print detailed analysis for each error")
        println("\nExample:")
        println("  julia analyze_extreme_errors.jl data/selected.h5 data/mc.h5 --nevents=100 --verbose=true")
        exit(0)
    end

    run_analysis(selected_file, mc_file; n_events=n_events, verbose=verbose)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
