# MST-based Extreme Finding Implementation Plan

## Overview

This document provides a complete implementation plan for adding MST (Minimum Spanning Tree) based extreme finding to the Petit package. Another Claude instance should be able to implement this by following the steps below.

## Problem Statement

The current extreme-finding algorithm (`find_extremes_combined_opt`) fails on U-shaped tracks. For xe137 single-electron events (e.g., event 12), both found extremes end up near the Bragg peak region instead of finding the true track start (low energy) and Bragg peak (high energy).

**Root Cause:**
1. The track graph (KDT or kNN) contains cycles
2. BFS/Dijkstra can take shortcuts across U-shaped bends
3. This causes both extremes to cluster in high-energy regions

**Solution:**
- MST removes all cycles, creating a tree
- Tree diameter algorithm (double BFS) is guaranteed to find the two farthest points
- Preliminary tests showed MST finds extremes (703, 339) vs original (126, 588) for event 12

## Implementation Tasks

### Task 1: Add `find_extremes_mst` function

**File:** `/Users/jjgomezcadenas/Projects/Petit/src/tracks.jl`

**Location:** After the `find_extremes_combined_opt` function (around line 2547, after the closing `end` of that function)

**Add this complete function:**

```julia
"""
    find_extremes_mst(track::Tracks)

Find track extremes using MST + tree diameter algorithm.
Guaranteed to find the two geometrically farthest points along the track topology.

The algorithm:
1. Compute MST of the track graph (removes all cycles)
2. Run BFS from arbitrary node to find one extreme
3. Run BFS from that extreme to find the other
4. The path between them is the tree diameter

# Arguments
- `track::Tracks`: Track struct containing voxels and graph

# Returns
- `(extreme1_idx, extreme2_idx, path, confidence, path_length)` where:
  - `extreme1_idx`: Index of first extreme voxel
  - `extreme2_idx`: Index of second extreme voxel
  - `path`: Vector of vertex indices along the path
  - `confidence`: Always 1.0 for MST (deterministic algorithm)
  - `path_length`: Physical length of path in mm
"""
function find_extremes_mst(track::Tracks)
    g = track.graph
    n = nv(g)

    # Handle trivial cases
    n == 0 && return (nothing, nothing, Int[], 0.0, 0.0)
    n == 1 && return (1, 1, [1], 1.0, 0.0)
    ne(g) == 0 && return (1, n, Int[], 0.0, 0.0)

    # Extract coordinates from voxels
    vox = track.voxels
    coords = hcat(vox.x, vox.y, vox.z)'  # 3×N matrix (each column is a voxel position)

    # Build distance matrix for existing edges only
    # This is a sparse representation - only edges in the graph have non-zero weights
    distmx = zeros(n, n)
    for e in edges(g)
        i, j = src(e), dst(e)
        d = sqrt(sum((coords[:, i] .- coords[:, j]).^2))
        distmx[i, j] = d
        distmx[j, i] = d
    end

    # Compute MST using Kruskal's algorithm from Graphs.jl
    # This removes all cycles, creating a tree with n-1 edges
    mst_edges = kruskal_mst(g, distmx)

    # Build MST as a new graph
    mst = SimpleGraph(n)
    for e in mst_edges
        add_edge!(mst, src(e), dst(e))
    end

    # Tree diameter algorithm: double BFS
    # Step 1: BFS from arbitrary node (vertex 1) to find farthest node
    dists1 = gdistances(mst, 1)
    extreme1 = argmax(dists1)

    # Step 2: BFS from extreme1 to find the actual farthest node (other extreme)
    dists2 = gdistances(mst, extreme1)
    extreme2 = argmax(dists2)

    # Find the path between extremes on the MST
    # Reuse the existing find_path_bfs_opt function
    path = find_path_bfs_opt(mst, extreme1, extreme2)

    # Calculate physical path length (sum of Euclidean distances along path)
    path_length = 0.0
    for i in 1:(length(path)-1)
        p1, p2 = path[i], path[i+1]
        path_length += sqrt(sum((coords[:, p1] .- coords[:, p2]).^2))
    end

    return (extreme1, extreme2, path, 1.0, path_length)
end
```

**Dependencies (already available):**
- `Graphs.jl`: `kruskal_mst`, `gdistances`, `SimpleGraph`, `edges`, `src`, `dst`, `nv`, `ne` (already imported)
- `find_path_bfs_opt`: Already exists in tracks.jl (line 1839)

---

### Task 2: Update `find_track_extremes` function signature

**File:** `/Users/jjgomezcadenas/Projects/Petit/src/tracks.jl`

**Location:** Line 433-439

**Current code:**
```julia
function find_track_extremes(track::Tracks; use_optimized::Bool=true)
    if use_optimized
        return find_track_extremes_opt(track)
    else
        return find_track_extremes_improved(track; method=:combined)
    end
end
```

**Replace with:**
```julia
function find_track_extremes(track::Tracks; use_optimized::Bool=true, method::Symbol=:combined)
    if method == :mst
        return find_extremes_mst(track)
    elseif use_optimized
        return find_track_extremes_opt(track)
    else
        return find_track_extremes_improved(track; method=method)
    end
end
```

This allows calling:
- `find_track_extremes(track)` - default behavior (combined_opt)
- `find_track_extremes(track; method=:mst)` - new MST method
- `find_track_extremes(track; method=:combined)` - explicit combined
- `find_track_extremes(track; use_optimized=false)` - legacy behavior

---

### Task 3: Export the new function

**File:** `/Users/jjgomezcadenas/Projects/Petit/src/Petit.jl`

**Location:** Line 45

**Current code:**
```julia
export find_track_extremes, find_track_extremes_legacy, find_track_extremes_improved, find_extremes_topology_based, find_extremes_curvature_based, find_extremes_combined, find_extremes_spatial_based
```

**Replace with:**
```julia
export find_track_extremes, find_track_extremes_legacy, find_track_extremes_improved, find_extremes_topology_based, find_extremes_curvature_based, find_extremes_combined, find_extremes_spatial_based, find_extremes_mst
```

---

### Task 4: Create test script

**File:** `/Users/jjgomezcadenas/Projects/Petit/scripts/test_extreme_methods.jl` (new file)

**Create this complete script:**

```julia
#!/usr/bin/env julia
"""
Test script to compare extreme-finding methods on selected events.

Usage:
    julia test_extreme_methods.jl [--input xe137_selected.h5] [--events 12,22,86]

This script compares the 'combined' and 'mst' methods for finding track extremes
and reports distances to MC ground truth endpoints.
"""

const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)

using DataFrames
using Graphs
using Printf
using Statistics
using LinearAlgebra
using ArgParse

include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

const DEFAULT_CMDIR = joinpath(ENV["DATA"], "HD5t/itaca")
const DEFAULT_EVENTS = [12, 22, 86, 190, 192, 292, 334, 341, 348]

"""
    test_event(hitsdf, event_id; kwargs...)

Test both extreme-finding methods on a single event.
Returns a Dict with results or nothing if event cannot be processed.
"""
function test_event(hitsdf, event_id;
                    σt=1.012, σl=0.0, voxel_size=2.024,
                    mcvox_size=2.024, max_distance=3.035)

    event_df = Petit.get_event(hitsdf, event_id)
    nrow(event_df) == 0 && return nothing

    # Get MC ground truth path
    mc_path = Petit.compute_mc_path(event_df, mcvox_size; is_double_beta=false)
    mc_start = [mc_path.x[1], mc_path.y[1], mc_path.z[1]]
    mc_end = [mc_path.x[end], mc_path.y[end], mc_path.z[end]]

    # Diffuse and voxelize
    event_mc = Petit.transform_hits_df(event_df)
    diffused_df = Petit.diffuse_xyz_image_mc(event_mc; sigma_t_mm=σt, sigma_l_mm=σl)
    voxels = Petit.voxelize_event(diffused_df, voxel_size)

    # Build track
    tracks = Petit.make_tracks(voxels; max_distance_mm=max_distance)
    length(tracks) != 1 && return nothing
    track = tracks[1]

    coords = hcat(track.voxels.x, track.voxels.y, track.voxels.z)
    energies = track.voxels.energy .* 1000  # Convert to keV

    results = Dict{String, Any}()
    results["event_id"] = event_id
    results["n_voxels"] = nrow(track.voxels)
    results["n_edges"] = ne(track.graph)

    # Test each method
    for (name, method) in [("combined", :combined), ("mst", :mst)]
        try
            extremes = Petit.find_track_extremes(track; method=method)
            e1, e2 = extremes[1], extremes[2]

            pos1 = coords[e1, :]
            pos2 = coords[e2, :]

            # Distance to MC endpoints (find best assignment)
            d1_start = norm(pos1 - mc_start)
            d1_end = norm(pos1 - mc_end)
            d2_start = norm(pos2 - mc_start)
            d2_end = norm(pos2 - mc_end)

            # Best assignment: minimize total distance
            assign1 = d1_start + d2_end  # e1→start, e2→end
            assign2 = d1_end + d2_start  # e1→end, e2→start

            if assign1 <= assign2
                d_to_start, d_to_end = d1_start, d2_end
            else
                d_to_start, d_to_end = d2_start, d1_end
            end

            results["$(name)_e1"] = e1
            results["$(name)_e2"] = e2
            results["$(name)_E1_keV"] = energies[e1]
            results["$(name)_E2_keV"] = energies[e2]
            results["$(name)_d_start_mm"] = d_to_start
            results["$(name)_d_end_mm"] = d_to_end
            results["$(name)_total_d_mm"] = d_to_start + d_to_end
            results["$(name)_confidence"] = extremes[4]
            results["$(name)_path_length"] = extremes[5]
        catch e
            results["$(name)_error"] = string(e)
        end
    end

    return results
end

function parse_commandline()
    s = ArgParseSettings(description="Compare extreme-finding methods")

    @add_arg_table! s begin
        "--input", "-i"
            help = "Input H5 file (relative to DATA/HD5t/itaca or absolute)"
            arg_type = String
            default = "xe137_selected.h5"
        "--events", "-e"
            help = "Comma-separated event IDs to test"
            arg_type = String
            default = ""
        "--cmdir"
            help = "Base directory for input files"
            arg_type = String
            default = DEFAULT_CMDIR
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()

    # Resolve input file path
    input_file = args["input"]
    if !isabspath(input_file)
        # Try current directory first, then cmdir
        if isfile(input_file)
            input_path = input_file
        elseif isfile(joinpath(args["cmdir"], input_file))
            input_path = joinpath(args["cmdir"], input_file)
        else
            error("Input file not found: $input_file")
        end
    else
        input_path = input_file
    end

    # Parse events
    events = if isempty(args["events"])
        DEFAULT_EVENTS
    else
        [parse(Int, s) for s in split(args["events"], ",")]
    end

    println("Loading data from: $input_path")
    hitsdf = Petit.load_data(basename(input_path), dirname(input_path))

    println("\n" * "="^80)
    println("EXTREME FINDING METHOD COMPARISON")
    println("="^80)

    all_results = []

    for event_id in events
        println("\n--- Event $event_id ---")
        result = test_event(hitsdf, event_id)

        if isnothing(result)
            println("  Skipped (no single track)")
            continue
        end

        push!(all_results, result)

        println("  Voxels: $(result["n_voxels"]), Edges: $(result["n_edges"])")

        for method in ["combined", "mst"]
            if haskey(result, "$(method)_e1")
                e1 = result["$(method)_e1"]
                e2 = result["$(method)_e2"]
                E1 = result["$(method)_E1_keV"]
                E2 = result["$(method)_E2_keV"]
                d_start = result["$(method)_d_start_mm"]
                d_end = result["$(method)_d_end_mm"]
                total = result["$(method)_total_d_mm"]
                path_len = result["$(method)_path_length"]

                println("  $method: extremes=($e1,$e2), E=($(round(E1,digits=1)),$(round(E2,digits=1)))keV")
                println("          d_start=$(round(d_start,digits=2))mm, d_end=$(round(d_end,digits=2))mm, total=$(round(total,digits=2))mm")
                println("          path_length=$(round(path_len,digits=2))mm")
            else
                println("  $method: ERROR - $(get(result, "$(method)_error", "unknown"))")
            end
        end
    end

    # Summary statistics
    println("\n" * "="^80)
    println("SUMMARY")
    println("="^80)

    for method in ["combined", "mst"]
        totals = [r["$(method)_total_d_mm"] for r in all_results if haskey(r, "$(method)_total_d_mm")]
        if !isempty(totals)
            println("$method:")
            println("  mean_total_d = $(round(mean(totals), digits=2))mm")
            println("  median_total_d = $(round(median(totals), digits=2))mm")
            println("  min_total_d = $(round(minimum(totals), digits=2))mm")
            println("  max_total_d = $(round(maximum(totals), digits=2))mm")
        end
    end

    # Improvement summary
    println("\n--- IMPROVEMENT (MST vs Combined) ---")
    improvements = Float64[]
    for r in all_results
        if haskey(r, "combined_total_d_mm") && haskey(r, "mst_total_d_mm")
            improvement = r["combined_total_d_mm"] - r["mst_total_d_mm"]
            push!(improvements, improvement)
            event = r["event_id"]
            println("Event $event: $(round(improvement, digits=2))mm $(improvement > 0 ? "(MST better)" : "(Combined better)")")
        end
    end

    if !isempty(improvements)
        n_better = sum(improvements .> 0)
        println("\nMST better in $n_better/$(length(improvements)) events")
        println("Average improvement: $(round(mean(improvements), digits=2))mm")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
```

---

## Verification Steps

After implementing all tasks:

1. **Run the test script:**
   ```bash
   cd /Users/jjgomezcadenas/Projects/Petit
   julia scripts/test_extreme_methods.jl --input scripts/xe137_selected.h5
   ```

2. **Expected output for event 12:**
   - Combined: d_start ~15mm, d_end ~1mm, total ~16mm
   - MST: d_start < 5mm, d_end < 5mm, total < 10mm
   - Energy asymmetry: One extreme ~500 keV (Bragg peak), one ~100-200 keV (track start)

3. **Success criteria:**
   - MST has lower total distance to MC endpoints
   - MST shows clear energy asymmetry between extremes
   - MST works better for majority of test events

## Notes

- MST computation is O(E log E), negligible for typical track sizes (~1000 voxels)
- Tree diameter (double BFS) is O(V) for trees
- The method is energy-agnostic - it finds geometric extremes based on topology
- Confidence is always 1.0 for MST since it's a deterministic algorithm
- The selected events file `xe137_selected.h5` should already exist in `scripts/` directory

## Test Events

Events selected for testing: **12, 22, 86, 190, 192, 292, 334, 341, 348**

These are xe137 single-electron events (2400-2500 keV) that likely contain U-shaped or curved tracks where the current algorithm fails.

## File Summary

| File | Action | Description |
|------|--------|-------------|
| `src/tracks.jl` | EDIT | Add `find_extremes_mst` function (~60 lines), update `find_track_extremes` |
| `src/Petit.jl` | EDIT | Add `find_extremes_mst` to exports |
| `scripts/test_extreme_methods.jl` | CREATE | Diagnostic script (~200 lines) |
