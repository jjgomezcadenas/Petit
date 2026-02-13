# Extreme Finding Algorithm Testing Plan

## Background

### Problem Statement
The current extreme-finding algorithm fails on U-shaped tracks. For xe137 single-electron events (e.g., event 12), both found extremes end up near the Bragg peak region instead of finding the true track start (low energy) and Bragg peak (high energy).

**Evidence from event 12:**
- Expected: One extreme at track start (low E), one at Bragg peak (high E ~500 keV)
- Actual: Both extremes near Bragg peak (d1=15.07mm, d2=1.15mm from MC endpoints, both ~500 keV)

### Root Cause Analysis
1. **Graph has shortcuts**: The radius graph (KDT) creates edges across U-shaped bends
2. **kNN doesn't fully fix it**: Even with kNN graph (fewer edges), the graph still has cycles
3. **Extreme-finding uses flawed traversal**: The `find_extremes_combined_opt` function uses BFS/Dijkstra on graphs with cycles, which can take shortcuts

### Solution Discovered
**MST + Tree Diameter** finds correct extremes:
- MST (Minimum Spanning Tree) removes all cycles, creating a tree
- Tree diameter algorithm (double BFS) is guaranteed to find the two farthest points
- Preliminary test showed MST finds extremes (703, 339) vs original (126, 588) - different endpoints!

## Test Dataset

### Selected Events
File: `xe137_selected.h5` (created from `xe137/electrons_2400_2500_15bar_100mum.next.h5`)

Events: **12, 22, 86, 190, 192, 292, 334, 341, 348**

These are xe137 single-electron events in the 2400-2500 keV range that likely contain U-shaped or curved tracks.

### File Location
```
scripts/xe137_selected.h5
```
or copy to:
```
$DATA/HD5t/itaca/xe137/xe137_selected.h5
```

## Implementation Tasks

### Task 1: Implement MST-based Extreme Finding

Add to `src/tracks.jl`:

```julia
"""
    find_extremes_mst(track::Tracks)

Find track extremes using MST + tree diameter algorithm.
Guaranteed to find the two geometrically farthest points along the track topology.

# Returns
- `(extreme1_idx, extreme2_idx, path, confidence, path_length)`
"""
function find_extremes_mst(track::Tracks)
    g = track.graph
    n = nv(g)

    # Handle trivial cases
    n == 0 && return (nothing, nothing, Int[], 0.0, 0.0)
    n == 1 && return (1, 1, [1], 1.0, 0.0)
    ne(g) == 0 && return (1, n, Int[], 0.0, 0.0)

    # Extract coordinates
    vox = track.voxels
    coords = hcat(vox.x, vox.y, vox.z)'  # 3×N matrix

    # Build distance matrix for existing edges
    distmx = zeros(n, n)
    for e in edges(g)
        i, j = src(e), dst(e)
        d = sqrt(sum((coords[:, i] .- coords[:, j]).^2))
        distmx[i, j] = d
        distmx[j, i] = d
    end

    # Compute MST using Kruskal's algorithm
    mst_edges = kruskal_mst(g, distmx)

    # Build MST graph
    mst = SimpleGraph(n)
    for e in mst_edges
        add_edge!(mst, src(e), dst(e))
    end

    # Tree diameter: double BFS
    # First BFS from arbitrary node
    dists1 = gdistances(mst, 1)
    extreme1 = argmax(dists1)

    # Second BFS from extreme1
    dists2 = gdistances(mst, extreme1)
    extreme2 = argmax(dists2)

    # Find path between extremes on MST
    path = find_path_bfs_opt(mst, extreme1, extreme2)

    # Calculate physical path length
    path_length = 0.0
    for i in 1:(length(path)-1)
        p1, p2 = path[i], path[i+1]
        path_length += sqrt(sum((coords[:, p1] .- coords[:, p2]).^2))
    end

    return (extreme1, extreme2, path, 1.0, path_length)  # confidence=1.0 for MST
end
```

**Dependencies**:
- `Graphs.jl`: `kruskal_mst`, `gdistances`, `SimpleGraph` (already imported)
- `find_path_bfs_opt`: already exists in tracks.jl

### Task 2: Add MST Method to find_track_extremes

Update `find_track_extremes` to support MST method:

```julia
function find_track_extremes(track::Tracks; use_optimized::Bool=true, method::Symbol=:combined)
    if method == :mst
        return find_extremes_mst(track)
    elseif use_optimized
        return find_track_extremes_opt(track)
    else
        return find_track_extremes_improved(track)
    end
end
```

### Task 3: Export New Function

Add to `src/Petit.jl`:
```julia
export find_extremes_mst
```

## Testing Protocol

### Test 1: Diagnostic Comparison (Manual)

Run this diagnostic script to compare all methods on selected events:

```julia
# File: scripts/test_extreme_methods.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DataFrames
using Graphs
using Printf

include(joinpath(@__DIR__, "..", "src", "Petit.jl"))
import .Petit

const SELECTED_FILE = "xe137_selected.h5"  # or full path
const EVENTS = [12, 22, 86, 190, 192, 292, 334, 341, 348]

function test_event(hitsdf, event_id; σt=1.012, σl=0.0, voxel_size=2.024,
                    mcvox_size=2.024, max_distance=3.035)

    event_df = Petit.get_event(hitsdf, event_id)
    nrow(event_df) == 0 && return nothing

    # MC path for ground truth
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
    energies = track.voxels.energy .* 1000  # keV

    results = Dict{String, Any}()
    results["event_id"] = event_id
    results["n_voxels"] = nrow(track.voxels)
    results["n_edges"] = ne(track.graph)

    # Test each method
    for (name, method) in [
        ("combined", :combined),
        ("mst", :mst)
    ]
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
        catch e
            results["$(name)_error"] = string(e)
        end
    end

    return results
end

function main()
    println("Loading data...")
    hitsdf = Petit.load_data(SELECTED_FILE, ".")  # Adjust path as needed

    println("\n" * "="^80)
    println("EXTREME FINDING METHOD COMPARISON")
    println("="^80)

    all_results = []

    for event_id in EVENTS
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

                println("  $method: extremes=($e1,$e2), E=($(round(E1,digits=1)),$(round(E2,digits=1)))keV")
                println("          d_start=$(round(d_start,digits=2))mm, d_end=$(round(d_end,digits=2))mm, total=$(round(total,digits=2))mm")
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
            println("$method: mean_total_d = $(round(mean(totals), digits=2))mm, median = $(round(median(totals), digits=2))mm")
        end
    end
end

main()
```

### Test 2: Full Pipeline Test

After implementing MST method, run the analysis script with both methods:

```bash
# Combined method (current)
julia itaca_single_track_analysis.jl -i 12 -e 348 \
    --input xe137_selected.h5 \
    --lmin 100 --lmax 100 \
    --print_level verbose \
    -o results_combined

# MST method (new) - requires adding --extreme-method flag
julia itaca_single_track_analysis.jl -i 12 -e 348 \
    --input xe137_selected.h5 \
    --lmin 100 --lmax 100 \
    --extreme-method mst \
    --print_level verbose \
    -o results_mst
```

### Test 3: Compare Results

Create comparison script:

```julia
# Compare CSV results
using CSV, DataFrames

df_combined = CSV.read("results_combined/analysis_results.csv", DataFrame)
df_mst = CSV.read("results_mst/analysis_results.csv", DataFrame)

comparison = DataFrame(
    event = df_combined.event,
    d1_combined = df_combined.d1_mm,
    d2_combined = df_combined.d2_mm,
    d1_mst = df_mst.d1_mm,
    d2_mst = df_mst.d2_mm,
    Eb1_combined = df_combined.Eb1_keV,
    Eb2_combined = df_combined.Eb2_keV,
    Eb1_mst = df_mst.Eb1_keV,
    Eb2_mst = df_mst.Eb2_keV
)

println(comparison)

# Check improvement
combined_total = df_combined.d1_mm .+ df_combined.d2_mm
mst_total = df_mst.d1_mm .+ df_mst.d2_mm
println("\nImprovement (negative = MST better):")
println(mst_total .- combined_total)
```

## Success Criteria

### For Single Electron (xe137) Events:

1. **Distance to MC endpoints**: MST should have smaller total distance (d_start + d_end)
2. **Energy asymmetry**: One extreme should have high energy (~400-600 keV, Bragg peak), one should have lower energy (track start)
3. **Blob analysis**: After using MST extremes, blob energies should be asymmetric (Eb1 >> Eb2 or vice versa)

### Specific Targets for Event 12:
- Current: d1=15.07mm, d2=1.15mm (wrong - both near Bragg peak)
- Expected with MST: d_start < 5mm, d_end < 5mm (both close to MC endpoints)
- Blob energies should be clearly asymmetric (one ~500 keV, one ~100-200 keV)

## Expected Outcomes

| Metric | Current (combined) | Expected (MST) |
|--------|-------------------|----------------|
| d_start (mm) | ~15 | < 5 |
| d_end (mm) | ~1 | < 5 |
| Total distance | ~16 | < 10 |
| Eb1 (keV) | ~500 | ~500 |
| Eb2 (keV) | ~460 | ~100-200 |
| Asymmetry | ~0.06 | > 0.3 |

## Implementation Order

1. **Add `find_extremes_mst` function** to `src/tracks.jl`
2. **Update `find_track_extremes`** to accept `method` parameter
3. **Export function** in `src/Petit.jl`
4. **Run diagnostic test** (Test 1) on selected events
5. **Verify improvement** in distance metrics
6. **Optionally**: Add `--extreme-method` flag to CLI scripts for A/B testing

## Files to Modify

| File | Changes |
|------|---------|
| `src/tracks.jl` | Add `find_extremes_mst`, update `find_track_extremes` |
| `src/Petit.jl` | Export `find_extremes_mst` |
| `scripts/test_extreme_methods.jl` | New diagnostic script |

## Notes

- The MST approach uses `Graphs.jl`'s `kruskal_mst` which is already available
- MST computation is O(E log E) which is negligible for typical track sizes (~1000 voxels)
- Tree diameter (double BFS) is O(V + E) = O(V) for trees
- The method is **energy-agnostic** - it finds geometric extremes, which should correspond to physical endpoints regardless of event type
