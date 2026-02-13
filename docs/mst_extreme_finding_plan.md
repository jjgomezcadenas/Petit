# MST-Based Track Extreme Finding: Implementation Plan

## Problem Statement

### Current Behavior
The track extreme-finding algorithm fails for U-shaped tracks. For example, event 12 of xe137 (single electron, 2400-2500 keV) shows:
- **Expected**: One extreme at track start (low energy), one at Bragg peak (high energy ~500 keV)
- **Actual**: Both extremes found near Bragg peak region (d1=15.07mm, d2=1.15mm from MC end, both ~500 keV)

### Root Cause Analysis

The problem originates in `build_tracks_kdtree` (`src/tracks.jl:168-220`):

```julia
function build_tracks_kdtree(event_data::DataFrame; max_distance::Float64 = 1.5, ...)
    # ...
    g = SimpleGraph(n_filtered)
    r = max_distance

    # Build graph using neighbors in range
    @inbounds for i in 1:n_filtered
        idxs = inrange(tree, coords[:, i], r, true)  # ALL neighbors within r
        for j in idxs
            j == i && continue
            j < i && continue
            add_edge!(g, i, j)  # Creates ALL possible edges
        end
    end
    # ...
end
```

This creates a **full connectivity graph** where every voxel is connected to all neighbors within `max_distance`. For U-shaped tracks, this creates **shortcut edges** across the bend:

```
U-shaped track topology:

    A -------- B              Physical track: A → B → D → C
    |          |              (electron starts at A, Bragg peak at C or D)
    |    X     |
    |          |              Full connectivity graph adds:
    C -------- D              - Diagonal shortcuts: A-D, B-C
                              - These shortcuts are WRONG
```

### Why Shortcuts Break Extreme Finding

All current extreme-finding methods (`find_extremes_topology_opt`, `find_extremes_spatial_opt`, etc.) use BFS/Dijkstra on this graph. With shortcuts:

1. **Path distortion**: Shortest path A→C is no longer A→B→D→C (length=3), but A→C direct (length=1 if shortcut exists)
2. **Wrong diameter**: Graph diameter is reduced, making distant points appear close
3. **False extremes**: When searching for "farthest" points, the algorithm may find points on the same arm of the U (both near Bragg peak) instead of true endpoints

For event 12:
- True extremes: track start (low E) and Bragg peak (high E)
- Found extremes: two points near Bragg peak (both ~500 keV)
- Cause: shortcuts made the track start appear closer than points on the Bragg peak arm

## Proposed Solution: Minimum Spanning Tree (MST)

### Why MST Works

**MST properties:**
1. Connects all vertices (no isolated voxels)
2. Uses minimum total edge weight
3. Contains no cycles (it's a tree)
4. **Critically: Removes redundant edges, including shortcuts**

With Euclidean distance as edge weight:
- Spine edges (along track): ~1-2 voxel spacing, short
- Shortcut edges (across bend): typically √2× to 2× longer

MST algorithm (Kruskal's):
1. Sort edges by weight (shortest first)
2. Add edges in order, skipping any that would create a cycle
3. Stop when all vertices connected

Result: Short spine edges added first → by the time shortcuts are considered, they would create cycles → **shortcuts rejected**.

```
MST of U-shaped track:

    A -------- B              MST preserves: A-B, B-D, D-C
    |          |              (the track spine)
    |          |
    |          |              MST removes: A-D, B-C
    C -------- D              (the shortcuts)
```

### Tree Diameter Algorithm

Once we have the MST (a tree), finding extremes becomes trivial:

**Algorithm:**
1. Start BFS from any vertex → find farthest vertex E1
2. Start BFS from E1 → find farthest vertex E2
3. E1 and E2 are the tree diameter endpoints (guaranteed correct for trees)

**Complexity:** O(V + E) - two BFS traversals

This works because in a tree, there's exactly one path between any two vertices. The diameter endpoints are the two vertices with maximum path length.

### Why MST Works for Both Event Types

**Single electron (xe137):**
- Track shape: start → curve → Bragg peak
- MST preserves spine, removes shortcuts across curve
- Diameter finds: start (low E) ↔ Bragg peak (high E) ✓

**Double-beta (bb0nu):**
- Track shape: Bragg peak 1 ← shared vertex → Bragg peak 2
- Track is relatively straight, few/no shortcuts to begin with
- Diameter finds: Bragg peak 1 ↔ Bragg peak 2 ✓

**Key insight:** The algorithm is energy-agnostic. It finds the geometrically farthest points along the track topology, which correspond to physical endpoints regardless of event type.

## Implementation Plan

### Phase 1: Core MST Function

Create new function `compute_mst_graph` in `src/tracks.jl`.

**Note**: `Graphs.jl` (already in Project.toml) has `kruskal_mst` built-in with distance matrix support:
```julia
kruskal_mst(g, distmx; minimize=true)
```

No additional dependencies needed.

```julia
using Graphs: kruskal_mst, SimpleGraph, add_edge!, nv, ne, src, dst, edges

"""
    compute_mst_graph(g::SimpleGraph, coords::Matrix{Float64}) -> SimpleGraph

Compute Minimum Spanning Tree from connectivity graph using Euclidean distances.

# Arguments
- `g::SimpleGraph`: Full connectivity graph from build_tracks_kdtree
- `coords::Matrix{Float64}`: 3×N coordinate matrix (x, y, z for each vertex)

# Returns
- `SimpleGraph`: MST as unweighted graph (subset of edges from g)
"""
function compute_mst_graph(g::SimpleGraph, coords::Matrix{Float64})
    n = nv(g)
    n <= 1 && return g  # Trivial case
    ne(g) == 0 && return g  # No edges

    # Build distance matrix (only for existing edges, others stay 0)
    distmx = zeros(n, n)
    for e in edges(g)
        i, j = src(e), dst(e)
        d = sqrt(sum((coords[:, i] .- coords[:, j]).^2))
        distmx[i, j] = d
        distmx[j, i] = d
    end

    # Compute MST using Kruskal's algorithm
    mst_edges = kruskal_mst(g, distmx)

    # Build unweighted MST graph
    mst_graph = SimpleGraph(n)
    for e in mst_edges
        add_edge!(mst_graph, src(e), dst(e))
    end

    return mst_graph
end
```

### Phase 2: Tree Diameter Function

Create `find_tree_diameter` function:

```julia
"""
    find_tree_diameter(tree::SimpleGraph) -> (Int, Int, Vector{Int}, Float64)

Find the diameter endpoints of a tree using double BFS.

# Returns
- `(extreme1, extreme2, path, confidence)`: Endpoint indices, path between them, confidence=1.0
"""
function find_tree_diameter(tree::SimpleGraph, coords::Matrix{Float64})
    n = nv(tree)
    n == 0 && return (nothing, nothing, Int[], 0.0)
    n == 1 && return (1, 1, [1], 1.0)

    # First BFS: find farthest from vertex 1
    dists1 = gdistances(tree, 1)
    extreme1 = argmax(dists1)

    # Second BFS: find farthest from extreme1
    dists2 = gdistances(tree, extreme1)
    extreme2 = argmax(dists2)

    # Find path between extremes
    path = find_path_bfs_opt(tree, extreme1, extreme2)

    # Calculate physical path length
    path_length = calculate_path_length_from_coords(coords, path)

    return (extreme1, extreme2, path, 1.0, path_length)
end
```

### Phase 3: Integration with Existing Code

#### Option A: New MST-based extreme finder (recommended)

Add `find_extremes_mst` as a new method:

```julia
"""
    find_extremes_mst(track::Tracks, coords::TrackCoords) -> Tuple

MST-based extreme finding. Removes graph shortcuts by computing MST first.
Most robust for curved/U-shaped tracks.
"""
function find_extremes_mst(track::Tracks, coords::TrackCoords)
    g = track.graph
    n = nv(g)

    # Extract raw coordinates
    coord_matrix = hcat(coords.x, coords.y, coords.z)'  # 3×N matrix

    # Compute MST
    mst = compute_mst_graph(g, coord_matrix)

    # Find diameter of MST
    return find_tree_diameter(mst, coord_matrix)
end
```

#### Option B: Replace graph in Tracks struct

Modify `build_tracks_kdtree` to store MST instead of full graph:

```julia
# In build_tracks_kdtree, after building g:
coord_matrix = Matrix(coords)  # 3×N
mst = compute_mst_graph(g, coord_matrix)

# Use mst instead of g for track
push!(tracks, Tracks(comp_data, mst, [comp], diffusion))
```

**Tradeoff**: Option A is safer (doesn't break existing code), Option B is cleaner (all methods benefit automatically).

**Recommendation**: Start with Option A, test thoroughly, then consider Option B.

### Phase 4: Update `find_extremes_combined_opt`

Integrate MST method into the combined strategy:

```julia
function find_extremes_combined_opt(track::Tracks, coords::TrackCoords;
                                     use_energy_weighting::Bool=true,
                                     use_mst::Bool=true)  # New parameter
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0)
    end

    # Check track density - high density suggests potential shortcuts
    avg_degree = 2 * ne(g) / n_vertices
    is_dense = avg_degree > 6.0

    # For dense tracks, MST is essential
    if is_dense && use_mst
        mst_result = find_extremes_mst(track, coords)
        if mst_result[4] >= 0.9  # High confidence
            return (mst_result[1], mst_result[2], mst_result[3], mst_result[4])
        end
    end

    # Continue with existing logic as fallback...
    # (existing code unchanged)
end
```

### Phase 5: Dependencies

**No new dependencies required.**

`Graphs.jl` (already in Project.toml v1.13.0) provides:
- `kruskal_mst(g, distmx)` - MST with distance matrix
- `gdistances(g, source)` - BFS distances for tree diameter
- `SimpleGraph`, `add_edge!`, `nv`, `ne`, `src`, `dst`, `edges`

## Testing Plan

### Unit Tests

1. **Linear track**: MST should equal original graph
2. **U-shaped track**: MST should remove diagonal shortcuts
3. **Star/branching track**: MST should preserve branches, remove redundant edges
4. **Dense cloud**: MST should create tree structure

### Integration Tests

1. **xe137 event 12**: Should find correct extremes (start and Bragg peak)
2. **bb0nu events**: Should still find both Bragg peaks correctly
3. **Random drift events**: Should work across drift lengths

### Validation Metrics

- `d1_mm` and `d2_mm`: Distance from MC endpoints
- `Eb1_keV` and `Eb2_keV`: Energy at found extremes
- For single electron: one extreme should have low E, one high E
- For bb0nu: both extremes should have high E

## Risk Assessment

### Low Risk
- MST computation is well-understood, O(E log E)
- Tree diameter is guaranteed correct for trees
- New function doesn't break existing code (Option A)

### Medium Risk
- Performance: MST adds ~O(E log E) overhead per track
- Edge cases: disconnected components (shouldn't happen with current pipeline)

### Mitigation
- Keep original methods as fallback
- Add `use_mst` parameter to enable/disable
- Profile performance on large events

## Implementation Order

1. **Add `SimpleWeightedGraphs` dependency** (if needed)
2. **Implement `compute_mst_graph`** with unit tests
3. **Implement `find_tree_diameter`** with unit tests
4. **Implement `find_extremes_mst`** integrated with TrackCoords
5. **Add to `find_extremes_combined_opt`** as primary method for dense tracks
6. **Test on xe137 event 12** - verify fix
7. **Test on bb0nu events** - verify no regression
8. **Optional: Implement Option B** (store MST in Tracks struct)

## Success Criteria

1. xe137 event 12 correctly identifies extremes (one low E, one high E)
2. bb0nu events maintain current accuracy (both high E extremes)
3. No performance regression >20% on typical events
4. All existing tests pass

## Files to Modify

| File | Changes |
|------|---------|
| `src/tracks.jl` | Add `compute_mst_graph`, `find_tree_diameter`, `find_extremes_mst`, update `find_extremes_combined_opt` |
| `src/Petit.jl` | Export new functions (optional) |

No changes to `Project.toml` - all required functions are in `Graphs.jl`.

## Appendix: Algorithm Complexity

| Operation | Complexity | Notes |
|-----------|------------|-------|
| Build full graph | O(N × k) | k = avg neighbors in range |
| Build weighted graph | O(E) | E = number of edges |
| Kruskal's MST | O(E log E) | Sorting edges |
| Tree diameter (2×BFS) | O(V + E_mst) | E_mst = V-1 for tree |
| **Total additional** | **O(E log E)** | Dominated by MST |

For typical tracks with ~1000 voxels and ~5000 edges, this is negligible.
