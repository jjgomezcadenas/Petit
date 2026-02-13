# kNN Track Building: Implementation Plan

## Summary

Add `build_tracks_knn` as an alternative to `build_tracks_kdtree` to fix U-shaped track extreme finding failures. This is an **upstream fix** that prevents shortcut edges from being created, rather than pruning them afterward.

## Comparison: kNN vs MST Approaches

| Aspect | kNN (this plan) | MST (alternative) |
|--------|-----------------|-------------------|
| **Where** | Graph construction | Post-processing |
| **Strategy** | Prevent shortcuts | Prune shortcuts |
| **Complexity** | O(N log N) for kNN queries | O(E log E) for MST |
| **Edge count** | ~N*k/2 (sparse, predictable) | N-1 (tree) |
| **Dependencies** | NearestNeighbors.jl (already have) | Graphs.jl (already have) |
| **Risk** | New graph topology | Changes existing graph |

**Recommendation**: kNN is simpler and attacks root cause directly.

---

## Implementation

### Phase 1: Core Function `build_tracks_knn`

Add to `src/tracks.jl`:

```julia
"""
    build_tracks_knn(event_data::DataFrame; k=10, max_distance=Inf,
                     energy_threshold=0.0, diffusion=DiffusionParams(),
                     mutual=false)

Build tracks using k-Nearest Neighbor graph instead of radius graph.

# Arguments
- `event_data::DataFrame`: Voxelized event data with x, y, z, energy columns
- `k::Int`: Number of nearest neighbors per voxel (default: 10)
- `max_distance::Float64`: Maximum edge length allowed (default: Inf, no limit)
- `energy_threshold::Float64`: Minimum voxel energy in MeV (default: 0.0)
- `diffusion::DiffusionParams`: Diffusion parameters for the tracks
- `mutual::Bool`: If true, keep only mutual kNN edges (default: false)

# Returns
- `Vector{Tracks}`: Tracks sorted by total energy (highest first)

# Notes
kNN graph prevents shortcut edges across U-shaped bends by limiting each
voxel to its k closest neighbors, rather than all neighbors within a radius.
"""
function build_tracks_knn(event_data::DataFrame;
                          k::Int = 10,
                          max_distance::Float64 = Inf,
                          energy_threshold::Float64 = 0.0,
                          diffusion::DiffusionParams = DiffusionParams(),
                          mutual::Bool = false)

    n_voxels = nrow(event_data)
    n_voxels == 0 && return Tracks[]

    # Filter by energy
    valid_indices = findall(e -> e >= energy_threshold, event_data.energy)
    isempty(valid_indices) && return Tracks[]

    filtered_data = event_data[valid_indices, :]
    n_filtered = nrow(filtered_data)

    # Trivial cases
    n_filtered == 1 && return [Tracks(filtered_data, SimpleGraph(1), [collect(1:1)], diffusion)]

    # Adjust k if we have fewer points
    k_actual = min(k, n_filtered - 1)
    k_actual < 1 && return [Tracks(filtered_data, SimpleGraph(n_filtered), [collect(1:n_filtered)], diffusion)]

    # Build KDTree
    coords = _coords_from_df(filtered_data)
    tree = KDTree(coords)

    # Build kNN graph
    g = SimpleGraph(n_filtered)

    if mutual
        # Mutual kNN: keep edge only if both endpoints have each other in their kNN
        # First, collect all kNN sets
        knn_sets = Vector{Set{Int}}(undef, n_filtered)
        @inbounds for i in 1:n_filtered
            idxs, dists = knn(tree, coords[:, i], k_actual + 1, true)
            # Filter by max_distance and exclude self
            neighbors = Set{Int}()
            for (idx, d) in zip(idxs, dists)
                idx == i && continue
                if d <= max_distance
                    push!(neighbors, idx)
                end
            end
            knn_sets[i] = neighbors
        end

        # Add edge only if mutual
        @inbounds for i in 1:n_filtered
            for j in knn_sets[i]
                j <= i && continue  # Avoid duplicates
                if i in knn_sets[j]
                    add_edge!(g, i, j)
                end
            end
        end
    else
        # Standard kNN: add edge if j is in kNN(i) OR i is in kNN(j)
        @inbounds for i in 1:n_filtered
            idxs, dists = knn(tree, coords[:, i], k_actual + 1, true)
            for (idx, d) in zip(idxs, dists)
                idx == i && continue
                idx < i && continue  # Will be added when processing idx
                if d <= max_distance
                    add_edge!(g, i, idx)
                end
            end
        end
    end

    # Connected components → tracks
    components = connected_components(g)
    tracks = Tracks[]

    for comp in components
        isempty(comp) && continue
        comp_data = filtered_data[comp, :]
        subg, _ = induced_subgraph(g, comp)
        push!(tracks, Tracks(comp_data, subg, [comp], diffusion))
    end

    return tracks
end
```

### Phase 2: Update `make_tracks` to Support kNN

Modify `make_tracks` in `src/tracks.jl`:

```julia
"""
    make_tracks(event_data; max_distance_mm=1.0, energy_threshold_kev=0.0,
                diffusion=DiffusionParams(), method="KDT", k=10)

Build and sort tracks by total energy (descending).

# Arguments
- `event_data::DataFrame`: Voxelized event data
- `max_distance_mm::Float64`: Max voxel connection distance in mm
- `energy_threshold_kev::Float64`: Min voxel energy in keV
- `diffusion::DiffusionParams`: Diffusion parameters for the tracks
- `method::String`: Graph construction method:
  - "KDT": Radius graph (all neighbors within max_distance)
  - "kNN": k-Nearest Neighbor graph (k closest neighbors)
  - "kNN_mutual": Mutual kNN (edge only if both are in each other's kNN)
- `k::Int`: Number of neighbors for kNN methods (default: 10)

# Returns
- `Vector{Tracks}`: Tracks sorted by total energy (highest first)
"""
function make_tracks(event_data::DataFrame;
                    max_distance_mm::Float64=1.0,
                    energy_threshold_kev::Float64=0.0,
                    diffusion::DiffusionParams=DiffusionParams(),
                    method::String="KDT",
                    k::Int=10)

    energy_threshold = energy_threshold_kev * 1e-3  # Convert keV to MeV

    if method == "KDT"
        tracks = build_tracks_kdtree(event_data;
                                     max_distance=max_distance_mm,
                                     energy_threshold=energy_threshold,
                                     diffusion=diffusion)
    elseif method == "kNN"
        tracks = build_tracks_knn(event_data;
                                  k=k,
                                  max_distance=max_distance_mm,
                                  energy_threshold=energy_threshold,
                                  diffusion=diffusion,
                                  mutual=false)
    elseif method == "kNN_mutual"
        tracks = build_tracks_knn(event_data;
                                  k=k,
                                  max_distance=max_distance_mm,
                                  energy_threshold=energy_threshold,
                                  diffusion=diffusion,
                                  mutual=true)
    else
        error("Unknown method: $method. Use 'KDT', 'kNN', or 'kNN_mutual'")
    end

    # Sort by total energy (descending)
    if length(tracks) > 0
        track_energies = [sum(track.voxels.energy) for track in tracks]
        sorted_indices = sortperm(track_energies, rev=true)
        tracks = tracks[sorted_indices]
    end

    return tracks
end
```

### Phase 3: Add CLI Options to Scripts

Update `scripts/itaca_single_track_analysis.jl` and `scripts/itaca_track_reco_mt.jl`:

```julia
# Add to argument parser:
"--graph-method", "-g"
    help = "Graph construction method: KDT (radius), kNN, kNN_mutual"
    arg_type = String
    default = "KDT"

"--knn-k"
    help = "Number of neighbors for kNN methods"
    arg_type = Int
    default = 10
```

Update call to `make_tracks`:
```julia
tracks = make_tracks(voxels;
                     max_distance_mm=max_distance,
                     energy_threshold_kev=energy_threshold,
                     diffusion=diffusion,
                     method=args["graph-method"],
                     k=args["knn-k"])
```

---

## Parameter Guidelines

### Choosing k

| Track density | Recommended k | Notes |
|---------------|---------------|-------|
| Sparse (large voxels) | 6-8 | Avoid disconnected components |
| Normal | 10 | Good default |
| Dense (small voxels) | 12-15 | More neighbors needed for connectivity |

**Rule of thumb**: k ≈ 2× expected coordination number in a chain-like track

### When to Use max_distance

- Set `max_distance = voxel_size * 2` as safety against very sparse regions
- Leave as `Inf` if voxel density is uniform

### When to Use Mutual kNN

Use `mutual=true` if:
- Standard kNN still shows shortcuts
- Track has very dense regions mixed with sparse regions
- Willing to accept potentially more disconnected components

---

## Testing Plan

### Test 1: xe137 Event 12 (U-shaped failure case)

```bash
julia scripts/itaca_single_track_analysis.jl \
    -i 12 -e 12 \
    --input xe137/electrons_2400_2500_15bar_100mum.next.h5 \
    --lmin 100 --lmax 100 \
    --graph-method kNN --knn-k 10 \
    --iterative --verbose
```

**Expected**: d1 should be large (track start), d2 small (Bragg peak), asymmetric blob energies.

### Test 2: bb0nu Events (Regression check)

```bash
julia scripts/itaca_single_track_analysis.jl \
    -i 1 -e 10 \
    --input bb0nu/bb_2450_2470_15bar_100mum.next.h5 \
    --lmin 100 --lmax 100 \
    --graph-method kNN --knn-k 10 \
    --iterative
```

**Expected**: Both blob energies high (~500 keV), similar to KDT method.

### Test 3: Compare Methods

Run same events with all three methods:
1. `--graph-method KDT` (current baseline)
2. `--graph-method kNN`
3. `--graph-method kNN_mutual`

Compare:
- Blob energies (Eb1, Eb2)
- Distances to MC endpoints (d1, d2)
- Track length
- Number of connected components

---

## Success Criteria

1. **xe137 event 12 fixed**: One extreme at track start, one at Bragg peak
2. **bb0nu unchanged**: Both extremes at Bragg peaks
3. **No disconnected tracks**: kNN with reasonable k keeps track connected
4. **Performance**: Similar or better than radius graph

---

## Files to Modify

| File | Changes |
|------|---------|
| `src/tracks.jl` | Add `build_tracks_knn`, update `make_tracks` |
| `src/Petit.jl` | Export `build_tracks_knn` |
| `scripts/itaca_single_track_analysis.jl` | Add `--graph-method`, `--knn-k` options |
| `scripts/itaca_track_reco_mt.jl` | Add `--graph-method`, `--knn-k` options |

---

## Implementation Order

1. Add `build_tracks_knn` function
2. Update `make_tracks` to support method selection
3. Test manually with xe137 event 12
4. If successful, add CLI options to scripts
5. Run regression tests on bb0nu
6. Document recommended parameters

---

## Appendix: Why kNN Prevents Shortcuts

Consider a U-shaped track:

```
    A---B---C
    |       |
    |       |
    D---E---F (Bragg peak)
```

**Radius graph** (max_distance covers the U width):
- A connects to: B, D, and possibly E or F (shortcuts!)
- Creates O(N²) edges in dense regions

**kNN graph** (k=4):
- A's 4 nearest: B, D, and 2 others along the spine
- A does NOT connect to F (not among 4 closest)
- Creates exactly N*k/2 edges

The shortcuts (A-F, B-E, C-D) require jumping across the U, which means larger distances. With kNN, these longer edges are excluded in favor of spine neighbors.
