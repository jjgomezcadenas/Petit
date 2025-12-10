# Plan: Refactor track_reco_mt.jl - Integrate KDE, Remove Double Smoothing

## Motivation

Current pipeline has double smoothing:
1. Central path smoothing with `filter_radius` (spatial)
2. KDE smoothing with `bandwidth` (1D)

This is redundant. We should do projection onto raw path, then KDE only.

## Goals

1. Remove central path smoothing (eliminate `filter_radius`)
2. Project RECO voxels onto raw path → compute arc-length `s` for each
3. Apply KDE to both RECO and MC voxels
4. Store KDE results in HDF5 for downstream analysis

---

## Implementation Steps

### Step 1: Update RecoResult struct

```julia
struct RecoResult
    event_id::Int
    track::Petit.Tracks
    path::DataFrame           # Raw path (x, y, z, s) - NOT smoothed
    track_length::Float64
    confidence::Float64
    mc_path::DataFrame        # MC path (x, y, z, energy, s)
    reco_s::Vector{Float64}   # Arc-length of each reco voxel
    kde_s::Vector{Float64}    # KDE evaluation points
    reco_kde_f::Vector{Float64}  # RECO energy density
    mc_kde_f::Vector{Float64}    # MC energy density
    kde_bandwidth::Float64    # Bandwidth used
end
```

### Step 2: Update process_single_event()

Replace:
```julia
# OLD: Reconstruct central path (smoothed)
central_path = Petit.reconstruct_central_path(track, walk_result.path_indices;
                                               filter_radius=filter_radius)
```

With:
```julia
# NEW: Get raw path (no smoothing)
path = Petit.get_raw_path(track, walk_result.path_indices)

# Compute arc-length for path
path.s = compute_arc_length(path)

# Project reco voxels onto path
reco_s, reco_E = project_voxels_to_path(track.voxels, path, path.s)

# Compute KDE for RECO
kde_s = range(0, path.s[end], length=n_kde_eval)
reco_kde_f, _ = energy_weighted_kde(reco_s, reco_E, kde_s; bandwidth=kde_bandwidth)

# Compute KDE for MC (mc_path already has s and energy)
mc_kde_f, _ = energy_weighted_kde(mc_path.s, mc_path.energy, kde_s; bandwidth=kde_bandwidth)
```

### Step 3: Add get_raw_path() function to Petit

In `src/tracks.jl` or similar:
```julia
function get_raw_path(track::Tracks, path_indices::Vector{Int})
    voxels = track.voxels
    path_df = DataFrame(
        x = voxels.x[path_indices],
        y = voxels.y[path_indices],
        z = voxels.z[path_indices]
    )
    return path_df
end
```

### Step 4: Update save_reco_results_to_hdf5()

Remove:
- `save_central_path_to_hdf5()` call

Add:
```julia
# Save raw path
g["path_data"] = Matrix(result.path)
g["path_columns"] = String.(names(result.path))

# Save projected voxel arc-lengths
g["reco_s"] = result.reco_s

# Save KDE results
g["kde_s"] = result.kde_s
g["reco_kde_f"] = result.reco_kde_f
g["mc_kde_f"] = result.mc_kde_f
g["kde_bandwidth"] = result.kde_bandwidth
```

### Step 5: Update CLI parameters

Remove:
- `--filter-radius` parameter

Add:
- `--kde-bandwidth=X` (default: 5.0 mm)
- `--n-kde-eval=N` (default: 200)

### Step 6: Update event_loop_reco_mt()

- Remove `filter_radius` from parameters and config printout
- Add `kde_bandwidth` and `n_kde_eval`
- Update metadata dict

### Step 7: Update main()

- Remove `filter_radius` variable and parsing
- Add `kde_bandwidth` and `n_kde_eval` variables and parsing

---

## Files to Modify

| File | Changes |
|------|---------|
| `scripts/track_reco_mt.jl` | All steps above |
| `src/tracks.jl` | Add `get_raw_path()` function |
| `src/Petit.jl` | Export `get_raw_path` |

## Files to Update (documentation)

| File | Changes |
|------|---------|
| `docs/track_reco_mt_summary.md` | Update structure, parameters |

---

## New HDF5 Structure Per Track

### Metadata
- `event_id`, `track_length`, `confidence`
- `kde_bandwidth`, `n_kde_eval` (NEW)

### Track Data (unchanged)
- `voxels`, `voxel_columns`
- `graph_edges`, `n_vertices`, `components`

### Path Data
- `path_data` (raw path: x, y, z, s) - replaces `central_path_data`
- `path_columns`
- `mc_path_data`, `mc_path_columns` (unchanged)

### Projected & KDE (NEW)
- `reco_s` - arc-length of each reco voxel
- `kde_s` - evaluation grid
- `reco_kde_f` - RECO energy density
- `mc_kde_f` - MC energy density

---

## Removed

- `filter_radius` parameter
- `central_path_data` (smoothed path)
- `reconstruct_central_path()` call

## Added

- `kde_bandwidth` parameter (default: 5.0 mm)
- `n_kde_eval` parameter (default: 200)
- `get_raw_path()` function
- `reco_s`, `kde_s`, `reco_kde_f`, `mc_kde_f` in output

---

## Testing

After implementation:
1. Run `track_reco_mt.jl` on test data
2. Verify HDF5 contains new fields
3. Update `test_kde_track.jl` to use pre-computed KDE (or keep for comparison)
4. Compare old vs new results visually
