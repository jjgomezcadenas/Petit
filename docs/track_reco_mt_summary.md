# Summary: `track_reco_mt.jl`

## Purpose

Multi-threaded track reconstruction script that processes MC events, reconstructs tracks, projects voxels onto paths, computes KDE energy profiles for both RECO and MC, and saves results to HDF5 files.

---

## Processing Pipeline

```text
Input HDF5 (MC hits) → Diffusion → Voxelization → Track Making → Path + KDE → Output HDF5
                    ↘ MC Path → MC KDE ↗
```

### Per-Event Processing

1. **Compute MC path** - Voxelize primary particle hits, compute arc-length
2. **Transform hits** - Convert energy to electron counts, remove MC-specific columns
3. **Diffuse** - Apply transverse (σt) and longitudinal (σl) diffusion
4. **Voxelize** - Group hits into voxels of configurable size
5. **Make tracks** - Build graph, find connected components
6. **Filter** - Keep only single-track events
7. **Walk track** - Find extremes via graph traversal
8. **Get raw path** - Extract path points (NO smoothing)
9. **Project voxels** - Map each reco voxel to arc-length s on path
10. **Compute RECO KDE** - Energy-weighted KDE along path
11. **Compute MC KDE** - Energy-weighted KDE along MC path

---

## Output Structure (HDF5)

For each reconstructed track:

### Metadata

| Field | Description |
|-------|-------------|
| `event_id` | Original event ID |
| `track_length` | Total path length (mm) |
| `confidence` | Path reconstruction confidence |
| `kde_bandwidth` | RECO KDE bandwidth used (mm) |
| `mc_kde_bandwidth` | MC KDE bandwidth used (mm) |

### Track Structure

| Field | Description |
|-------|-------------|
| `voxels` | Voxel DataFrame (x, y, z, energy, electrons) |
| `graph_edges` | Track graph connectivity |
| `n_vertices` | Number of graph vertices |
| `components` | Connected components |

### Path Data

| Field | Description |
|-------|-------------|
| `path_data` | Raw path: x, y, z, s (unsmoothed) |
| `mc_path_data` | MC truth path: x, y, z, energy, s |

### Projected Voxels & KDE

| Field | Description |
|-------|-------------|
| `reco_s` | Arc-length of each reco voxel on path |
| `kde_s` | KDE evaluation points (shared grid) |
| `reco_kde_f` | RECO energy density f(s) |
| `mc_kde_f` | MC energy density f(s) |

**Note:** MC truth extremes available as `mc_path_data[1,:]` and `mc_path_data[end,:]`.

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `σt` | -1 (auto) | Transverse diffusion (mm). If <0, computed from ion formula |
| `σl` | 0.0 | Longitudinal diffusion (mm) |
| `voxel_scale` | 1.0 | Voxel size = σt × voxel_scale |
| `voxel_distance_scale` | 1.5 | Max connection distance = voxel_size × scale |
| `energy_threshold_ions` | 10 | Minimum voxel energy |
| `ldrft` | 100 cm | Drift length for diffusion calculation |
| `mcvox` | 1.0 mm | MC path voxel size |
| `kde_bandwidth` | 5.0 mm | RECO KDE smoothing bandwidth |
| `mc_kde_bandwidth` | 1.0 mm | MC KDE smoothing bandwidth |
| `n_kde_eval` | 200 | Number of KDE evaluation points |

---

## Multi-Threading

- Splits event range across `nthreads` threads
- Each thread processes independent event subset
- Each thread writes its own output file: `{output_base}_th_{i}.h5`

---

## Usage

```bash
julia -t 8 track_reco_mt.jl /path/to/data/ input.h5 output_base \
    --nthreads=8 --ievt=1 --levt=1000 --kde-bandwidth=5.0
```

### Command Line Options

```text
Required arguments:
  cmdir           Directory containing the input file
  input_file      Name of the HDF5 input file
  output_base     Base name for output files (no .h5 extension)

Optional arguments:
  --ievt=N                   First event to process (default: 1)
  --levt=N                   Last event to process (default: -1, all)
  --nthreads=N               Number of threads to use (default: 1)
  --ldrft=X                  Drift length in cm (default: 100.0)
  --sigmat=X                 Transverse diffusion in mm (default: -1, compute from ion)
  --sigmal=X                 Longitudinal diffusion in mm (default: 0.0)
  --voxel-scale=X            Voxel scale factor (default: 1.0)
  --voxel-distance-scale=X   Voxel distance scale (default: 1.5)
  --mcvox=X                  MC path voxel size in mm (default: 1.0)
  --kde-bandwidth=X          RECO KDE bandwidth in mm (default: 5.0)
  --mc-kde-bandwidth=X       MC KDE bandwidth in mm (default: 1.0)
  --n-kde-eval=N             Number of KDE evaluation points (default: 200)
```

---

## Key Data Structures

### `RecoResult`

```julia
struct RecoResult
    event_id::Int
    track::Petit.Tracks
    path::DataFrame              # Raw path (x, y, z, s) - NOT smoothed
    track_length::Float64
    confidence::Float64
    mc_path::DataFrame           # MC path (x, y, z, energy, s)
    reco_s::Vector{Float64}      # Arc-length of each reco voxel
    kde_s::Vector{Float64}       # KDE evaluation points
    reco_kde_f::Vector{Float64}  # RECO energy density
    mc_kde_f::Vector{Float64}    # MC energy density
    kde_bandwidth::Float64       # RECO KDE bandwidth
    mc_kde_bandwidth::Float64    # MC KDE bandwidth
end
```

### Path DataFrame

```julia
# Columns: x, y, z, s
# Raw path from track walk (no smoothing applied)
# s = cumulative arc-length from start
```

### KDE Output

```julia
# kde_s: evaluation points from 0 to track_length
# reco_kde_f: f(s) = Σᵢ Eᵢ × K((s - sᵢ)/h_reco)  for RECO voxels (h_reco = kde_bandwidth)
# mc_kde_f: f(s) = Σᵢ Eᵢ × K((s - sᵢ)/h_mc)  for MC path (h_mc = mc_kde_bandwidth)
```

**Note:** Separate bandwidths allow less smoothing for MC (default 1 mm) than RECO (default 5 mm).

---

## Connection to Downstream Analysis

The HDF5 output from `track_reco_mt.jl` is read by:

```julia
results, metadata = Petit.read_reco_results_from_hdf5(reco_file)
```

This provides pre-computed KDE profiles for:

- `test_kde_track.jl` - Bragg peak feature analysis (uses `reco_kde_f`, `mc_kde_f`)
- `batch_blob_analysis.jl` - Blob energy analysis
- RECO vs MC comparison studies

---

## Key Functions

| Function | Description |
|----------|-------------|
| `process_single_event()` | Full reconstruction + KDE pipeline for one event |
| `compute_mc_path()` | Voxelizes MC hits, computes arc-length |
| `get_raw_path()` | Extracts unsmoothed path from track walk |
| `project_voxels_to_path()` | Maps voxels to arc-length coordinates |
| `energy_weighted_kde()` | Computes KDE with energy weights |
| `analysis_loop_reco_mt()` | Thread worker function |
| `event_loop_reco_mt()` | Main entry point, coordinates threads |
| `save_reco_results_to_hdf5()` | Writes results to HDF5 |

---

## Design Notes

### Why No Path Smoothing?

Previous versions used `filter_radius` to smooth the central path before projection. This created double smoothing:

1. Spatial smoothing of path (filter_radius)
2. KDE smoothing of energy density (bandwidth)

The new design projects voxels onto the **raw path** and applies only KDE smoothing. This:

- Avoids redundant smoothing
- Gives cleaner control via bandwidth parameters
- Produces sharper, more accurate energy profiles

### Separate MC Bandwidth

MC and RECO use different KDE bandwidths:

- **RECO** (`kde_bandwidth`, default 5 mm): larger smoothing to handle reconstruction noise
- **MC** (`mc_kde_bandwidth`, default 1 mm): minimal smoothing since MC is truth data

This allows fair comparison while preserving MC resolution.

---

## Related Files

- `scripts/test_kde_track.jl` - KDE feature analysis (uses pre-computed KDE)
- `scripts/batch_blob_analysis.jl` - Blob energy computation
- `src/Petit.jl` - Core reconstruction functions
- `src/kde.jl` - KDE implementation
- `src/tracks.jl` - Track data structures
- `src/track_io.jl` - HDF5 I/O functions
- `docs/track_reco_mt_refactoring_plan.md` - Implementation plan
