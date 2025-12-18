# Track Modules Overview

The track analysis functionality in Petit is split into four modules, each handling a specific aspect of track reconstruction and analysis.

## Module Hierarchy

```
voxels.jl                    # euclidean_distance
    │
    └── track_building.jl    # Tracks, DiffusionParams, build_tracks
            │
            └── track_extreme_finding.jl  # find_track_extremes, walk_track_from_extremes
                    │
                    ├── track_analysis.jl      # diagnostics, utilities
                    │
                    └── track_blob_analysis.jl # blob energy analysis
```

## Required Inclusion Order

When including these modules in `Petit.jl`, they must be loaded in dependency order:

```julia
include("voxels.jl")
include("track_building.jl")
include("track_extreme_finding.jl")
include("track_analysis.jl")
include("track_blob_analysis.jl")
```

---

## Module Summaries

### [track_building.jl](track_building.md)

**Purpose:** Construct `Tracks` objects from raw voxel data.

**Key exports:**
- `Tracks` - Main track data structure
- `DiffusionParams` - Track parameters
- `build_tracks` - Build tracks from hits
- `make_tracks` - Convenience wrapper

**Methods:**
- KDTree radius graph (default)
- k-Nearest Neighbor graph
- Mutual kNN graph

---

### [track_extreme_finding.jl](track_extreme_finding.md)

**Purpose:** Find track endpoints using multiple algorithms.

**Key exports:**
- `find_track_extremes` - Main entry point
- `walk_track_from_extremes` - Get path between endpoints
- `TrackCoords` - Efficient coordinate storage

**Algorithms:**
- Topology-based (degree-1 vertices)
- Curvature-based (low curvature = endpoint)
- Spatial-based (coordinate extremes)
- Energy-weighted (Bragg peak detection)
- Edge-weighted Dijkstra
- MST diameter fallback

---

### [track_analysis.jl](track_analysis.md)

**Purpose:** Diagnostic utilities for track quality assessment.

**Key exports:**
- `track_positions` - Extract all coordinates
- `track_energies_keV` - Get track energies
- `diagnose_path_efficiency` - Path tortuosity
- `diagnose_skeleton_coverage` - Energy coverage
- `diagnose_endpoint_stability` - Algorithm comparison

---

### [track_blob_analysis.jl](track_blob_analysis.md)

**Purpose:** Analyze energy blobs at track endpoints for signal/background discrimination.

**Key exports:**
- `find_blob_energies` - Blob analysis with smoothed path
- `energy_in_spheres_around_extremes` - Fixed radius blobs
- `energy_in_variable_spheres_around_extremes` - Adaptive radius
- `blob_asymmetry` - Batch asymmetry analysis

**Key metric:** `asymmetry = (Eb1 - Eb2) / (Eb1 + Eb2)`

---

## Complete Analysis Pipeline

```julia
using Petit

# 1. Load and build tracks
hitsdf = load_hits("data.h5")
tracks = make_tracks(hitsdf, event_id; max_distance_mm=2.0, method="kNN")

# 2. Select main track (highest energy)
track = tracks[1]

# 3. Find extremes
walk = walk_track_from_extremes(track)
println("Track length: ", walk.total_length, " mm")
println("Confidence: ", walk.confidence)

# 4. Analyze blobs
blobs = energy_in_spheres_around_extremes(track, walk, 10.0)
Eb1 = blobs.blob1_energy * 1e3  # keV
Eb2 = blobs.blob2_energy * 1e3  # keV
asymmetry = (Eb1 - Eb2) / (Eb1 + Eb2)

println("Eb1: ", Eb1, " keV")
println("Eb2: ", Eb2, " keV")
println("Asymmetry: ", asymmetry)

# 5. Diagnose quality
coords = extract_coords(track)
η, L, D = diagnose_path_efficiency(coords, walk.path_indices)
println("Path efficiency η: ", η)

f, _, _ = diagnose_skeleton_coverage(track, coords, walk.path_indices; R_cover=3.0)
println("Energy coverage: ", round(f*100, digits=1), "%")
```

---

## Selection Criteria

Typical cuts for double-beta decay analysis:

| Cut | Purpose |
|-----|---------|
| `confidence > 0.7` | Reliable endpoint finding |
| `asymmetry < 0.4` | Two Bragg peaks (signal-like) |
| `total_length > 100 mm` | Minimum track length |
| `η < 2.0` | Not too tortuous |
| `coverage > 0.8` | Path captures most energy |
