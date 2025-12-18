# Track Extreme Finding Module

Functions for finding track endpoints (extremes) using topology, curvature, energy, and spatial analysis.

## Dependencies

```julia
using LinearAlgebra
using Statistics
using DataFrames
using Graphs
using SparseArrays
```

**Internal dependencies:**
- `Tracks` from `track_building.jl`
- `euclidean_distance` from `voxels.jl`

---

## Types

### `TrackCoords`

Pre-extracted coordinates for faster access (avoids DataFrame column overhead).

```julia
struct TrackCoords
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end
```

---

## Main API Functions

### `find_track_extremes`

```julia
find_track_extremes(track::Tracks;
                    use_energy_weighting::Bool=true,
                    use_edge_energy_weighting::Bool=true,
                    use_mst_fallback::Bool=false,
                    dense_track_threshold::Float64=6.0)
    -> (extreme1_idx, extreme2_idx, path, confidence)
```

Find track endpoints using combined topology/curvature/energy analysis.

**Arguments:**
- `track`: Track object
- `use_energy_weighting`: Enable energy-weighted extreme finding (default: true)
- `use_edge_energy_weighting`: Enable edge-energy-weighted Dijkstra method (default: true)
- `use_mst_fallback`: Enable MST-based fallback method (default: false)
- `dense_track_threshold`: Average degree threshold for dense tracks (default: 6.0)

**Returns:** Tuple of (extreme1_idx, extreme2_idx, path, confidence)
- `extreme1_idx`, `extreme2_idx`: Vertex indices of endpoints
- `path`: Vector of vertex indices connecting the endpoints
- `confidence`: Score from 0 to 1 indicating reliability

---

### `walk_track_from_extremes`

```julia
walk_track_from_extremes(track::Tracks) -> NamedTuple
```

Walk through a track from one endpoint to the other.

**Returns:** NamedTuple with:
- `extremes`: (start_voxel, end_voxel) DataFrameRows
- `path_indices`: Vertex indices along the path
- `path_voxels`: DataFrame of voxels in path order
- `total_length`: Path length in mm
- `confidence`: Confidence score (0-1)

---

### `find_track_extremes` (convenience)

```julia
find_track_extremes(trk::Vector{Tracks}; i=1)
```

Convenience function: find extremes for i-th track in a vector.

---

## Algorithm Functions

### `find_extremes_topology`

```julia
find_extremes_topology(track::Tracks, coords::TrackCoords)
    -> (extreme1, extreme2, path, confidence, path_length)
```

Topology-based extreme finding. Looks for degree-1 vertices (natural endpoints).
Includes path coverage check to avoid false positives from adjacent degree-1 vertices.

---

### `find_extremes_curvature`

```julia
find_extremes_curvature(track::Tracks, coords::TrackCoords)
    -> (extreme1, extreme2, path, confidence, path_length)
```

Curvature-based extreme finding. Track endpoints typically have low curvature (straight segments).

---

### `find_extremes_spatial`

```julia
find_extremes_spatial(track::Tracks, coords::TrackCoords)
    -> (extreme1, extreme2, path, confidence, path_length)
```

Spatial-based extreme finding for dense tracks. Uses spatial extremes (min/max coordinates) and minimum degree vertices as candidates.

---

### `find_extremes_energy_weighted`

```julia
find_extremes_energy_weighted(track::Tracks, coords::TrackCoords;
                              min_coverage::Float64=0.6)
    -> (extreme1, extreme2, path, confidence, path_length)
```

Energy-weighted extreme finding. Designed for Bragg peak detection where track endpoints have high energy deposits.

**Confidence calibration based on:**
1. Energy at endpoints (Bragg peak detection)
2. Path coverage (path_length / track_extent)
3. Endpoint separation quality

---

### `find_extremes_edge_energy_weighted_opt`

```julia
find_extremes_edge_energy_weighted_opt(track::Tracks, coords::TrackCoords;
                                       epsE::Float64=1e-6,
                                       α::Float64=1.0,
                                       β::Float64=1.0)
    -> (extreme1, extreme2, path, confidence, path_length)
```

Energy-weighted edge-cost extreme finding using Dijkstra's algorithm.

**Algorithm:**
1. Build sparse weight matrix: `w_ij = (d_ij^α) / ((0.5*(E_i + E_j) + epsE)^β)`
2. Double-sweep using weighted Dijkstra distances
3. Return endpoints and the weighted-shortest path

This suppresses "geometric shortcuts" through low-energy regions.

---

### `find_extremes_mst_diameter`

```julia
find_extremes_mst_diameter(track::Tracks, coords::TrackCoords)
    -> (extreme1, extreme2, path, confidence, path_length)
```

MST-based fallback method:
1. Build MST from track graph
2. Find diameter endpoints by double BFS
3. Return BFS path and Euclidean path length

---

### `find_extremes_combined`

```julia
find_extremes_combined(track::Tracks, coords::TrackCoords;
                       use_energy_weighting::Bool=true,
                       use_edge_energy_weighting::Bool=true,
                       use_mst_fallback::Bool=false,
                       dense_track_threshold::Float64=6.0)
    -> (extreme1, extreme2, path, confidence)
```

Main dispatcher that combines all methods with early-exit strategies.

**Strategy:**
1. Always try topology first (fast, often sufficient)
2. For dense tracks (avg_degree > threshold): prefer energy-weighted methods
3. For sparse tracks: compare topology vs curvature
4. Pick best result by path length and confidence

---

## Utility Functions

### `extract_coords`

```julia
extract_coords(track::Tracks) -> TrackCoords
```

Extract coordinates from track voxels once for efficient reuse.

### `find_path_bfs`

```julia
find_path_bfs(g::SimpleGraph, start_vertex::Int, end_vertex::Int) -> Vector{Int}
```

Optimized BFS using parent pointers. Returns path as vector of vertex indices.

### `calculate_path_length`

```julia
calculate_path_length(coords::TrackCoords, path::Vector{Int}) -> Float64
```

Calculate Euclidean path length using pre-extracted coordinates.

### `calculate_path_length_from_coords`

```julia
calculate_path_length_from_coords(coords::AbstractMatrix{<:Real},
                                  path::AbstractVector{Int}) -> Float64
```

Calculate path length from a 3xN coordinate matrix.

### `calculate_vertex_curvatures`

```julia
calculate_vertex_curvatures(track::Tracks, coords::TrackCoords) -> Vector{Float64}
```

Compute curvature at each vertex. Low curvature indicates straight segments (potential endpoints).

---

## Usage Example

```julia
using Petit

# Build tracks
tracks = make_tracks(event_data; max_distance_mm=2.0)
track = tracks[1]

# Find extremes
extreme1, extreme2, path, confidence = find_track_extremes(track)
println("Confidence: ", confidence)

# Get detailed walk result
walk = walk_track_from_extremes(track)
println("Path length: ", walk.total_length, " mm")
println("Start voxel energy: ", walk.extremes[1].energy)
println("End voxel energy: ", walk.extremes[2].energy)
```
