# Track Analysis Module

Utility functions and diagnostics for track analysis.

## Dependencies

```julia
using DataFrames
using Graphs
using NearestNeighbors
```

**Internal dependencies:**
- `Tracks` from `track_building.jl`
- `TrackCoords`, `calculate_path_length_from_coords` from `track_extreme_finding.jl`
- `euclidean_distance` from `voxels.jl`

---

## Data Extraction Functions

### `track_positions`

```julia
track_positions(tracks::Vector{Tracks}) -> (X, Y, Z)
```

Extract all voxel positions from a vector of tracks.

**Returns:** Tuple of three vectors (X, Y, Z) containing all coordinates.

---

### `track_energies_keV`

```julia
track_energies_keV(tracks::Vector{Tracks}) -> Vector{Float64}
```

Compute total energy (keV) for each track in the collection.

**Returns:** Vector of energies in keV, one per track.

---

## Diagnostic Functions

### `diagnose_path_efficiency`

```julia
diagnose_path_efficiency(coords::TrackCoords, path::Vector{Int})
    -> (η, L_path, D_end)
```

Compute path efficiency ratio.

**Returns:**
- `η = L_path / D_end`: Ratio of path length to endpoint distance
- `L_path`: Euclidean length along the vertex sequence
- `D_end`: Straight-line distance between endpoints

**Interpretation:**
- `η ≈ 1`: Nearly straight path
- `η > 1`: Curved path (larger values = more tortuous)
- `η >> 1`: Highly curved or coiled track

---

### `diagnose_endpoint_degrees`

```julia
diagnose_endpoint_degrees(g::SimpleGraph, u::Int, v::Int) -> (deg_u, deg_v)
```

Get the graph degree of two endpoint vertices.

**Use case:** Degree-1 vertices are natural track endpoints. Higher degrees suggest branching or dense regions.

---

### `diagnose_skeleton_coverage`

```julia
diagnose_skeleton_coverage(track::Tracks, coords::TrackCoords, path::Vector{Int};
                           R_cover::Float64,
                           energy_col::Symbol=:energy)
    -> (f, E_in, E_tot)
```

Compute fraction of total energy within distance `R_cover` of the skeleton path.

**Algorithm:** Uses KDTree over skeleton vertices and nearest-neighbor queries to find voxels within `R_cover` of any path vertex.

**Returns:**
- `f`: Fraction of energy covered (E_in / E_tot)
- `E_in`: Energy within coverage radius
- `E_tot`: Total track energy

**Use case:** Validates that the reconstructed path captures most of the track's energy.

---

### `diagnose_endpoint_stability`

```julia
diagnose_endpoint_stability(res1, res2, coords::TrackCoords) -> Δ
```

Compute endpoint displacement between two solutions.

**Arguments:**
- `res1`, `res2`: Tuples with endpoints as first two elements (u, v, ...)
- `coords`: Track coordinates

**Returns:** Maximum displacement (mm) between corresponding endpoints, allowing for endpoint swapping to minimize total displacement.

**Use case:** Compare results from different extreme-finding algorithms to assess stability.

---

### `diagnose_track_extent`

```julia
diagnose_track_extent(coords::TrackCoords,
                      vertices::AbstractVector{Int}=1:length(coords.x))
    -> extent_mm
```

Compute spatial extent (diagonal of bounding box) for selected vertices.

**Returns:** Extent in mm (sqrt of sum of squared ranges in x, y, z).

**Use case:** Condition shortcut detection - η near 1 is only suspicious if extent is large.

---

## Usage Examples

### Basic Track Analysis

```julia
using Petit

tracks = make_tracks(event_data; max_distance_mm=2.0)

# Get all positions for plotting
X, Y, Z = track_positions(tracks)

# Get energies
energies = track_energies_keV(tracks)
println("Track energies: ", energies, " keV")
```

### Path Quality Diagnostics

```julia
track = tracks[1]
coords = extract_coords(track)
extreme1, extreme2, path, confidence = find_track_extremes(track)

# Check path efficiency
η, L_path, D_end = diagnose_path_efficiency(coords, path)
println("Path efficiency η = ", η)
println("Path length: ", L_path, " mm")
println("Endpoint distance: ", D_end, " mm")

# Check endpoint degrees
deg1, deg2 = diagnose_endpoint_degrees(track.graph, extreme1, extreme2)
println("Endpoint degrees: ", deg1, ", ", deg2)

# Check skeleton coverage
f, E_in, E_tot = diagnose_skeleton_coverage(track, coords, path; R_cover=3.0)
println("Energy coverage at R=3mm: ", round(f*100, digits=1), "%")
```

### Comparing Algorithms

```julia
# Get results from two methods
result1 = find_extremes_topology(track, coords)
result2 = find_extremes_energy_weighted(track, coords)

# Check stability
displacement = diagnose_endpoint_stability(result1, result2, coords)
println("Endpoint displacement: ", displacement, " mm")

if displacement < 5.0
    println("Methods agree well")
else
    println("Methods disagree - check track quality")
end
```
