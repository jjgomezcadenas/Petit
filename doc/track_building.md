# Track Building Module

Functions for constructing `Tracks` objects from raw voxel data.

## Dependencies

```julia
using DataFrames
using Graphs
using NearestNeighbors
```

**Internal dependencies:** `get_event` from `event_and_hits.jl`

---

## Types

### `DiffusionParams`

Diffusion and voxelization parameters for a track.

```julia
struct DiffusionParams
    ldrift::Float64        # Drift length (mm)
    sigma_t::Float64       # Transverse diffusion (mm)
    sigma_l::Float64       # Longitudinal diffusion (mm)
    voxel_size::Float64    # Voxel size (mm)
    max_distance::Float64  # Max distance for clustering (mm)
    energy_threshold::Float64  # Energy threshold (keV)
    nbins_df::Int          # Number of bins for diffusion histogram
    nsigma_df::Float64     # Number of sigmas for histogram padding
end
```

**Default constructor:** `DiffusionParams()` creates parameters for MC tracks with no diffusion.

### `Tracks`

A track structure representing connected voxels in a particle track.

```julia
struct Tracks
    voxels::DataFrame              # Voxel data (x, y, z, energy, electrons)
    graph::SimpleGraph{Int}        # Graph connecting nearby voxels
    components::Vector{Vector{Int}} # Connected component indices
    diffusion::DiffusionParams     # Diffusion and voxelization parameters
end
```

---

## Functions

### `build_tracks`

```julia
build_tracks(hitsdf::DataFrame, event_id::Int;
             max_distance::Float64=1.5,
             energy_threshold::Float64=0.0,
             diffusion::DiffusionParams=DiffusionParams(),
             method::String="KDT",
             k::Int=10) -> Vector{Tracks}
```

Build tracks from hits DataFrame for a specific event.

**Arguments:**
- `hitsdf`: Hits DataFrame with event_id column
- `event_id`: Event ID to process
- `max_distance`: Max distance (mm) to connect voxels (default: 1.5)
- `energy_threshold`: Min energy (MeV) to include voxel (default: 0.0)
- `diffusion`: Diffusion parameters for the tracks
- `method`: Graph construction method:
  - `"KDT"` (default): Radius graph via KDTree
  - `"kNN"`: k-Nearest Neighbor graph
  - `"kNN_mutual"`: Mutual kNN (edge only if both vertices are in each other's kNN)
- `k`: Number of neighbors for kNN methods (default: 10)

**Returns:** Vector of Tracks (connected components)

---

### `build_tracks_kdtree`

```julia
build_tracks_kdtree(event_data::DataFrame;
                    max_distance::Float64=1.5,
                    energy_threshold::Float64=0.0,
                    diffusion::DiffusionParams=DiffusionParams()) -> Vector{Tracks}
```

Build tracks using radius graph via KDTree. Connects all voxel pairs within `max_distance`.

---

### `build_tracks_knn`

```julia
build_tracks_knn(event_data::DataFrame;
                 k::Int=10,
                 max_distance::Float64=Inf,
                 energy_threshold::Float64=0.0,
                 diffusion::DiffusionParams=DiffusionParams(),
                 mutual::Bool=false) -> Vector{Tracks}
```

Build tracks using k-Nearest Neighbor graph.

**Why use kNN?** kNN graphs prevent shortcut edges across U-shaped bends by limiting each voxel to its k closest neighbors, rather than all neighbors within a radius. This fixes extreme-finding failures on curved tracks.

**Arguments:**
- `k`: Number of nearest neighbors per voxel (default: 10)
- `max_distance`: Maximum edge length allowed (default: Inf)
- `mutual`: If true, keep only mutual kNN edges (default: false)

---

### `make_tracks`

```julia
make_tracks(hitsdf::DataFrame, event_id::Int;
            max_distance_mm::Float64=10.0,
            energy_threshold_kev::Float64=0.0,
            diffusion::DiffusionParams=DiffusionParams(),
            method::String="KDT",
            k::Int=10) -> Vector{Tracks}

make_tracks(event_data::DataFrame;
            max_distance_mm::Float64=1.0,
            energy_threshold_kev::Float64=0.0,
            diffusion::DiffusionParams=DiffusionParams(),
            method::String="KDT",
            k::Int=10) -> Vector{Tracks}
```

Convenience wrapper that:
1. Accepts parameters in mm and keV units
2. Sorts resulting tracks by total energy (descending)

---

## Usage Example

```julia
using Petit

# Load hits data
hitsdf = load_hits("data.h5")

# Build tracks for event 42
tracks = make_tracks(hitsdf, 42; max_distance_mm=2.0, method="kNN", k=8)

# Access the highest-energy track
main_track = tracks[1]
println("Voxels: ", nrow(main_track.voxels))
println("Total energy: ", sum(main_track.voxels.energy), " MeV")
```
