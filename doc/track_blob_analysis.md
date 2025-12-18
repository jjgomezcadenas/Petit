# Track Blob Analysis Module

Functions for analyzing energy blobs at track endpoints.

## Dependencies

```julia
using DataFrames
using NearestNeighbors
using Statistics
```

**Internal dependencies:**
- `Tracks` from `track_building.jl`
- `walk_track_from_extremes` from `track_extreme_finding.jl`
- `euclidean_distance` from `voxels.jl`

---

## Overview

In double-beta decay experiments, the two electrons deposit energy along their tracks with characteristic "Bragg peaks" at the endpoints. Blob analysis quantifies the energy in spherical regions around track endpoints to distinguish signal (two high-energy blobs) from background (typically one high-energy blob).

**Key metrics:**
- `Eb1`: Energy of the higher-energy blob
- `Eb2`: Energy of the lower-energy blob
- `asymmetry = (Eb1 - Eb2) / (Eb1 + Eb2)`: Signal events have low asymmetry

---

## Main Functions

### `find_blob_energies`

```julia
find_blob_energies(track::Tracks, central_path::DataFrame;
                   radius::Float64) -> NamedTuple
```

Sum energy in spheres around endpoints of a reconstructed central path.

**Arguments:**
- `track`: Must have `track.voxels` DataFrame with x, y, z, energy columns
- `central_path`: Smoothed path DataFrame with x, y, z columns
- `radius`: Sphere radius in mm

**Returns:** NamedTuple with:
- `blob1`: (x, y, z, energy, voxels_idx, endpoint_index) - higher energy blob
- `blob2`: (x, y, z, energy, voxels_idx, endpoint_index) - lower energy blob
- `Eb1`: Energy of blob1 in keV
- `Eb2`: Energy of blob2 in keV
- `asymmetry`: (Eb1 - Eb2) / (Eb1 + Eb2)

---

### `energy_in_spheres_around_extremes`

```julia
energy_in_spheres_around_extremes(track::Tracks, walk_result, radius::Float64)
    -> NamedTuple

energy_in_spheres_around_extremes(track::Tracks, radius::Float64)
    -> NamedTuple  # convenience version
```

Calculate energy within spheres of given radius around track endpoints from `walk_track_from_extremes`.

**Returns:** NamedTuple with:
- `start_sphere_energy`, `end_sphere_energy`: Energy at each endpoint
- `start_voxel_count`, `end_voxel_count`: Voxel counts
- `start_center`, `end_center`: Coordinates (x, y, z)
- `blob1_energy`, `blob2_energy`: Sorted by energy (blob1 >= blob2)
- `blob1_voxel_count`, `blob2_voxel_count`
- `blob1_center`, `blob2_center`

---

### `energy_blobs_from_path`

```julia
energy_blobs_from_path(walk_result, n::Int;
                       energy_col::Symbol=:energy) -> NamedTuple
```

Calculate blob energies using n voxels from each end of the walk path.

**Arguments:**
- `walk_result`: Result from `walk_track_from_extremes`
- `n`: Number of voxels to include from each end
- `energy_col`: Column to use (:energy or :electrons)

**Returns:** NamedTuple with:
- `eb1`, `eb2`: Blob energies (eb1 >= eb2)
- `n1`, `n2`: Number of voxels in each blob
- `start_energy`, `end_energy`: Energy at start/end of path

**Note:** This method uses path ordering rather than spatial spheres, which can be more robust for curved tracks.

---

### `energy_in_variable_spheres_around_extremes`

```julia
energy_in_variable_spheres_around_extremes(track::Tracks, walk_result;
                                           seed_radius::Float64=3.0,
                                           step::Float64=1.0,
                                           max_radius::Float64=10.0,
                                           threshold::Float64=0.05) -> NamedTuple
```

Adaptive radius expansion starting from `seed_radius`.

**Algorithm:**
1. Start with seed_radius
2. Calculate energy within sphere
3. Expand radius by step
4. If relative change < threshold, stop
5. If radius >= max_radius, stop
6. Repeat

**Returns:** NamedTuple with:
- `blob1_energy`, `blob2_energy`
- `blob1_radius`, `blob2_radius`: Final radii used
- `blob1_voxel_count`, `blob2_voxel_count`
- `blob1_center`, `blob2_center`
- `blob1_history`, `blob2_history`: Vector of (radius, energy) pairs

---

## Batch Analysis Functions

### `blob_asymmetry`

```julia
blob_asymmetry(tracks::Vector{Tracks};
               i0::Int=1, il::Int=10,
               r0::Int=5, rl::Int=15) -> (mean_asymmetry, std_asymmetry)
```

Compute blob asymmetry over a range of events and sphere radii.

**Arguments:**
- `tracks`: Vector of track objects
- `i0`, `il`: First and last event indices
- `r0`, `rl`: Min and max sphere radii in mm

**Returns:** Tuple of vectors (mean, std) for each radius value.

---

### `blob_asymmetry_from_path`

```julia
blob_asymmetry_from_path(tracks::Vector{Tracks};
                         i0::Int=1, il::Int=10,
                         n0::Int=1, nl::Int=10,
                         energy_col::Symbol=:energy) -> (mean_asymmetry, std_asymmetry)
```

Compute blob asymmetry using path voxels directly (n voxels from each end).

**Arguments:**
- `n0`, `nl`: Min and max number of voxels from each end

---

### `blob_analysis_vs_radius`

```julia
blob_analysis_vs_radius(track, walk_result;
                        r0::Int=5, rl::Int=15) -> NamedTuple

blob_analysis_vs_radius(track; r0::Int=5, rl::Int=15) -> NamedTuple
```

Compute blob energies, asymmetry, and track length for a single track across radii.

**Returns:** NamedTuple with:
- `EB1`, `EB2`: Vectors of blob energies (MeV) for each radius
- `DB`: Vector of asymmetries
- `DL`: Track length (mm)
- `radii`: Vector of radius values used

---

## Usage Examples

### Basic Blob Analysis

```julia
using Petit

tracks = make_tracks(event_data; max_distance_mm=2.0)
track = tracks[1]

# Simple blob analysis at fixed radius
walk = walk_track_from_extremes(track)
blobs = energy_in_spheres_around_extremes(track, walk, 10.0)

println("Blob1 energy: ", blobs.blob1_energy * 1e3, " keV")
println("Blob2 energy: ", blobs.blob2_energy * 1e3, " keV")

asymmetry = abs(blobs.blob1_energy - blobs.blob2_energy) /
            (blobs.blob1_energy + blobs.blob2_energy)
println("Asymmetry: ", asymmetry)
```

### Radius Scan

```julia
# Scan across radii
results = blob_analysis_vs_radius(track; r0=5, rl=20)

# Find optimal radius (minimum asymmetry)
min_idx = argmin(results.DB)
optimal_radius = results.radii[min_idx]
println("Optimal radius: ", optimal_radius, " mm")
println("Asymmetry at optimal: ", results.DB[min_idx])
```

### Adaptive Blob Finding

```julia
# Let the algorithm find optimal blob sizes
adaptive = energy_in_variable_spheres_around_extremes(track, walk;
                                                       seed_radius=3.0,
                                                       max_radius=15.0)

println("Blob1: ", adaptive.blob1_energy * 1e3, " keV at R=", adaptive.blob1_radius, " mm")
println("Blob2: ", adaptive.blob2_energy * 1e3, " keV at R=", adaptive.blob2_radius, " mm")
```

### Batch Analysis

```julia
# Analyze many events
tracks = [make_tracks(get_event(hitsdf, i))[1] for i in 1:100]

mean_asym, std_asym = blob_asymmetry(tracks; i0=1, il=100, r0=8, rl=12)

# Plot asymmetry vs radius
using Plots
plot(8:12, mean_asym, ribbon=std_asym,
     xlabel="Radius (mm)", ylabel="Asymmetry",
     label="Mean +/- std")
```

### Signal vs Background Discrimination

```julia
# For signal (ββ0ν), expect two Bragg peaks → low asymmetry
# For background (single electron), expect one Bragg peak → high asymmetry

function classify_event(track; radius=10.0, threshold=0.4)
    blobs = energy_in_spheres_around_extremes(track, radius)
    eb1, eb2 = blobs.blob1_energy, blobs.blob2_energy
    asymmetry = (eb1 - eb2) / (eb1 + eb2)

    return asymmetry < threshold ? :signal : :background
end
```
