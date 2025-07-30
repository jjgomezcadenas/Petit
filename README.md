# Petit.jl

A Julia package for analyzing voxelized particle physics data, particularly for PET (Positron Emission Tomography) applications.

## Features

- **Voxel Analysis**: Convert hit data into voxel representations with energy aggregation
- **Graph-based Clustering**: Build connected graphs of voxels based on spatial proximity and energy thresholds
- **Statistical Analysis**: Histogram energy distributions, distances, and nearest neighbor analysis
- **Visualization**: Integrated plotting capabilities for 2D/3D data visualization

## Key Functions

### Voxelization
- `voxelize_hits()`: Convert particle hits into voxel representation

### Graph Construction
- `build_vgraph()`: Build connected graphs of voxels using optimized algorithms from Graphs.jl

### Analysis Tools
- `histogram_voxel_energy()`: Histogram energy distribution of voxels
- `histogram_voxel_distances()`: Analyze all pairwise distances between voxels
- `histogram_closest_distance()`: Distribution of nearest neighbor distances

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/jjgomezcadenas/Petit.git")
```

## Usage

```julia
using Petit

# Voxelize hit data
voxel_data = voxelize_hits(grouped_hits, 1.0)  # 1mm voxel size

# Build connected voxel graphs
graphs = build_vgraph(grouped_voxels, event_id; max_distance=1.5, energy_threshold=0.1)

# Analyze energy distribution
energy_hist = histogram_voxel_energy(grouped_voxels)
```

## Dependencies

- DataFrames.jl
- Graphs.jl  
- StatsBase.jl
- Plots.jl
- HDF5.jl (for data I/O)

## License

MIT License
