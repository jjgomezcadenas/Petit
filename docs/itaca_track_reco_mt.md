# ITACA Track Reconstruction - Multi-threaded

## Overview

`itaca_track_reco_mt.jl` is a multi-threaded batch processing script that combines:

- **Analysis pipeline** from `itaca_single_track_analysis.jl` (diffusion, voxelization, track finding, blob analysis, KDE peaks)
- **Multi-threading and HDF5 persistence** from `track_reco_mt.jl`

It processes Monte Carlo events in parallel, selects single-track events, performs full ITACA analysis, and saves:

1. **CSV summary file** with analysis metrics for all events
2. **HDF5 files** with full track data for events passing cuts

No plotting is performed (batch mode). Default print level is `quiet`.

## Location

```
scripts/itaca_track_reco_mt.jl
```

## Architecture: Source of Each Component

### From `itaca_single_track_analysis.jl` (ALL analysis logic)

#### Helper Functions (to copy/adapt)

| Function | Purpose |
|----------|---------|
| `get_sigma(particle_type, ldrft; ...)` | Compute σt, σl based on particle type (ion vs electron) |
| `get_voxel_size_and_distance(ldrft, σt)` | Dynamic voxel sizing based on drift length |
| `get_energy_threshold(particle_type; ...)` | Energy threshold conversion for ions |

#### Event Processing Pipeline (from `main()` loop)

| Step | Function | Description |
|------|----------|-------------|
| 1 | `Petit.get_event(hitsdf, nevent)` | Load event hits |
| 2 | `Petit.compute_mc_path(event_df, mcvox_size)` | Compute MC truth path |
| 3 | `Petit.transform_hits_df(event_df)` | Add electron counts |
| 4 | `Petit.diffuse_xyz_image_mc(event_mc; σt, σl, ...)` | Apply Gaussian diffusion |
| 5 | `Petit.voxelize_event(diffused_df, voxel_size)` | Create 3D voxels |
| 6 | `Petit.make_tracks(voxels; max_distance, eth, ...)` | Graph-based track finding |
| 7 | Filter: `length(tracks) == 1` | Keep only single-track events |
| 8 | `Petit.walk_track_from_extremes(track)` | Find path through track |
| 9 | `Petit.get_raw_path(track, path_indices)` | Extract path DataFrame |
| 10 | `Petit.compute_extreme_distances(path, mc_path)` | RECO vs MC distances |
| 11 | `Petit.get_reco_kde(track, path; bandwidth, n_eval)` | RECO energy density |
| 12 | `Petit.get_mc_kde(mc_path; bandwidth, n_eval)` | MC energy density |
| 13 | `Petit.find_peaks(kde_f, kde_s; prom_scale)` | Peak detection |
| 14 | `kde_peaks(peaks, kde_f)` | Extract peak1/peak2 info (from `itaca_aux.jl`) |
| 15 | `Petit.find_blob_energies(track, path; radius)` | Blob energy analysis |

#### CSV Output Structure

Same 14 columns as `analysis_results.csv` from `itaca_single_track_analysis.jl`.

### From `track_reco_mt.jl` (ONLY multi-threading + persistence)

#### Multi-threading Pattern

| Component | Description |
|-----------|-------------|
| `Petit.load_and_validate_input()` | Load HDF5 and validate event range |
| `Petit.get_optimal_threads()` | Determine optimal thread count |
| `Petit.split_events_for_threads()` | Divide events across threads |
| `Threads.@threads` loop | Parallel processing pattern |
| Thread-local results collection | Each thread accumulates its own results |

#### HDF5 Persistence (adapted)

| Function | Description |
|----------|-------------|
| `save_path_to_hdf5(path, group)` | Save path DataFrame to HDF5 group |
| `save_mc_path_to_hdf5(mc_path, group)` | Save MC path to HDF5 group |
| `save_itaca_results_to_hdf5(results, path, metadata)` | Main save function (extended for blob/peak data) |

## Algorithm Step-by-Step

```
1. INITIALIZATION
   │
   │  [From itaca_single_track_analysis.jl]
   ├── Parse CLI arguments
   ├── get_sigma() → compute σt, σl based on particle type
   ├── get_voxel_size_and_distance() → voxel_size, max_distance
   ├── get_energy_threshold() → eth
   └── Create Petit.DiffusionParams

2. SETUP MULTI-THREADING
   │
   │  [From track_reco_mt.jl]
   ├── Petit.load_and_validate_input()
   ├── Petit.get_optimal_threads()
   └── Petit.split_events_for_threads()

3. PARALLEL PROCESSING
   │
   │  [MT pattern from track_reco_mt.jl]
   │  [Analysis from itaca_single_track_analysis.jl]
   │
   Threads.@threads for each thread range:
   │
   │  FOR each event in range:
   │  │
   │  │  ┌─────────────────────────────────────────────────────┐
   │  │  │ All from itaca_single_track_analysis.jl main() loop │
   │  │  └─────────────────────────────────────────────────────┘
   │  │  ├── Petit.get_event()
   │  │  ├── Petit.compute_mc_path()
   │  │  ├── Petit.transform_hits_df()
   │  │  ├── Petit.diffuse_xyz_image_mc()
   │  │  ├── Petit.voxelize_event()
   │  │  ├── Petit.make_tracks()
   │  │  ├── if length(tracks) != 1 → skip event
   │  │  ├── Petit.walk_track_from_extremes()
   │  │  ├── Petit.get_raw_path()
   │  │  ├── Petit.compute_extreme_distances()
   │  │  ├── Petit.get_reco_kde()
   │  │  ├── Petit.get_mc_kde()
   │  │  ├── Petit.find_peaks()
   │  │  ├── kde_peaks()
   │  │  ├── Petit.find_blob_energies()
   │  │  └── Store result in ItacaResult struct
   │  │
   │  END FOR
   │
   │  ┌──────────────────────────────┐
   │  │ From track_reco_mt.jl       │
   │  └──────────────────────────────┘
   │  └── save_itaca_results_to_hdf5() → write thread's HDF5 file
   │
   END THREADS

4. MERGE RESULTS
   ├── Collect CSV rows from all threads
   ├── Write combined <output_base>_analysis.csv
   └── Write <output_base>_metadata.csv

5. SUMMARY
   └── Print totals (events processed, single-track count)
```

## Data Structures

### ItacaResult Struct

Combines `RecoResult` from `track_reco_mt.jl` with blob/peak data from `itaca_single_track_analysis.jl`:

```julia
struct ItacaResult
    # ─────────────────────────────────────────────────────────
    # From track_reco_mt.jl RecoResult
    # ─────────────────────────────────────────────────────────
    event_id::Int
    track::Petit.Tracks
    path::DataFrame              # Raw path (x, y, z, s)
    track_length::Float64
    confidence::Float64
    mc_path::DataFrame           # MC path (x, y, z, energy, s)
    reco_s::Vector{Float64}      # Arc-length of voxels on path
    kde_s::Vector{Float64}       # RECO KDE evaluation points
    reco_kde_f::Vector{Float64}  # RECO energy density f(s)
    mc_kde_s::Vector{Float64}    # MC KDE evaluation points
    mc_kde_f::Vector{Float64}    # MC energy density f(s)
    kde_bandwidth::Float64
    mc_kde_bandwidth::Float64
    d1::Float64                  # RECO-MC extreme distance 1
    d2::Float64                  # RECO-MC extreme distance 2

    # ─────────────────────────────────────────────────────────
    # NEW from itaca_single_track_analysis.jl
    # ─────────────────────────────────────────────────────────
    Eb1::Float64                 # High-energy blob (keV)
    Eb2::Float64                 # Low-energy blob (keV)
    asymmetry::Float64           # (Eb1-Eb2)/(Eb1+Eb2)
    blob1_pos::NTuple{3,Float64} # Blob1 position (x, y, z)
    blob2_pos::NTuple{3,Float64} # Blob2 position (x, y, z)
    blob_radius::Float64         # Blob radius used (mm)
    n_peaks::Int                 # Number of KDE peaks
    peak1_left::Float64          # Peak1 left edge (mm)
    peak1_right::Float64         # Peak1 right edge (mm)
    peak1_prom::Float64          # Peak1 prominence (normalized)
    peak2_left::Float64          # Peak2 left edge (mm)
    peak2_right::Float64         # Peak2 right edge (mm)
    peak2_prom::Float64          # Peak2 prominence (normalized)
end
```

## Usage

```bash
julia -t <nthreads> scripts/itaca_track_reco_mt.jl <cmdir> <input_file> [options]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `cmdir` | Directory containing the input HDF5 file |
| `input_file` | Name of the input HDF5 file |

### Optional Arguments

#### Output Control

| Option | Description | Default |
|--------|-------------|---------|
| `--outdir=PATH` | Output directory (created if needed) | `results` |
| `--outbase=NAME` | Base name for output files | `itaca` |
| `--print_level=LEVEL` | `quiet`, `verbose`, `very_verbose` | `quiet` |

#### Event Selection

| Option | Description | Default |
|--------|-------------|---------|
| `--ievt=N` | First event to process | 1 |
| `--levt=N` | Last event to process (-1 = all) | -1 |
| `--nthreads=N` | Number of threads | 1 |

#### Physics Parameters

| Option | Description | Default |
|--------|-------------|---------|
| `--particle=TYPE` | Particle type: `ion` or `electron` | `ion` |
| `--ldrft=X` | Drift length (cm) | 100.0 |
| `--tK=X` | Temperature (K) | 297.0 |
| `--edrift=X` | Drift field (V/cm) | 500.0 |
| `--Pbar=X` | Pressure (bar) | 15.0 |
| `--dt=X` | Transverse diffusion coeff (mm/√cm) | 3.5 |
| `--dl=X` | Longitudinal diffusion coeff (mm/√cm) | 0.9 |

#### Analysis Parameters

| Option | Description | Default |
|--------|-------------|---------|
| `--eth-ion=X` | Energy threshold for ions | 10.0 |
| `--nbins=N` | Diffusion histogram bins | 100 |
| `--nsigma=X` | Diffusion extent (sigma) | 3.0 |
| `--nkde=N` | KDE evaluation points | 200 |
| `--Rb=X` | Blob radius (mm) | 10.0 |

## Output Files

All output files are created in `--outdir` with `--outbase` prefix.

### 1. CSV Summary File: `<outdir>/<outbase>_analysis.csv`

One row per single-track event with all analysis metrics.

| Column | Description | Units |
|--------|-------------|-------|
| `event` | Event number | - |
| `thread_id` | Thread that processed this event | - |
| `track_length_mm` | Reconstructed track length | mm |
| `confidence` | Path-finding confidence score | - |
| `Eb1_keV` | High-energy blob energy | keV |
| `Eb2_keV` | Low-energy blob energy | keV |
| `asymmetry` | Blob asymmetry: (Eb1-Eb2)/(Eb1+Eb2) | - |
| `d1_mm` | Distance: RECO extreme 1 to MC | mm |
| `d2_mm` | Distance: RECO extreme 2 to MC | mm |
| `n_peaks` | Number of KDE peaks found | - |
| `peak1_left` | Leftmost peak left edge | mm |
| `peak1_right` | Leftmost peak right edge | mm |
| `peak1_prom` | Leftmost peak prominence (normalized) | - |
| `peak2_left` | Rightmost peak left edge | mm |
| `peak2_right` | Rightmost peak right edge | mm |
| `peak2_prom` | Rightmost peak prominence (normalized) | - |

### 2. HDF5 Track Files: `<outdir>/<outbase>_th_<N>.h5`

One file per thread containing full track data for offline analysis.

#### File Attributes (metadata)

| Attribute | Description |
|-----------|-------------|
| `input_file` | Source HDF5 file name |
| `thread_id` | Thread identifier |
| `first_event_processed` | First event in this thread's range |
| `last_event_processed` | Last event in this thread's range |
| `events_processed` | Total events attempted by this thread |
| `total_tracks_saved` | Number of single-track events saved |
| `ldrift_cm` | Drift length |
| `sigma_t_mm` | Transverse diffusion |
| `sigma_l_mm` | Longitudinal diffusion |
| `voxel_size_mm` | Voxel size used |
| `max_distance_mm` | Max voxel connection distance |
| `energy_threshold_kev` | Energy threshold |
| `kde_bandwidth_mm` | RECO KDE bandwidth |
| `mc_kde_bandwidth_mm` | MC KDE bandwidth |
| `blob_radius_mm` | Blob analysis radius |
| `n_kde_eval` | Number of KDE evaluation points |

#### Per-Track Data: `batch_1/track_<idx>/`

Each track group contains:

##### Event Metadata

| Dataset | Type | Description |
|---------|------|-------------|
| `event_id` | Int | Original event number |
| `track_length` | Float64 | Track length (mm) |
| `confidence` | Float64 | Path-finding confidence |

##### Voxel Data

| Dataset | Type | Description |
|---------|------|-------------|
| `voxels` | Matrix | Voxel data (N × 4): x, y, z, energy |
| `voxel_columns` | String[] | Column names: ["x", "y", "z", "energy"] |
| `n_vertices` | Int | Number of vertices in track graph |
| `graph_edges` | Matrix | Edge list (E × 2): source, destination |
| `components` | Matrix | Connected component indices |

##### Path Data

| Dataset | Type | Description |
|---------|------|-------------|
| `path_data` | Matrix | Raw path (M × 4): x, y, z, s |
| `path_columns` | String[] | Column names: ["x", "y", "z", "s"] |
| `reco_s` | Float64[] | Arc-length of each voxel projected onto path |

##### MC Path Data

| Dataset | Type | Description |
|---------|------|-------------|
| `mc_path_data` | Matrix | MC path (P × 5): x, y, z, energy, s |
| `mc_path_columns` | String[] | Column names |

##### KDE Data

| Dataset | Type | Description |
|---------|------|-------------|
| `kde_s` | Float64[] | RECO KDE evaluation points (0 to track_length) |
| `reco_kde_f` | Float64[] | RECO energy density f(s) |
| `mc_kde_s` | Float64[] | MC KDE evaluation points |
| `mc_kde_f` | Float64[] | MC energy density f(s) |
| `kde_bandwidth` | Float64 | RECO KDE bandwidth used |
| `mc_kde_bandwidth` | Float64 | MC KDE bandwidth used |

##### Extreme Distances

| Dataset | Type | Description |
|---------|------|-------------|
| `d1` | Float64 | Distance: RECO extreme 1 to matched MC |
| `d2` | Float64 | Distance: RECO extreme 2 to matched MC |

##### Blob Analysis

| Dataset | Type | Description |
|---------|------|-------------|
| `Eb1_keV` | Float64 | High-energy blob energy (keV) |
| `Eb2_keV` | Float64 | Low-energy blob energy (keV) |
| `asymmetry` | Float64 | Blob energy asymmetry |
| `blob1_pos` | Float64[3] | Blob1 position (x, y, z) |
| `blob2_pos` | Float64[3] | Blob2 position (x, y, z) |
| `blob_radius` | Float64 | Blob radius used (mm) |

##### KDE Peak Analysis

| Dataset | Type | Description |
|---------|------|-------------|
| `n_peaks` | Int | Number of peaks found |
| `peak1_left` | Float64 | Peak1 left edge (mm) |
| `peak1_right` | Float64 | Peak1 right edge (mm) |
| `peak1_prom` | Float64 | Peak1 normalized prominence |
| `peak2_left` | Float64 | Peak2 left edge (mm) |
| `peak2_right` | Float64 | Peak2 right edge (mm) |
| `peak2_prom` | Float64 | Peak2 normalized prominence |

### 3. Metadata File: `<outdir>/<outbase>_metadata.csv`

Run parameters and processing statistics.

| Parameter | Description |
|-----------|-------------|
| `input_file` | Source file |
| `particle_type` | ion or electron |
| `ldrft_cm` | Drift length |
| `sigma_t_mm` | Transverse diffusion |
| `sigma_l_mm` | Longitudinal diffusion |
| `voxel_size_mm` | Voxel size |
| `max_distance_mm` | Max connection distance |
| `energy_threshold_keV` | Energy threshold |
| `kde_bandwidth_mm` | KDE bandwidth |
| `blob_radius_mm` | Blob radius |
| `n_events_processed` | Total events attempted |
| `n_single_track` | Events with exactly 1 track |
| `nthreads` | Threads used |
| `timestamp` | Completion time |

## Examples

### Basic usage with 8 threads

```bash
julia -t 8 scripts/itaca_track_reco_mt.jl \
    /data/HD5t/itaca/ \
    bb0nu_15bar.h5 \
    --outdir=results/bb0nu \
    --outbase=analysis \
    --nthreads=8 \
    --ievt=1 --levt=1000
```

### Ion readout at 50 cm drift

```bash
julia -t 4 scripts/itaca_track_reco_mt.jl \
    /data/HD5t/itaca/ \
    bb0nu_15bar.h5 \
    --outdir=results/ion_study \
    --outbase=ion_50cm \
    --nthreads=4 \
    --particle=ion \
    --ldrft=50.0 \
    --Rb=12.0
```

### Electron readout

```bash
julia -t 8 scripts/itaca_track_reco_mt.jl \
    /data/HD5t/itaca/ \
    bb0nu_15bar.h5 \
    --outdir=results/electron \
    --outbase=run01 \
    --nthreads=8 \
    --particle=electron \
    --Pbar=15.0
```

#### Output structure for this example

Given `--outdir=results/electron`, `--outbase=run01`, and `--nthreads=8`, the following files are created:

```
results/electron/                       # --outdir (created if needed)
├── run01_analysis.csv                  # --outbase + _analysis.csv
├── run01_metadata.csv                  # --outbase + _metadata.csv
├── run01_th_1.h5                       # --outbase + _th_N.h5
├── run01_th_2.h5
├── run01_th_3.h5
├── run01_th_4.h5
├── run01_th_5.h5
├── run01_th_6.h5
├── run01_th_7.h5
└── run01_th_8.h5
```

**CSV file** (`results/electron/run01_analysis.csv`):
```csv
event,thread_id,track_length_mm,confidence,Eb1_keV,Eb2_keV,asymmetry,d1_mm,d2_mm,n_peaks,peak1_left,peak1_right,peak1_prom,peak2_left,peak2_right,peak2_prom
42,1,198.5,0.95,412.3,389.1,0.029,2.1,3.4,2,0.0,12.5,0.85,185.2,198.5,0.91
57,1,201.2,0.93,425.0,401.8,0.028,1.8,2.9,2,0.0,11.8,0.82,189.4,201.2,0.88
...
```

**HDF5 file structure** (`results/electron/run01_th_1.h5`):
```
run01_th_1.h5
│
├── [attributes]                        # File-level metadata
│   ├── input_file = "bb0nu_15bar.h5"
│   ├── thread_id = 1
│   ├── first_event_processed = 1
│   ├── last_event_processed = 125      # (if 1000 events / 8 threads)
│   ├── events_processed = 125
│   ├── total_tracks_saved = 98         # Single-track events
│   ├── ldrift_cm = 100.0
│   ├── sigma_t_mm = 3.5                # Electron diffusion
│   ├── sigma_l_mm = 0.9
│   ├── voxel_size_mm = 7.0
│   ├── max_distance_mm = 10.5
│   ├── blob_radius_mm = 10.0
│   └── ...
│
└── batch_1/
    ├── track_1/
    │   ├── event_id = 42
    │   ├── track_length = 198.5
    │   ├── confidence = 0.95
    │   │
    │   ├── voxels = [N × 4 matrix]     # x, y, z, energy
    │   ├── voxel_columns = ["x", "y", "z", "energy"]
    │   ├── n_vertices = 156
    │   ├── graph_edges = [E × 2 matrix]
    │   ├── components = [1 × 156 matrix]
    │   │
    │   ├── path_data = [M × 4 matrix]  # x, y, z, s
    │   ├── path_columns = ["x", "y", "z", "s"]
    │   ├── reco_s = [N floats]         # Voxel arc-lengths
    │   │
    │   ├── mc_path_data = [P × 5 matrix]
    │   ├── mc_path_columns = ["x", "y", "z", "energy", "s"]
    │   │
    │   ├── kde_s = [200 floats]        # 0 to track_length
    │   ├── reco_kde_f = [200 floats]   # Energy density
    │   ├── mc_kde_s = [200 floats]
    │   ├── mc_kde_f = [200 floats]
    │   ├── kde_bandwidth = 14.0
    │   ├── mc_kde_bandwidth = 14.0
    │   │
    │   ├── d1 = 2.1
    │   ├── d2 = 3.4
    │   │
    │   ├── Eb1_keV = 412.3
    │   ├── Eb2_keV = 389.1
    │   ├── asymmetry = 0.029
    │   ├── blob1_pos = [x1, y1, z1]
    │   ├── blob2_pos = [x2, y2, z2]
    │   ├── blob_radius = 10.0
    │   │
    │   ├── n_peaks = 2
    │   ├── peak1_left = 0.0
    │   ├── peak1_right = 12.5
    │   ├── peak1_prom = 0.85
    │   ├── peak2_left = 185.2
    │   ├── peak2_right = 198.5
    │   └── peak2_prom = 0.91
    │
    ├── track_2/
    │   └── ... (event 57)
    │
    └── track_98/
        └── ... (last single-track event for thread 1)
```

**Metadata file** (`results/electron/run01_metadata.csv`):
```csv
parameter,value
input_file,bb0nu_15bar.h5
particle_type,electron
ldrft_cm,100.0
sigma_t_mm,3.5
sigma_l_mm,0.9
voxel_size_mm,7.0
max_distance_mm,10.5
energy_threshold_keV,0.25
kde_bandwidth_mm,14.0
blob_radius_mm,10.0
n_events_processed,1000
n_single_track,784
nthreads,8
timestamp,2025-12-10T15:30:00
```

### Verbose output for debugging

```bash
julia -t 2 scripts/itaca_track_reco_mt.jl \
    /data/HD5t/itaca/ \
    bb0nu_15bar.h5 \
    --outdir=results/debug \
    --outbase=test \
    --nthreads=2 \
    --ievt=1 --levt=10 \
    --print_level=verbose
```

## Dependencies

- **Petit module** (`src/Petit.jl`)
- **itaca_aux.jl** (`src/itaca_aux.jl`) - for `kde_peaks()`
- **Packages**: HDF5, DataFrames, CSV, Statistics, Graphs, Dates

## Comparison with Parent Scripts

| Feature | `track_reco_mt.jl` | `itaca_single_track_analysis.jl` | `itaca_track_reco_mt.jl` |
|---------|-------------------|----------------------------------|--------------------------|
| Multi-threading | Yes | No | Yes |
| Diffusion params | Fixed formula | `get_sigma()` by particle type | `get_sigma()` by particle type |
| Voxel sizing | Fixed scale | Dynamic by drift length | Dynamic by drift length |
| Blob analysis | No | Yes | Yes |
| KDE peaks | No | Yes | Yes |
| CSV summary | No | Yes | Yes |
| HDF5 tracks | Yes | No | Yes |
| Plots | No | Optional | No |
| Default print | verbose | verbose | quiet |

## See Also

- `track_reco_mt.jl` - Multi-threaded reconstruction (provides MT + HDF5 pattern)
- `itaca_single_track_analysis.jl` - Single-threaded with plots (provides all analysis logic)
- `pluto/plot_analysis.jl` - Functions to analyze the CSV output
