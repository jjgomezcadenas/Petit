# ITACA Single Track Analysis

## Overview

`itaca_single_track_analysis.jl` is a command-line Julia script for analyzing single particle tracks in a xenon gas TPC detector. It processes Monte Carlo events, simulates detector response (diffusion), reconstructs tracks, and performs blob/KDE analysis to characterize energy deposits at track extremes.

This analysis is relevant for neutrinoless double beta decay (0νββ) searches where identifying the topology of electron tracks is crucial for background rejection.

## Location

```
scripts/itaca_single_track_analysis.jl
```

## Dependencies

- **Petit module** (`src/Petit.jl`) - Core physics and analysis functions
- **itaca_aux.jl** (`src/itaca_aux.jl`) - Helper functions for printing/displaying results
- **External packages**: Printf, DataFrames, CSV, Glob, ArgParse, Plots, Dates

## Data Requirements

Input data must be HDF5 files located in `$DATA/HD5t/itaca/`. The default expects:
```
$DATA/HD5t/itaca/bb0nu/bb0nu_15bar_p1.h5
```

## Processing Pipeline

The script performs the following steps for each event:

### 1. Diffusion Parameter Calculation
Computes transverse (σ_t) and longitudinal (σ_l) diffusion based on:
- **Particle type**: `ion` or `electron`
- **Drift length** (cm)
- **Detector parameters**: temperature (K), drift field (V/cm), pressure (bar)

For ions: uses thermal diffusion model (`sigma_t_ion_mm`)
For electrons: uses parameterized diffusion (`sigma_t_mm`, `sigma_l_mm`)

### 2. Voxel Size Determination
Dynamically calculates voxel size based on diffusion:
- For drift > 50 cm: `voxel_size = 2.0 × σ_t`, `max_distance = 1.5 × voxel_size`
- For drift ≤ 50 cm: `voxel_size = 3.0 × σ_t`, `max_distance = 2.0 × voxel_size`

### 3. MC Path Computation
Computes the "true" Monte Carlo path through the event hits for comparison with reconstruction.

### 4. Hit Transformation and Diffusion
- Transforms MC hits to detector coordinates
- Applies Gaussian diffusion to simulate detector response
- Creates a diffused image of the event

### 5. Voxelization
Groups diffused hits into 3D voxels of the computed size.

### 6. Track Making
Uses Dijkstra-based graph traversal to identify connected tracks:
- Applies energy threshold cut
- Uses `max_distance` for voxel connectivity
- Selects only events with exactly **one track**

### 7. Track Walking
`walk_track_from_extremes()` identifies:
- Track extremes (endpoints)
- Optimal path through the track
- Path confidence metric

### 8. Extreme Distance Calculation
Compares reconstructed (RECO) track extremes with MC truth:
- Computes d₁ and d₂ (distances between RECO and MC endpoints)
- Uses optimal pairing (direct or crossed) to minimize total distance

### 9. KDE Analysis
Kernel Density Estimation of energy distribution along the track:
- Computes longitudinal energy profile f(s) vs arc length s
- Bandwidth: `2.0 × voxel_size`
- Evaluates on `n_kde_eval` points (default: 200)

### 10. Peak Detection
Identifies peaks in the KDE energy profile:
- Uses prominence-based peak finding
- Extracts leftmost (peak1) and rightmost (peak2) peaks
- Records left edge, right edge, and normalized prominence

### 11. Blob Analysis
Finds energy deposits within spheres at track extremes:
- Uses configurable blob radius (default: 10 mm)
- Computes blob energies Eb1 (high) and Eb2 (low)
- Calculates asymmetry: `(Eb1 - Eb2) / (Eb1 + Eb2)`

## Command Line Interface

```bash
julia scripts/itaca_single_track_analysis.jl [OPTIONS]
```

### Event Selection
| Option | Description | Default |
|--------|-------------|---------|
| `-i, --ievent` | Initial event number | 1 |
| `-e, --levent` | Last event number | 1 |
| `-f, --input` | Input file (relative to DATA/HD5t/itaca) | `bb0nu/bb0nu_15bar_p1.h5` |

### Physics Parameters
| Option | Description | Default |
|--------|-------------|---------|
| `-p, --particle` | Particle type: `ion` or `electron` | `ion` |
| `-l, --ldrft` | Drift length (cm) | 100.0 |
| `--dt` | Transverse diffusion (mm/√cm) | 3.5 |
| `--dl` | Longitudinal diffusion (mm/√cm) | 0.9 |
| `--tK` | Temperature (K) | 297.0 |
| `--edrift` | Drift field (V/cm) | 500.0 |
| `--Pbar` | Pressure (bar) | 15.0 |

### Analysis Parameters
| Option | Description | Default |
|--------|-------------|---------|
| `--eth-ion` | Energy threshold for ions (keV) | 10.0 |
| `--eth` | Energy threshold (keV) | 10.0 |
| `--nbins` | Number of bins for diffusion | 100 |
| `--nsigma` | Number of sigma for diffusion | 3.0 |
| `--nkde` | Number of KDE evaluation points | 200 |
| `--Rb` | Blob radius (mm) | 10.0 |

### Output Options
| Option | Description | Default |
|--------|-------------|---------|
| `-o, --outdir` | Output directory | `results` |
| `--interactive` | Enable interactive mode (display + wait) | false |
| `--nplot` | Save plots every N events (0 = disabled) | 0 |
| `--print_level` | Output verbosity: `quiet`, `verbose`, or `very_verbose` | `verbose` |

#### Print Level Modes

- **`very_verbose`**: Full detailed output including:
  - Diffusion parameters at startup
  - Processing step headers for each event (`### Computing MC path`, `### Transforming hits`, etc.)
  - Number of tracks found
  - Peak information with details
  - Pretty-printed results table for each event

- **`verbose`** (default): Moderate output showing:
  - Diffusion parameters at startup
  - Event number with track count (e.g., `Event 5: 1 track` or `Event 5: 2 tracks - rejected`)
  - Pretty-printed results table for each successful event

- **`quiet`**: Minimal output for batch processing:
  - Only prints event number every `--nplot` events (together with plot saving)
  - Shows rejection message if event has multiple tracks
  - Prints "done" when event completes successfully

## Output Files

### analysis_results.csv
Main results file with columns:

| Column | Description | Units |
|--------|-------------|-------|
| `event` | Event number | - |
| `track_length_mm` | Reconstructed track length | mm |
| `Eb1_keV` | High-energy blob energy | keV |
| `Eb2_keV` | Low-energy blob energy | keV |
| `asymmetry` | Blob energy asymmetry | - |
| `d1_mm` | RECO-MC distance (extreme 1) | mm |
| `d2_mm` | RECO-MC distance (extreme 2) | mm |
| `n_peaks` | Number of KDE peaks found | - |
| `peak1_left` | Leftmost peak left edge | mm |
| `peak1_right` | Leftmost peak right edge | mm |
| `peak1_prom` | Leftmost peak prominence (normalized) | - |
| `peak2_left` | Rightmost peak left edge | mm |
| `peak2_right` | Rightmost peak right edge | mm |
| `peak2_prom` | Rightmost peak prominence (normalized) | - |

### metadata.csv
Run parameters and processing statistics:
- All input parameters
- `n_events_processed`: Total events attempted
- `n_events_single_track`: Events with exactly one track
- `timestamp`: Analysis completion time

### plots/ (optional)
When `--nplot N` is set (N > 0), saves PNG plots every N events:
1. `event_X_01_mc_hits.png` - Original MC hits
2. `event_X_02_mc_path.png` - MC path overlay
3. `event_X_03_diffused.png` - Diffused event
4. `event_X_04_voxels.png` - Voxelized event
5. `event_X_05_track_paths.png` - RECO vs MC paths
6. `event_X_06_kde_comparison.png` - KDE and histogram comparison
7. `event_X_07_kde_peaks.png` - KDE with peak markers
8. `event_X_08_blobs.png` - Track with blob spheres

## Usage Examples

### Basic analysis of 10 events
```bash
julia scripts/itaca_single_track_analysis.jl -i 1 -e 10 -o results/test
```

### Ion analysis at 50 cm drift
```bash
julia scripts/itaca_single_track_analysis.jl \
    -i 1 -e 100 \
    -p ion -l 50.0 \
    -o results/ion_50cm
```

### Electron analysis with plots
```bash
julia scripts/itaca_single_track_analysis.jl \
    -i 1 -e 50 \
    -p electron \
    --nplot 10 \
    -o results/electron_with_plots
```

### Interactive debugging (single event)
```bash
julia scripts/itaca_single_track_analysis.jl \
    -i 5 -e 5 \
    --interactive \
    -o results/debug
```

### Custom blob radius and energy threshold
```bash
julia scripts/itaca_single_track_analysis.jl \
    -i 1 -e 100 \
    --Rb 15.0 \
    --eth 5.0 \
    -o results/custom_params
```

### Quiet mode batch processing with periodic plots
```bash
julia scripts/itaca_single_track_analysis.jl \
    -i 1 -e 1000 \
    --print_level quiet \
    --nplot 100 \
    -o results/batch_run
```

### Very verbose mode for debugging
```bash
julia scripts/itaca_single_track_analysis.jl \
    -i 1 -e 10 \
    --print_level very_verbose \
    -o results/debug_verbose
```

## Key Physics Concepts

### Blob Energy Asymmetry
For 0νββ events, both electrons deposit energy at their endpoints (Bragg peaks). The asymmetry metric helps distinguish signal from background:
- **Signal (0νββ)**: Two similar blobs → low asymmetry
- **Background (single electron)**: One blob much larger → high asymmetry

### KDE Peak Analysis
The longitudinal energy profile reveals Bragg peak structure:
- **Two prominent peaks**: Characteristic of 0νββ topology
- **Peak prominence**: Indicates how well-defined the Bragg peak is
- **Peak width**: Related to diffusion and track complexity

### Extreme Distance Metrics
Measures reconstruction quality:
- **d₁, d₂**: Distance between RECO and MC track endpoints
- **Small values**: Good reconstruction fidelity
- **Large values**: May indicate reconstruction issues or complex topologies

## Helper Functions

### Script-local functions (in `itaca_single_track_analysis.jl`)

| Function | Purpose |
|----------|---------|
| `get_sigma()` | Compute diffusion parameters |
| `get_voxel_size_and_distance()` | Determine voxelization parameters |
| `get_energy_threshold()` | Get energy threshold by particle type |
| `print_event_results()` | Pretty-print CSV row contents (verbose mode) |
| `print_blobs()` | Print blob analysis results |
| `print_peaks()` | Print peak analysis results |
| `plot_kde_comparison()` | Generate KDE comparison plot |
| `save_plot()` | Save plot to file |
| `save_and_continue()` | Save plot and wait for user (interactive) |
| `save_metadata()` | Save run parameters to CSV |

### Functions from `itaca_aux.jl`

| Function | Purpose |
|----------|---------|
| `kde_peaks()` | Extract peak1/peak2 from KDE analysis |
| `diffusion_params_print()` | Print diffusion parameters |
| `walk_result_print()` | Print walk result summary |
| `path_print()` | Print path summary |
| `extreme_distances_print()` | Print extreme distances |

## See Also

- `batch_itaca_analysis_mt.jl` - Multi-threaded batch version
- `itaca_resolution_single_track.jl` - Resolution studies
- `src/kde.jl` - KDE implementation
- `src/itaca_functions.jl` - Core ITACA functions
