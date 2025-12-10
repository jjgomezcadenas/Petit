# Summary: `test_kde_track.jl`

## Purpose

Interactive analysis script for testing Bragg peak discrimination using Kernel Density Estimation (KDE) along the reconstructed central path of particle tracks. Computes multiple features to distinguish signal (2 electrons → 2 Bragg peaks) from background (1 electron → 1 Bragg peak).

---

## Analysis Pipeline

```text
Reco HDF5 → Load Track → Project onto Path → KDE → Feature Extraction → Plots
```

### Steps

1. **Load reconstructed track** from HDF5 (output of `track_reco_mt.jl`)
2. **Analyze energy profile** using `Petit.analyze_track_energy_profile()`
3. **Compute discrimination features** at both endpoints
4. **Generate diagnostic plots** (15 subplots)

---

## Input

Reconstructed track HDF5 file from `track_reco_mt.jl` containing:

- `track.voxels` - Voxel DataFrame with energy
- `central_path` - Smoothed path points with arc-length
- `mc_path` - MC truth path (if available)

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `bandwidth` | 5.0 mm | KDE smoothing bandwidth |
| `rblob` | 15.0 mm | Blob integration radius at endpoints |
| `binwidth` | 3.0 mm | Histogram bin width for energy analysis |

---

## Computed Features

### Phase 1: Core Discriminants

| Feature | Range | Signal (2β) | Background (1e) |
|---------|-------|-------------|-----------------|
| **Energy asymmetry A_E** | [0, 1] | → 0 (symmetric) | → 1 (asymmetric) |
| **FWHM asymmetry A_FWHM** | [0, 1] | → 0 | → 1 |
| **Peak height asymmetry A_h** | [0, 1] | → 0 | → 1 |
| **Peak-to-Valley ratio PV** | [1, ∞] | High (deep valley) | Low |
| **Max bin energy** | keV | Similar at both ends | One end dominant |
| **Max/next-to-max ratio** | [1, ∞] | Moderate | High at one end |

### Phase 2: Shape Analysis

| Feature | Description | Bragg Signature |
|---------|-------------|-----------------|
| **Skewness γ₁** | Distribution asymmetry | Positive (tail toward center) |
| **Kurtosis γ₂** | Peakedness (excess) | High (>0, sharper than normal) |

### Phase 3: Statistical Measures

| Feature | Range | Concentrated | Spread |
|---------|-------|--------------|--------|
| **Entropy H** | [0, log(n)] | Low | High |
| **Normalized entropy** | [0, 1] | → 0 | → 1 |
| **Entropy asymmetry A_H** | [-1, 1] | ≈ 0 (symmetric) | ≠ 0 |

---

## Feature Formulas

### Asymmetry Metrics

```text
A = (max - min) / (max + min)   for A_E, A_FWHM, A_h
```

All asymmetries in [0, 1]: 0 = symmetric, 1 = completely asymmetric

### Peak-to-Valley Ratio

```text
V = min(f(s))  for s ∈ [rblob, L - rblob]   (valley depth)
PV = (h_start + h_end) / (2 × V)
```

High PV indicates deep valley between two peaks (double Bragg signature)

### Weighted Skewness and Kurtosis

```text
μ = Σ(sᵢ × Eᵢ) / Σ Eᵢ                    (energy-weighted mean)
σ² = Σ(Eᵢ × (sᵢ - μ)²) / Σ Eᵢ            (energy-weighted variance)
γ₁ = Σ(Eᵢ × (sᵢ - μ)³) / (Σ Eᵢ × σ³)    (skewness)
γ₂ = Σ(Eᵢ × (sᵢ - μ)⁴) / (Σ Eᵢ × σ⁴) - 3 (excess kurtosis)
```

### Entropy

```text
pᵢ = Eᵢ / Σ Eᵢ                           (energy probability)
H = -Σ pᵢ log(pᵢ)                         (Shannon entropy)
H_norm = H / log(n_bins)                  (normalized to [0,1])
```

---

## Output

### Console Output

- Track statistics (voxels, path points, length)
- Detected KDE maxima with prominences
- Integrated energy at each endpoint
- All computed features with interpretation

### Plot Output

`kde_track_test.png` saved in the same directory as input file.

15 subplots organized in 9 rows:

| Row | Content |
|-----|---------|
| 1 | KDE energy density f(s) with peak markers |
| 2 | Energy histogram vs arc-length |
| 3 | Integrated blob energy + FWHM comparison |
| 4 | Start/End region energy histograms |
| 5 | Asymmetry metrics (A_E, A_FWHM, A_h) |
| 6 | Max bin energies + Max/next ratio |
| 7 | Peak-to-Valley visualization + PV ratio |
| 8 | Skewness + Kurtosis comparison |
| 9 | Entropy + Entropy asymmetry |

---

## Usage

```bash
julia test_kde_track.jl <reco_file.h5> [options]
```

### Command Line Options

```text
Required:
  reco_file.h5       Reconstructed track file (from track_reco_mt.jl)

Optional:
  --event=N          Event index (1-based, default: 1)
  --bandwidth=H      KDE bandwidth in mm (default: 5.0)
  --rblob=R          Blob radius in mm (default: 15.0)
  --binwidth=W       Histogram bin width in mm (default: 3.0)
```

### Examples

```bash
# Analyze first event with defaults
julia test_kde_track.jl output_th_1.h5

# Analyze event 5 with custom parameters
julia test_kde_track.jl output_th_1.h5 --event=5 --bandwidth=3.0 --rblob=20.0
```

---

## Key Functions

| Function | Description |
|----------|-------------|
| `main()` | Main analysis and plotting function |
| `compute_fwhm()` | Calculate FWHM of KDE peak |
| `get_max_and_next()` | Extract max and next-to-max bin values |
| `compute_skew_kurt()` | Energy-weighted skewness and kurtosis |
| `compute_entropy()` | Shannon entropy from energy distribution |

### External (from Petit)

| Function | Description |
|----------|-------------|
| `Petit.read_reco_results_from_hdf5()` | Load reconstructed tracks |
| `Petit.analyze_track_energy_profile()` | KDE analysis along path |

---

## Interpretation Guide

### Double Bragg (Signal)

- Low asymmetries: A_E, A_FWHM, A_h → 0
- High PV ratio (> 2.5)
- Similar energy at both endpoints
- Similar FWHM at both endpoints
- Low, similar entropy at both endpoints

### Single Bragg (Background)

- High asymmetries: A_E, A_FWHM, A_h → 1
- Low PV ratio (< 1.5)
- Energy concentrated at one endpoint
- High max/next ratio at one endpoint
- One endpoint with low entropy, other with high

---

## Related Files

- `scripts/track_reco_mt.jl` - Produces input HDF5 files
- `src/kde.jl` - KDE implementation
- `src/itaca_functions.jl` - Track analysis functions
- `docs/bragg_peak_features.md` - Feature documentation
