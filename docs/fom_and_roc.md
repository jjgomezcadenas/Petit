# Figure of Merit and ROC Analysis Module

**File**: `src/fom_and_roc.jl`
**Module**: `Petit`

## Overview

This module provides functions for evaluating signal/background discrimination using energy cuts. It computes efficiency curves, figure of merit optimization, and ROC curves for bb0nu (signal) vs xe137 (background) analysis.

All functions take vectors directly (not DataFrames), making them flexible and easy to use with any data source.

---

## Quick Start

```julia
using Petit
using CSV, DataFrames

# Load analysis results from ITACA output
df_signal = CSV.read("bb0nu_analysis.csv", DataFrame)
df_background = CSV.read("xe137_analysis.csv", DataFrame)

# Extract Eb2 vectors
eb2_sig = df_signal.Eb2
eb2_bkg = df_background.Eb2

# Quick summary - prints optimal cut and metrics
print_fom_summary(eb2_sig, eb2_bkg; signal_label="bb0nu", background_label="xe137")

# Generate all plots in one call
p = plot_all_fom_analysis(eb2_sig, eb2_bkg; xlabel="Eb2 (keV)")
savefig(p, "fom_analysis.png")
```

---

## Typical Workflow

### Step 1: Load Data

```julia
using Petit
using CSV, DataFrames

# Load ITACA analysis results
df_bb0nu = CSV.read("/path/to/bb0nu_analysis.csv", DataFrame)
df_xe137 = CSV.read("/path/to/xe137_analysis.csv", DataFrame)

# Check available columns
println(names(df_bb0nu))
# [:event_id, :track_length_mm, :Eb1, :Eb2, :asymmetry, :n_peaks, ...]
```

### Step 2: Extract Discriminating Variable

```julia
# For Eb2-based discrimination
eb2_sig = df_bb0nu.Eb2
eb2_bkg = df_xe137.Eb2

println("Signal events: $(length(eb2_sig))")
println("Background events: $(length(eb2_bkg))")
```

### Step 3: Find Optimal Cut

```julia
# Compute FOM and find optimal cut
fom = figure_of_merit(eb2_sig, eb2_bkg; cuts=range(0, 800, length=201))

println("Optimal Eb2 cut: $(fom.optimal_cut) keV")
println("Signal efficiency: $(round(100*fom.optimal_signal_eff, digits=1))%")
println("Background efficiency: $(round(100*fom.optimal_background_eff, digits=1))%")
println("FOM: $(round(fom.optimal_fom, digits=3))")
```

### Step 4: Evaluate ROC Performance

```julia
# Compute ROC curve
roc = roc_curve(eb2_sig, eb2_bkg; cuts=range(0, 800, length=201))

println("ROC AUC: $(round(roc.auc, digits=4))")

# AUC interpretation:
# - AUC = 1.0: Perfect separation
# - AUC = 0.9: Excellent discrimination
# - AUC = 0.8: Good discrimination
# - AUC = 0.5: No discrimination (random)
```

### Step 5: Generate Plots

```julia
# Individual plots
p1 = plot_efficiency_vs_cut(eb2_sig, eb2_bkg; xlabel="Eb2 (keV)")
savefig(p1, "efficiency_vs_eb2.png")

p2 = plot_fom_vs_cut(eb2_sig, eb2_bkg; xlabel="Eb2 (keV)")
savefig(p2, "fom_vs_eb2.png")

p3 = plot_roc(eb2_sig, eb2_bkg)
savefig(p3, "roc_curve.png")

# Or all at once in 2x2 layout
p_all = plot_all_fom_analysis(eb2_sig, eb2_bkg;
                               xlabel="Eb2 (keV)",
                               signal_label="bb0nu",
                               background_label="xe137")
savefig(p_all, "fom_analysis_complete.png")
```

---

## Functions Reference

### 1. `efficiency_vs_cut`

Computes efficiency as a function of cut threshold for a single dataset.

```julia
efficiency_vs_cut(data::AbstractVector{<:Real};
                  cuts::AbstractVector{<:Real} = range(0, 1000, length=101)) -> NamedTuple
```

**Parameters:**

- `data`: Vector of values (e.g., Eb2 energies)
- `cuts`: Vector of threshold values to scan (default: 0 to 1000 in 101 steps)

**Returns:** NamedTuple with fields:

- `cuts`: Vector of cut values
- `efficiency`: Vector of efficiencies (fraction of events passing data < cut)
- `n_total`: Total number of events
- `n_passing`: Vector of events passing each cut

**Example:**

```julia
# Single dataset efficiency scan
result = efficiency_vs_cut(df_bb0nu.Eb2; cuts=range(0, 600, length=61))

# Plot efficiency curve
plot(result.cuts, result.efficiency,
     xlabel="Eb2 (keV)", ylabel="Efficiency",
     title="bb0nu Efficiency vs Eb2 Cut")
```

---

### 2. `compute_efficiencies`

Computes signal and background efficiencies for a range of cuts.

```julia
compute_efficiencies(signal::AbstractVector{<:Real},
                     background::AbstractVector{<:Real};
                     cuts::AbstractVector{<:Real} = range(0, 1000, length=101)) -> NamedTuple
```

**Returns:** NamedTuple with fields:

- `cuts`: Vector of cut values
- `signal_eff`: Vector of signal efficiencies
- `background_eff`: Vector of background efficiencies
- `n_signal`: Total signal events
- `n_background`: Total background events

**Example:**

```julia
eff = compute_efficiencies(eb2_sig, eb2_bkg; cuts=range(0, 800, length=81))

# Find cut for 90% signal efficiency
idx_90 = findfirst(e -> e >= 0.90, eff.signal_eff)
cut_90 = eff.cuts[idx_90]
bkg_at_90 = eff.background_eff[idx_90]
println("At 90% signal efficiency: cut=$(cut_90) keV, bkg_eff=$(round(100*bkg_at_90, digits=1))%")
```

---

### 3. `figure_of_merit`

Computes figure of merit (FOM = signal_eff / sqrt(background_eff)) vs cut.

```julia
figure_of_merit(signal::AbstractVector{<:Real},
                background::AbstractVector{<:Real};
                cuts::AbstractVector{<:Real} = range(0, 1000, length=101),
                epsilon::Float64 = 1e-10) -> NamedTuple
```

**Returns:** NamedTuple with fields:

- `cuts`: Vector of cut values
- `fom`: Vector of FOM values
- `signal_eff`, `background_eff`: Efficiency vectors
- `optimal_cut`: Cut value that maximizes FOM
- `optimal_fom`: Maximum FOM value
- `optimal_signal_eff`, `optimal_background_eff`: Efficiencies at optimal cut

**Example:**

```julia
fom = figure_of_merit(eb2_sig, eb2_bkg)

# Access all results
println("Optimal cut: $(fom.optimal_cut)")
println("Max FOM: $(fom.optimal_fom)")
println("Signal eff at optimal: $(fom.optimal_signal_eff)")
println("Background eff at optimal: $(fom.optimal_background_eff)")

# Plot FOM curve manually
plot(fom.cuts, fom.fom, xlabel="Eb2 (keV)", ylabel="FOM")
vline!([fom.optimal_cut], linestyle=:dash, label="Optimal")
```

---

### 4. `roc_curve`

Computes ROC curve: background rejection vs signal efficiency.

```julia
roc_curve(signal::AbstractVector{<:Real},
          background::AbstractVector{<:Real};
          cuts::AbstractVector{<:Real} = range(0, 1000, length=101)) -> NamedTuple
```

**Returns:** NamedTuple with fields:

- `signal_eff`: Vector of signal efficiencies (x-axis for ROC)
- `background_rejection`: Vector of background rejection (1 - background_eff)
- `background_eff`: Vector of background efficiencies
- `cuts`: Vector of cut values (for reference)
- `auc`: Area Under the ROC Curve

**Example:**

```julia
roc = roc_curve(eb2_sig, eb2_bkg)

# AUC is the key metric
println("AUC: $(roc.auc)")

# Find operating point for specific signal efficiency
target_eff = 0.85
idx = argmin(abs.(roc.signal_eff .- target_eff))
println("At $(100*target_eff)% signal efficiency:")
println("  Background rejection: $(round(100*roc.background_rejection[idx], digits=1))%")
println("  Cut value: $(roc.cuts[idx])")
```

---

### 5. `plot_efficiency_vs_cut`

Plots signal and background efficiency vs cut.

```julia
plot_efficiency_vs_cut(signal, background;
                       cuts = range(0, 1000, length=101),
                       xlabel = "Cut Value",
                       title = "Efficiency vs Cut",
                       signal_label = "Signal",
                       background_label = "Background") -> Plot
```

**Example:**

```julia
p = plot_efficiency_vs_cut(eb2_sig, eb2_bkg;
                           xlabel="Eb2 (keV)",
                           title="Efficiency vs Eb2 Cut",
                           signal_label="bb0nu",
                           background_label="xe137")
savefig(p, "efficiency.png")
```

---

### 6. `plot_fom_vs_cut`

Plots figure of merit vs cut with optimal point marked.

```julia
plot_fom_vs_cut(signal, background;
                cuts = range(0, 1000, length=101),
                xlabel = "Cut Value",
                title = "Figure of Merit vs Cut",
                mark_optimal = true) -> Plot
```

**Example:**

```julia
# With optimal point marked (default)
p = plot_fom_vs_cut(eb2_sig, eb2_bkg; xlabel="Eb2 (keV)")

# Without optimal point marker
p = plot_fom_vs_cut(eb2_sig, eb2_bkg; xlabel="Eb2 (keV)", mark_optimal=false)
```

---

### 7. `plot_roc`

Plots ROC curve (background rejection vs signal efficiency).

```julia
plot_roc(signal, background;
         cuts = range(0, 1000, length=101),
         title = "ROC Curve",
         show_auc = true,
         show_diagonal = true) -> Plot
```

**Example:**

```julia
# Full ROC with AUC and diagonal reference
p = plot_roc(eb2_sig, eb2_bkg)

# Clean ROC without annotations
p = plot_roc(eb2_sig, eb2_bkg; show_auc=false, show_diagonal=false)
```

---

### 8. `plot_all_fom_analysis`

Creates a 2x2 layout with all analysis plots.

```julia
plot_all_fom_analysis(signal, background;
                      cuts = range(0, 1000, length=101),
                      size = (1000, 800),
                      xlabel = "Cut Value",
                      signal_label = "Signal",
                      background_label = "Background") -> Plot
```

**Layout:**

```
┌─────────────────────┬─────────────────────┐
│ Efficiency vs Cut   │  FOM vs Cut         │
├─────────────────────┼─────────────────────┤
│ ROC Curve           │  Distributions      │
└─────────────────────┴─────────────────────┘
```

**Example:**

```julia
p = plot_all_fom_analysis(eb2_sig, eb2_bkg;
                          cuts=range(0, 800, length=101),
                          xlabel="Eb2 (keV)",
                          signal_label="bb0nu",
                          background_label="xe137",
                          size=(1200, 900))
savefig(p, "complete_fom_analysis.png")
```

---

### 9. `print_fom_summary`

Prints summary statistics to console.

```julia
print_fom_summary(signal, background;
                  cuts = range(0, 1000, length=101),
                  signal_label = "Signal",
                  background_label = "Background")
```

**Example Output:**

```
═══════════════════════════════════════════════════════════════
                    FOM Analysis Summary
═══════════════════════════════════════════════════════════════
bb0nu events:       10000
xe137 events:       50000
───────────────────────────────────────────────────────────────
Optimal cut:        450.0
───────────────────────────────────────────────────────────────
At optimal cut:
  Signal efficiency:      0.8523 (85.23%)
  Background efficiency:  0.1247 (12.47%)
  Background rejection:   0.8753 (87.53%)
  Figure of Merit:        2.413
───────────────────────────────────────────────────────────────
ROC AUC:            0.9234
═══════════════════════════════════════════════════════════════
```

---

## Advanced Usage

### Using Different Discriminating Variables

The functions work with any variable, not just Eb2:

```julia
# Asymmetry-based discrimination
asym_sig = df_bb0nu.asymmetry
asym_bkg = df_xe137.asymmetry

# Note: asymmetry ranges from -1 to 1
fom_asym = figure_of_merit(asym_sig, asym_bkg; cuts=range(-1, 1, length=101))
println("Optimal asymmetry cut: $(fom_asym.optimal_cut)")

p = plot_all_fom_analysis(asym_sig, asym_bkg;
                          cuts=range(-1, 1, length=101),
                          xlabel="Asymmetry")
```

### Comparing Multiple Variables

```julia
# Compare Eb2 vs asymmetry discrimination
roc_eb2 = roc_curve(df_bb0nu.Eb2, df_xe137.Eb2)
roc_asym = roc_curve(df_bb0nu.asymmetry, df_xe137.asymmetry; cuts=range(-1, 1, length=101))

println("Eb2 AUC: $(round(roc_eb2.auc, digits=4))")
println("Asymmetry AUC: $(round(roc_asym.auc, digits=4))")

# Plot both ROC curves
p = plot(roc_eb2.signal_eff, roc_eb2.background_rejection,
         label="Eb2 (AUC=$(round(roc_eb2.auc, digits=3)))",
         xlabel="Signal Efficiency", ylabel="Background Rejection")
plot!(p, roc_asym.signal_eff, roc_asym.background_rejection,
      label="Asymmetry (AUC=$(round(roc_asym.auc, digits=3)))")
plot!(p, [0,1], [0,1], linestyle=:dash, color=:gray, label="Random")
```

### Custom Cut Ranges

```julia
# Fine scan around expected optimal region
fine_cuts = range(300, 500, length=201)
fom_fine = figure_of_merit(eb2_sig, eb2_bkg; cuts=fine_cuts)
println("Fine-tuned optimal: $(fom_fine.optimal_cut) keV")

# Coarse scan for quick overview
coarse_cuts = range(0, 1000, length=21)
p = plot_efficiency_vs_cut(eb2_sig, eb2_bkg; cuts=coarse_cuts, xlabel="Eb2 (keV)")
```

### Applying Pre-Selection Cuts

```julia
# Apply track length cut before FOM analysis
min_length = 150.0  # mm
max_length = 250.0  # mm

# Filter DataFrames
df_bb0nu_sel = filter(row -> min_length < row.track_length_mm < max_length, df_bb0nu)
df_xe137_sel = filter(row -> min_length < row.track_length_mm < max_length, df_xe137)

# Compute FOM on filtered data
fom_sel = figure_of_merit(df_bb0nu_sel.Eb2, df_xe137_sel.Eb2)
println("After track length cut:")
println("  Signal events: $(length(df_bb0nu_sel.Eb2)) ($(round(100*length(df_bb0nu_sel.Eb2)/length(df_bb0nu.Eb2), digits=1))%)")
println("  Optimal Eb2 cut: $(fom_sel.optimal_cut) keV")
```

### Extracting Working Points

```julia
eff = compute_efficiencies(eb2_sig, eb2_bkg; cuts=range(0, 800, length=801))

# Find cuts for specific signal efficiencies
for target in [0.95, 0.90, 0.85, 0.80]
    idx = findfirst(e -> e >= target, eff.signal_eff)
    if idx !== nothing
        cut = eff.cuts[idx]
        bkg = eff.background_eff[idx]
        rej = 1 - bkg
        println("Sig eff $(Int(100*target))%: cut=$(round(cut, digits=1)) keV, bkg rej=$(round(100*rej, digits=1))%")
    end
end
```

---

## Physics Background

### Cut Convention

Events pass if `value < cut_threshold`. For Eb2 (second blob energy):

- **bb0nu signal**: One electron escapes with neutrinos, so one blob has low energy
- **xe137 background**: Single electron deposits energy more uniformly
- Lower Eb2 cut → higher signal purity, lower efficiency

### Figure of Merit

The FOM formula `signal_eff / sqrt(background_eff)` is appropriate when:

- Background scales linearly with efficiency (typical for rare event searches)
- You want to optimize sensitivity to signal

Alternative FOMs (not implemented but easy to compute):

```julia
# Significance: S / sqrt(S + B)
# Purity: S / (S + B)
# Custom: modify the efficiency vectors as needed
```

### ROC AUC Interpretation

| AUC Range | Interpretation |
|-----------|----------------|
| 0.95-1.00 | Excellent discrimination |
| 0.90-0.95 | Very good |
| 0.80-0.90 | Good |
| 0.70-0.80 | Fair |
| 0.50-0.70 | Poor |
| 0.50 | No discrimination (random) |

---

## Notes

1. **Dependencies**: Functions use only `Plots.jl` (already imported by Petit module)

2. **Edge cases**:
   - Zero background efficiency: epsilon (1e-10) added to avoid division by zero
   - Empty vectors: functions return NaN/empty results

3. **Performance**: For large datasets (>100k events), consider using coarser cut ranges first, then refine around the optimal region

4. **Saving results**: All functions return NamedTuples that can be easily converted to DataFrames for saving:
   ```julia
   fom = figure_of_merit(eb2_sig, eb2_bkg)
   df_fom = DataFrame(cuts=fom.cuts, fom=fom.fom,
                      signal_eff=fom.signal_eff,
                      background_eff=fom.background_eff)
   CSV.write("fom_scan.csv", df_fom)
   ```
