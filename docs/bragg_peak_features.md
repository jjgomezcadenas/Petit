# Bragg Peak Discrimination Features

## Problem Statement
Distinguish **signal** (2 electrons → 2 Bragg peaks) from **background** (1 electron → 1 Bragg peak) in double-beta decay experiments.

## Current Features (in test_kde_track.jl)
- Integrated energy at each endpoint (E_start, E_end)
- Energy concentration (max bin / total)
- FWHM of KDE peaks
- Number of significant bins

---

## Proposed New Features

### 1. Asymmetry Metrics
```
Energy asymmetry:      A_E = (E_max - E_min) / (E_max + E_min)
FWHM asymmetry:        A_FWHM = (FWHM_max - FWHM_min) / (FWHM_max + FWHM_min)
Peak height asymmetry: A_h = (h_max - h_min) / (h_max + h_min)
```
- All values in range [0, 1]
- **Single Bragg**: A → 1 (highly asymmetric)
- **Double Bragg**: A → 0 (symmetric)

### 2. Peak Shape Features
| Feature | Description | Single vs Double |
|---------|-------------|------------------|
| **Skewness** at each endpoint | Asymmetry of local distribution | Different signatures |
| **Kurtosis** | Peakedness (sharp vs flat) | Bragg = high kurtosis |
| **Curvature at peak** | -f''(s_max) | Sharp peak = high curvature |
| **Rise/fall ratio** | Slope before vs after peak | Bragg has characteristic shape |

### 3. Two-Peak Discriminants
```
Min peak height:      h_min = min(h_start, h_end)
Peak product:         P = h_start × h_end
Geometric mean:       √(E_start × E_end)
```
- **Double Bragg**: Both peaks significant → high h_min, high P
- **Single Bragg**: One peak weak → low h_min, low P

### 4. Valley Analysis (between endpoints)
```
Valley depth:         V = min(f(s)) for s ∈ [rblob, L - rblob]
Peak-to-valley:       PV = (h_start + h_end) / (2 × V)
Central flatness:     σ of f(s) in middle 50% of track
```
- **Double Bragg**: Deep valley, high PV ratio
- **Single Bragg**: Gradual decline, lower PV

### 5. Energy Fraction Features
```
Blob fraction:        F_blob = (E_start + E_end) / E_total
Endpoint dominance:   D = max(E_start, E_end) / E_total
Balance:              B = min(E_start, E_end) / max(E_start, E_end)
```
- **Double Bragg**: B → 1, moderate D
- **Single Bragg**: B → 0, high D

### 6. Gradient Features
```
Max gradient:         |df/ds|_max
Gradient at endpoints: df/ds at s=0 and s=L
Integral of |df/ds|:  Total variation of f(s)
```

### 7. Statistical Measures
```
Entropy:              H = -Σ p_i log(p_i)  where p_i = E_i / E_total
Gini coefficient:     Energy inequality along track
Chi² to Bragg model:  Fit quality to theoretical curve
```

### 8. Composite Score (ML-ready)
```
S_2β = w₁(1-|A_E|) + w₂·h_min + w₃·B + w₄/PV + ...
```

---

## Implementation Priority

### Phase 1 (Core discriminants)
1. **Asymmetry metrics** (A_E, A_FWHM, A_h) ✅ IMPLEMENTED
2. **Max bin analysis** ✅ IMPLEMENTED
   - Max bin energy at each endpoint
   - Next-to-max bin energy
   - Ratio max/next-to-max (higher = more concentrated = Bragg-like)
3. **Balance ratio** (B)
4. **Peak-to-valley ratio** (PV) ✅ IMPLEMENTED
   - Valley depth V = min(f(s)) in middle region
   - PV = (h_start + h_end) / (2 × V)
   - Central flatness (σ of middle 50%)
   - High PV → deep valley → double Bragg signature
5. **Minimum peak height** (h_min)

### Phase 2 (Shape analysis)
5. Curvature at peaks
6. **Skewness/Kurtosis at endpoints** ✅ IMPLEMENTED
   - Weighted skewness γ₁ = E[(X-μ)³] / σ³
   - Weighted excess kurtosis γ₂ = E[(X-μ)⁴] / σ⁴ - 3
   - High kurtosis → sharp peak (Bragg-like)
   - Positive skewness → energy tails toward track center
7. Gradient features

### Phase 3 (Advanced)
8. **Entropy** ✅ IMPLEMENTED
   - H = -Σ pᵢ log(pᵢ) where pᵢ = Eᵢ / E_total
   - Computed for full track, start region, end region
   - Normalized entropy H_norm ∈ [0,1] (0=concentrated, 1=uniform)
   - Entropy asymmetry A_H for signal/background discrimination
9. Chi² fit to Bragg model
10. Composite ML score

---

## Related Files
- `scripts/test_kde_track.jl` - Current test script
- `src/kde.jl` - KDE implementation
- `src/itaca_functions.jl` - Track analysis functions
