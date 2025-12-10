# Note: This file is included after Petit module is loaded
# Functions use Petit.* via the including script's scope
using DataFrames
using Markdown



  """
      walk_result_md(walk_result)

  Return markdown summary of walk_result for Pluto display.
  """
  function walk_result_md(walk_result)
      if isempty(walk_result.path_indices)
          return Markdown.parse("**Walk Result:** No valid path found")
      end

      ext1, ext2 = walk_result.extremes

      md = """
      #### Walk Result Summary
      
      | Property | Value |
      |----------|-------|
      | **Path length** | $(round(walk_result.total_length, digits=2)) mm |
      | **Confidence** | $(round(walk_result.confidence, digits=3)) |
      | **N voxels in path** | $(length(walk_result.path_indices)) |
      
      ##### Extreme 1 (start)
      | x | y | z | energy |
      |---|---|---|--------|
      | $(round(ext1.x, digits=2)) | $(round(ext1.y, digits=2)) | $(round(ext1.z, digits=2)) | $(round(ext1.energy*1e+3, digits=2)) |
      
      ##### Extreme 2 (end)
      | x | y | z | energy |
      |---|---|---|--------|
      | $(round(ext2.x, digits=2)) | $(round(ext2.y, digits=2)) | $(round(ext2.z, digits=2)) | $(round(ext2.energy*1e+3, digits=2)) |
      """

      return Markdown.parse(md)
  end


  """
      path_md(path::DataFrame)

  Return markdown summary of path DataFrame for Pluto display.
  """
  function path_md(path::DataFrame)
      if nrow(path) == 0
          return Markdown.parse("**Path:** Empty")
      end

      # Basic stats
      track_length = path.s[end]
      n_points = nrow(path)

      # Spatial extent
      dx = maximum(path.x) - minimum(path.x)
      dy = maximum(path.y) - minimum(path.y)
      dz = maximum(path.z) - minimum(path.z)

      # Start and end positions
      p_start = path[1, :]
      p_end = path[end, :]

      # Average step size
      avg_step = n_points > 1 ? track_length / (n_points - 1) : 0.0

      md = """
      #### Path Summary
      
      | Property | Value |
      |----------|-------|
      | **Track length** | $(round(track_length, digits=2)) mm |
      | **N points** | $(n_points) |
      | **Avg step** | $(round(avg_step, digits=3)) mm |
      
      #### Spatial Extent
      | Axis | Range (mm) |
      |------|------------|
      | X | $(round(dx, digits=2)) |
      | Y | $(round(dy, digits=2)) |
      | Z | $(round(dz, digits=2)) |
      
      #### Start Point (s=0)
      | x | y | z |
      |---|---|---|
      | $(round(p_start.x, digits=2)) | $(round(p_start.y, digits=2)) | $(round(p_start.z, digits=2)) |
      
      #### End Point (s=$(round(track_length, digits=2)))
      | x | y | z |
      |---|---|---|
      | $(round(p_end.x, digits=2)) | $(round(p_end.y, digits=2)) | $(round(p_end.z, digits=2)) |
      """

      return Markdown.parse(md)
  end


  function diffusion_params_md(dfpars)
      md = """
      | Parameter | Value |
      |-----------|-------|
      | Drift length | $(dfpars.ldrift) cm |
      | σ_t (transverse) | $(round(dfpars.sigma_t, digits=3)) mm |
      | σ_l (longitudinal) | $(round(dfpars.sigma_l, digits=3)) mm |
      | Voxel size | $(round(dfpars.voxel_size, digits=3)) mm |
      | Max distance | $(round(dfpars.max_distance, digits=3)) mm |
      | Energy threshold | $(round(dfpars.energy_threshold, digits=3)) keV |
      | Diffusion bins | $(dfpars.nbins_df) |
      | Diffusion nsigma | $(dfpars.nsigma_df) |
      """
      return Markdown.parse(md)
  end


  
"""
      extreme_distances_md(dists)

Return markdown summary of extreme distances for Pluto display.
Takes either the NamedTuple from compute_extreme_distances or individual d1, d2 values.
"""
function extreme_distances_md(dists::NamedTuple)
    if dists.pairing == :none
        return Markdown.parse("**Extreme Distances:** No valid paths")
    end

    pairing_str = dists.pairing == :direct ?
        "reco₁→mc₁, reco₂→mc₂" : "reco₁→mc₂, reco₂→mc₁"

    md = """
    ## Extreme Distances (RECO vs MC)
    
    | Metric | Value |
    |--------|-------|
    | **d₁** | $(round(dists.d1, digits=2)) mm |
    | **d₂** | $(round(dists.d2, digits=2)) mm |
    | **Total** | $(round(dists.total, digits=2)) mm |
    | **Pairing** | $(pairing_str) |
    
    *Optimal pairing minimizes total distance.*
    """

    return Markdown.parse(md)
end

# Convenience method for d1, d2 directly
function extreme_distances_md(d1::Real, d2::Real)
    if isnan(d1) || isnan(d2)
        return Markdown.parse("**Extreme Distances:** Not available")
    end

    md = """
    #### Extreme Distances (RECO vs MC)

    | Metric | Value |
    |--------|-------|
    | **d₁** | $(round(d1, digits=2)) mm |
    | **d₂** | $(round(d2, digits=2)) mm |
    | **Total** | $(round(d1 + d2, digits=2)) mm |
    """

    return Markdown.parse(md)
end


# ============================================================
# Print versions (console output with println)
# ============================================================

"""
    walk_result_print(walk_result)

Print summary of walk_result to console.
"""
function walk_result_print(walk_result)
    if isempty(walk_result.path_indices)
        println("Walk Result: No valid path found")
        return
    end

    ext1, ext2 = walk_result.extremes

    println("═══════════════════════════════════════")
    println("         Walk Result Summary")
    println("═══════════════════════════════════════")
    println("Path length:      $(round(walk_result.total_length, digits=2)) mm")
    println("Confidence:       $(round(walk_result.confidence, digits=3))")
    println("N voxels in path: $(length(walk_result.path_indices))")
    println("───────────────────────────────────────")
    println("Extreme 1 (start):")
    println("  x=$(round(ext1.x, digits=2)), y=$(round(ext1.y, digits=2)), z=$(round(ext1.z, digits=2)), E=$(round(ext1.energy*1e+3, digits=2))")
    println("───────────────────────────────────────")
    println("Extreme 2 (end):")
    println("  x=$(round(ext2.x, digits=2)), y=$(round(ext2.y, digits=2)), z=$(round(ext2.z, digits=2)), E=$(round(ext2.energy*1e+3, digits=2))")
    println("═══════════════════════════════════════")
end


"""
    path_print(path::DataFrame)

Print summary of path DataFrame to console.
"""
function path_print(path::DataFrame)
    if nrow(path) == 0
        println("Path: Empty")
        return
    end

    # Basic stats
    track_length = path.s[end]
    n_points = nrow(path)

    # Spatial extent
    dx = maximum(path.x) - minimum(path.x)
    dy = maximum(path.y) - minimum(path.y)
    dz = maximum(path.z) - minimum(path.z)

    # Start and end positions
    p_start = path[1, :]
    p_end = path[end, :]

    # Average step size
    avg_step = n_points > 1 ? track_length / (n_points - 1) : 0.0

    println("═══════════════════════════════════════")
    println("            Path Summary")
    println("═══════════════════════════════════════")
    println("Track length: $(round(track_length, digits=2)) mm")
    println("N points:     $(n_points)")
    println("Avg step:     $(round(avg_step, digits=3)) mm")
    println("───────────────────────────────────────")
    println("Spatial Extent:")
    println("  ΔX = $(round(dx, digits=2)) mm")
    println("  ΔY = $(round(dy, digits=2)) mm")
    println("  ΔZ = $(round(dz, digits=2)) mm")
    println("───────────────────────────────────────")
    println("Start Point (s=0):")
    println("  x=$(round(p_start.x, digits=2)), y=$(round(p_start.y, digits=2)), z=$(round(p_start.z, digits=2))")
    println("───────────────────────────────────────")
    println("End Point (s=$(round(track_length, digits=2))):")
    println("  x=$(round(p_end.x, digits=2)), y=$(round(p_end.y, digits=2)), z=$(round(p_end.z, digits=2))")
    println("═══════════════════════════════════════")
end


"""
    diffusion_params_print(dfpars)

Print diffusion parameters to console.
"""
function diffusion_params_print(dfpars)
    println("═══════════════════════════════════════")
    println("        Diffusion Parameters")
    println("═══════════════════════════════════════")
    println("Drift length:      $(dfpars.ldrift) cm")
    println("σ_t (transverse):  $(round(dfpars.sigma_t, digits=3)) mm")
    println("σ_l (longitudinal):$(round(dfpars.sigma_l, digits=3)) mm")
    println("Voxel size:        $(round(dfpars.voxel_size, digits=3)) mm")
    println("Max distance:      $(round(dfpars.max_distance, digits=3)) mm")
    println("Energy threshold:  $(round(dfpars.energy_threshold, digits=3)) keV")
    println("Diffusion bins:    $(dfpars.nbins_df)")
    println("Diffusion nsigma:  $(dfpars.nsigma_df)")
    println("═══════════════════════════════════════")
end


"""
    extreme_distances_print(dists::NamedTuple)

Print extreme distances summary to console.
"""
function extreme_distances_print(dists::NamedTuple)
    if dists.pairing == :none
        println("Extreme Distances: No valid paths")
        return
    end

    pairing_str = dists.pairing == :direct ?
        "reco₁→mc₁, reco₂→mc₂" : "reco₁→mc₂, reco₂→mc₁"

    println("═══════════════════════════════════════")
    println("    Extreme Distances (RECO vs MC)")
    println("═══════════════════════════════════════")
    println("d₁:      $(round(dists.d1, digits=2)) mm")
    println("d₂:      $(round(dists.d2, digits=2)) mm")
    println("Total:   $(round(dists.total, digits=2)) mm")
    println("Pairing: $(pairing_str)")
    println("───────────────────────────────────────")
    println("(Optimal pairing minimizes total distance)")
    println("═══════════════════════════════════════")
end

# Convenience method for d1, d2 directly
function extreme_distances_print(d1::Real, d2::Real)
    if isnan(d1) || isnan(d2)
        println("Extreme Distances: Not available")
        return
    end

    println("═══════════════════════════════════════")
    println("    Extreme Distances (RECO vs MC)")
    println("═══════════════════════════════════════")
    println("d₁:    $(round(d1, digits=2)) mm")
    println("d₂:    $(round(d2, digits=2)) mm")
    println("Total: $(round(d1 + d2, digits=2)) mm")
    println("═══════════════════════════════════════")
end


# ============================================================
# KDE Peak Analysis
# ============================================================

"""
    kde_peaks(peaks, kde_f)

Extract peak1 (leftmost) and peak2 (rightmost) from KDE peaks.
Prominence is normalized by the KDE range.

# Arguments
- `peaks`: Result from Petit.find_peaks()
- `kde_f`: KDE values (for normalization)

# Returns
NamedTuple with:
- `peak1_left`, `peak1_right`, `peak1_prom`: Leftmost peak info
- `peak2_left`, `peak2_right`, `peak2_prom`: Rightmost peak info

Values are 0.0 if peak not found.
"""
function kde_peaks(peaks, kde_f)
    # Default values (0 = not found)
    peak1_left, peak1_right, peak1_prom = 0.0, 0.0, 0.0
    peak2_left, peak2_right, peak2_prom = 0.0, 0.0, 0.0

    f_range = maximum(kde_f) - minimum(kde_f)

    if length(peaks.indices) >= 1
        # Find leftmost peak (peak1)
        left_idx = argmin(peaks.positions)
        peak1_left = peaks.lefts[left_idx]
        peak1_right = peaks.rights[left_idx]
        peak1_prom = f_range > 0 ? peaks.proms[left_idx] / f_range : 0.0

        if length(peaks.indices) >= 2
            # Find rightmost peak (peak2)
            right_idx = argmax(peaks.positions)
            peak2_left = peaks.lefts[right_idx]
            peak2_right = peaks.rights[right_idx]
            peak2_prom = f_range > 0 ? peaks.proms[right_idx] / f_range : 0.0
        end
    end

    (peak1_left=peak1_left, peak1_right=peak1_right, peak1_prom=peak1_prom,
     peak2_left=peak2_left, peak2_right=peak2_right, peak2_prom=peak2_prom)
end