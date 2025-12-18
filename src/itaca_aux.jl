# Note: This file is included after Petit module is loaded
# Functions use Petit.* via the including script's scope
using DataFrames
using Markdown
# Petit module is already defined via include() in the parent script
# We use .Petit to reference it (local module, not a package)

function get_sigma(particle_type, ldrft; 
				   dt = 3.5, dl = 0.9, 
				   tK = 297.0, edrift = 500.0, Pbar=15.0)
	
	if particle_type == "ion"
		σt =  Petit.sigma_t_ion_mm(tK, ldrft, edrift)
		σl =  0.0
	else
		σt =  Petit.sigma_t_mm(ldrft, Pbar; dtmm=dt)
		σl =  Petit.sigma_l_mm(ldrft, Pbar; dlmm=dl)
	end
	σt, σl
end


function get_energy_threshold(particle_type; 
							  energy_threshold_ions =  10.0,
							  energy_threshold_keV =  10.0)
	f = 1e+5/2.5 # ions per MeV
	fkeV = f*1e-3 # ions per keV

	if particle_type == "ion"
		energy_threshold_keV = energy_threshold_ions/fkeV
	end
	energy_threshold_keV
end

function get_voxel_size_and_distance(ldrft, σt)
	if ldrft > 50.0
		voxel_scale = 1.5 
		voxel_distance_scale = 1.5
	else
		voxel_scale = 3.0 
		voxel_distance_scale = 2.0
	end
	
	voxel_size = σt * voxel_scale
	mcvox_size = 0.5
	max_distance = voxel_size * voxel_distance_scale
	(voxel_size, mcvox_size, max_distance)
end


function length_and_energy_of_tracks(tracks)
	LT = [length(track.voxels.energy) for track in tracks]
	E = [sum(track.voxels.energy)*1e+3 for track in tracks]
	LT, E
end

function print_diagnostics(label, res, track, coords; R_cover=3.0)
    u, v, path = res[1], res[2], res[3]
    println("\n--- $label ---")
    println("endpoints = ($u, $v)")

    # A) path efficiency
    η, L_path, D_end = Petit.diagnose_path_efficiency(coords, path)
    println("path_efficiency η = L_path/D_end = $(round(η, digits=3))  (L_path=$(round(L_path,digits=2)) mm, D_end=$(round(D_end,digits=2)) mm)")

    # B) endpoint degrees (interior-ness)
    deg_u, deg_v = Petit.diagnose_endpoint_degrees(track.graph, u, v)
    println("deg(u)=$deg_u, deg(v)=$deg_v")

    # C) skeleton coverage (how much energy lies near skeleton)
    f, Ein, Etot = Petit.diagnose_skeleton_coverage(track, coords, path; R_cover=R_cover, energy_col=:energy)
    println("coverage f = Ein/Etot = $(round(f,digits=3))  (Ein=$(round(Ein,digits=3)), Etot=$(round(Etot,digits=3)))")
end

"""
    walk_result_md(walk_result)

Return markdown summary of walk_result for Pluto display.
"""
function walk_result_md(walk_result)
    if isempty(walk_result.path_indices)
        return Markdown.parse("**Walk:** No valid path")
    end

    ext1, ext2 = walk_result.extremes
    r = x -> round(x, digits=1)

    md = """
    **Walk:** L=$(r(walk_result.total_length))mm, N=$(length(walk_result.path_indices)), conf=$(round(walk_result.confidence, digits=2)) |
    **Ext1** ($(r(ext1.x)),$(r(ext1.y)),$(r(ext1.z))) E=$(r(ext1.energy*1e+3))keV |
    **Ext2** ($(r(ext2.x)),$(r(ext2.y)),$(r(ext2.z))) E=$(r(ext2.energy*1e+3))keV
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

    r = x -> round(x, digits=1)
    L = path.s[end]
    n = nrow(path)
    dx = maximum(path.x) - minimum(path.x)
    dy = maximum(path.y) - minimum(path.y)
    dz = maximum(path.z) - minimum(path.z)
    p1, p2 = path[1, :], path[end, :]

    md = """
    **Path:** L=$(r(L))mm, N=$(n) | Δ=($(r(dx)),$(r(dy)),$(r(dz)))mm |
    Start ($(r(p1.x)),$(r(p1.y)),$(r(p1.z))) → End ($(r(p2.x)),$(r(p2.y)),$(r(p2.z)))
    """

    return Markdown.parse(md)
end


function diffusion_params_md(dfpars)
    r = x -> round(x, digits=2)
    md = """
    **Diffusion:** L=$(dfpars.ldrift)cm, σt=$(r(dfpars.sigma_t))mm, σl=$(r(dfpars.sigma_l))mm |
    voxel=$(r(dfpars.voxel_size))mm, dist=$(r(dfpars.max_distance))mm, Eth=$(r(dfpars.energy_threshold))keV
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
        return Markdown.parse("**Distances:** No valid paths")
    end

    r = x -> round(x, digits=1)
    pairing = dists.pairing == :direct ? "direct" : "crossed"

    md = """
    **Distances:** d₁=$(r(dists.d1)), d₂=$(r(dists.d2)), total=$(r(dists.total))mm | $(pairing)
    """

    return Markdown.parse(md)
end

# Convenience method for d1, d2 directly
function extreme_distances_md(d1::Real, d2::Real)
    if isnan(d1) || isnan(d2)
        return Markdown.parse("**Distances:** Not available")
    end

    r = x -> round(x, digits=1)

    md = """
    **Distances:** d₁=$(r(d1)), d₂=$(r(d2)), total=$(r(d1 + d2))mm
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


"""
    kde_peaks_md(pk)

Return compact markdown summary of KDE peaks for Pluto display.
"""
function kde_peaks_md(pk)
    r = x -> round(x, digits=1)

    if pk.peak1_prom == 0.0 && pk.peak2_prom == 0.0
        return Markdown.parse("**KDE Peaks:** None found")
    end

    p1 = pk.peak1_prom > 0 ? "P1=[$(r(pk.peak1_left)),$(r(pk.peak1_right))] prom=$(round(pk.peak1_prom, digits=2))" : "P1=none"
    p2 = pk.peak2_prom > 0 ? "P2=[$(r(pk.peak2_left)),$(r(pk.peak2_right))] prom=$(round(pk.peak2_prom, digits=2))" : "P2=none"

    md = """
    **KDE Peaks:** $(p1) | $(p2)
    """

    return Markdown.parse(md)
end


"""
    blobs_md(blobs)

Return compact markdown summary of blob analysis for Pluto display.
"""
function blobs_md(blobs)
    r = x -> round(x, digits=1)

    md = """
    **Blobs:** Eb1=$(r(blobs.Eb1))keV, Eb2=$(r(blobs.Eb2))keV, asym=$(round(blobs.asymmetry, digits=2)) |
    B1 ($(r(blobs.blob1.x)),$(r(blobs.blob1.y)),$(r(blobs.blob1.z))) n=$(length(blobs.blob1.voxels_idx)) |
    B2 ($(r(blobs.blob2.x)),$(r(blobs.blob2.y)),$(r(blobs.blob2.z))) n=$(length(blobs.blob2.voxels_idx))
    """

    return Markdown.parse(md)
end