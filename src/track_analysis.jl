# Track Analysis Utilities
# Utility functions and diagnostics for track analysis

using DataFrames
using Graphs
using NearestNeighbors



"""
    track_positions(tracks::Vector{Tracks})

Extract all voxel positions from a vector of tracks.

# Returns
- `(X, Y, Z)`: Vectors of all x, y, z coordinates
"""
function track_positions(tracks::Vector{Tracks})
	X = Float64[]
	Y = Float64[]
	Z = Float64[]
	for i in 1:length(tracks)
		append!(X, tracks[i].voxels.x)
		append!(Y, tracks[i].voxels.y)
		append!(Z, tracks[i].voxels.z)
	end
	X, Y, Z
end


"""
    track_energies_keV(tracks::Vector{Tracks})

Compute total energy (keV) for each track in the collection.

# Returns
- `Vector{Float64}`: Energy in keV for each track
"""
function track_energies_keV(tracks::Vector{Tracks})
	E = Float64[]
	for i in 1:length(tracks)
		energy_kev = 1e+3 * sum(tracks[i].voxels.energy)
		push!(E, energy_kev)
	end
	E
end


"""
    track_stats(tracks::Vector{Tracks})

Compute summary statistics for a collection of tracks.

# Returns
NamedTuple with:
- `n_voxels::Vector{Int}`: Number of voxels in each track
- `energies_keV::Vector{Float64}`: Total energy (keV) for each track
"""
function track_stats(tracks::Vector{Tracks})
    n_voxels = [nrow(t.voxels) for t in tracks]
    energies_keV = [1e+3 * sum(t.voxels.energy) for t in tracks]
    (n_voxels=n_voxels, energies_keV=energies_keV)
end


"""
    walk_result_print(walk_result; track=nothing, method="")

Print summary of walk_result to console.

# Arguments
- `walk_result`: Result from walk_track_from_extremes
- `track`: Optional Tracks object to include graph info (edges, vertices)
- `method`: Graph construction method string (e.g., "KDT", "kNN", "kNN_mutual")
"""
function walk_result_print(walk_result; track=nothing, method::String="")
    if isempty(walk_result.path_indices)
        println("Walk Result: No valid path found")
        return
    end

    ext1, ext2 = walk_result.extremes
    title = isempty(method) ? "Walk Result Summary" : "Walk Result Summary ($method)"

    println("═══════════════════════════════════════")
    println("         $title")
    println("═══════════════════════════════════════")
    if !isnothing(track)
        method_str = isempty(method) ? "" : " ($method)"
        println("Graph$method_str:       $(ne(track.graph)) edges, $(nv(track.graph)) vertices")
    end
    println("Path length:      $(round(walk_result.total_length, digits=2)) mm")
    println("Confidence:       $(round(walk_result.confidence, digits=3))")
    println("N voxels in path: $(length(walk_result.path_indices))")
    println("───────────────────────────────────────")
    println("Extreme 1 (start):")
    println("  x=$(round(ext1.x, digits=2)), y=$(round(ext1.y, digits=2)), z=$(round(ext1.z, digits=2)), E=$(round(ext1.energy*1e+3, digits=2)) keV")
    println("───────────────────────────────────────")
    println("Extreme 2 (end):")
    println("  x=$(round(ext2.x, digits=2)), y=$(round(ext2.y, digits=2)), z=$(round(ext2.z, digits=2)), E=$(round(ext2.energy*1e+3, digits=2)) keV")
    println("═══════════════════════════════════════")
end


"""
    extreme_distances_print(dists::NamedTuple; method="")

Print extreme distances summary to console.

# Arguments
- `dists`: NamedTuple from compute_extreme_distances with d1, d2, total, pairing
- `method`: Graph construction method string (e.g., "KDT", "kNN", "kNN_mutual")
"""
function extreme_distances_print(dists::NamedTuple; method::String="")
    if dists.pairing == :none
        println("Extreme Distances: No valid paths")
        return
    end

    pairing_str = dists.pairing == :direct ?
        "reco₁→mc₁, reco₂→mc₂" : "reco₁→mc₂, reco₂→mc₁"
    title = isempty(method) ? "Extreme Distances (RECO vs MC)" : "Extreme Distances ($method)"

    println("═══════════════════════════════════════")
    println("    $title")
    println("═══════════════════════════════════════")
    println("d₁:      $(round(dists.d1, digits=2)) mm")
    println("d₂:      $(round(dists.d2, digits=2)) mm")
    println("Total:   $(round(dists.total, digits=2)) mm")
    println("Pairing: $(pairing_str)")
    println("═══════════════════════════════════════")
end

# Convenience method for d1, d2 directly
function extreme_distances_print(d1::Real, d2::Real; method::String="")
    if isnan(d1) || isnan(d2)
        println("Extreme Distances: Not available")
        return
    end

    title = isempty(method) ? "Extreme Distances (RECO vs MC)" : "Extreme Distances ($method)"

    println("═══════════════════════════════════════")
    println("    $title")
    println("═══════════════════════════════════════")
    println("d₁:    $(round(d1, digits=2)) mm")
    println("d₂:    $(round(d2, digits=2)) mm")
    println("Total: $(round(d1 + d2, digits=2)) mm")
    println("═══════════════════════════════════════")
end


"""
    path_print(path::DataFrame; method="")

Print summary of path DataFrame to console.

# Arguments
- `path`: DataFrame from get_raw_path with x, y, z, s columns
- `method`: Graph construction method string (e.g., "KDT", "kNN", "kNN_mutual")
"""
function path_print(path::DataFrame; method::String="")
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

    title = isempty(method) ? "Path Summary" : "Path Summary ($method)"

    println("═══════════════════════════════════════")
    println("            $title")
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
    blobs_print(blobs; method="")

Print blob analysis results to console.

# Arguments
- `blobs`: NamedTuple from find_blob_energies with Eb1, Eb2, asymmetry, blob1, blob2
- `method`: Graph construction method string (e.g., "KDT", "kNN", "kNN_mutual")

Note: blobs.Eb1 and blobs.Eb2 are already in keV.
"""
function blobs_print(blobs; method::String="")
    title = isempty(method) ? "Blob Analysis" : "Blob Analysis ($method)"

    println("═══════════════════════════════════════")
    println("          $title")
    println("═══════════════════════════════════════")
    println("Blob1 (high energy):")
    println("  position = ($(round(blobs.blob1.x, digits=2)), $(round(blobs.blob1.y, digits=2)), $(round(blobs.blob1.z, digits=2))) mm")
    println("  energy   = $(round(blobs.Eb1, digits=1)) keV")
    println("───────────────────────────────────────")
    println("Blob2 (low energy):")
    println("  position = ($(round(blobs.blob2.x, digits=2)), $(round(blobs.blob2.y, digits=2)), $(round(blobs.blob2.z, digits=2))) mm")
    println("  energy   = $(round(blobs.Eb2, digits=1)) keV")
    println("───────────────────────────────────────")
    println("Asymmetry  = $(round(blobs.asymmetry, digits=3))")
    println("═══════════════════════════════════════")
end


"""
    diffusion_params_print(dfpars)

Print diffusion parameters to console.

# Arguments
- `dfpars`: DiffusionParams struct with ldrift, sigma_t, sigma_l, voxel_size, max_distance, energy_threshold, nbins_df, nsigma_df
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


#####################
# Diagnosis utilities #
#####################

"""
    diagnose_path_efficiency(coords, path) -> (η, L_path, D_end)

η = L_path / D_end, where:
- L_path is Euclidean length along the vertex sequence `path`
- D_end is straight-line distance between endpoints
"""
function diagnose_path_efficiency(coords::TrackCoords, path::Vector{Int})
    isempty(path) && return (Inf, 0.0, 0.0)
    coord_matrix = hcat(coords.x, coords.y, coords.z)'
    L = calculate_path_length_from_coords(coord_matrix, path)
    u = first(path); v = last(path)
    D = euclidean_distance(coords.x[u], coords.y[u], coords.z[u],
                           coords.x[v], coords.y[v], coords.z[v])
    η = (D > 0) ? (L / D) : Inf
    return (η, L, D)
end


"""
    diagnose_endpoint_degrees(g, u, v) -> (deg_u, deg_v)
"""
diagnose_endpoint_degrees(g::SimpleGraph, u::Int, v::Int) = (degree(g, u), degree(g, v))


"""
    diagnose_skeleton_coverage(track, coords, path; R_cover, energy_col=:energy) -> (f, E_in, E_tot)

Fraction of total energy within distance `R_cover` of the skeleton path.
Uses a KDTree over skeleton vertices and nearest-neighbor queries.
"""
function diagnose_skeleton_coverage(track::Tracks, coords::TrackCoords, path::Vector{Int};
                                   R_cover::Float64,
                                   energy_col::Symbol = :energy)

    vox = track.voxels
    n = nrow(vox)
    (n == 0 || isempty(path)) && return (0.0, 0.0, 0.0)

    # Total energy
    Evec = Float64.(getproperty(vox, energy_col))
    Etot = sum(Evec)
    Etot == 0 && return (0.0, 0.0, 0.0)

    # Build KDTree over skeleton points, filtering out invalid indices and NaN coordinates
    valid_path = Int[]
    for idx in path
        if idx >= 1 && idx <= length(coords.x)
            x, y, z = coords.x[idx], coords.y[idx], coords.z[idx]
            if !isnan(x) && !isnan(y) && !isnan(z)
                push!(valid_path, idx)
            end
        end
    end

    isempty(valid_path) && return (0.0, 0.0, Etot)

    Ns = length(valid_path)
    skel = Matrix{Float64}(undef, 3, Ns)
    @inbounds for (k, idx) in enumerate(valid_path)
        skel[1,k] = coords.x[idx]
        skel[2,k] = coords.y[idx]
        skel[3,k] = coords.z[idx]
    end
    tree = KDTree(skel)

    Ein = 0.0
    @inbounds for i in 1:n
        x, y, z = coords.x[i], coords.y[i], coords.z[i]
        # Skip points with NaN coordinates
        (isnan(x) || isnan(y) || isnan(z)) && continue
        idxs, dists = knn(tree, [x, y, z], 1, true)
        if !isempty(dists) && dists[1] <= R_cover
            Ein += Evec[i]
        end
    end

    f = Ein / Etot
    return (f, Ein, Etot)
end


"""
    diagnose_endpoint_stability(f1, f2, coords; match_by=:minsum) -> Δ

Compute endpoint displacement between two solutions (u1,v1) and (u2,v2) given the same `coords`.
Allows swapping endpoints to minimize total displacement.
- `f1`, `f2` are tuples with endpoints as their first two elements (u,v,...).
"""
function diagnose_endpoint_stability(res1, res2, coords::TrackCoords)
    u1, v1 = res1[1], res1[2]
    u2, v2 = res2[1], res2[2]

    (u1 === nothing || v1 === nothing || u2 === nothing || v2 === nothing) && return Inf

    function dist(a,b)
        euclidean_distance(coords.x[a], coords.y[a], coords.z[a],
                           coords.x[b], coords.y[b], coords.z[b])
    end

    d_noswap = max(dist(u1,u2), dist(v1,v2))
    d_swap   = max(dist(u1,v2), dist(v1,u2))
    return min(d_noswap, d_swap)
end


"""
    diagnose_track_extent(coords, vertices=1:length(coords.x)) -> extent_mm

Compute a crude spatial extent (diagonal of bounding box) for selected vertices.
Useful for conditioning shortcut alarms (η near 1 is only suspicious if extent is large).
"""
function diagnose_track_extent(coords::TrackCoords, vertices::AbstractVector{Int}=collect(1:length(coords.x)))
    isempty(vertices) && return 0.0
    xs = coords.x[vertices]; ys = coords.y[vertices]; zs = coords.z[vertices]
    dx = maximum(xs) - minimum(xs)
    dy = maximum(ys) - minimum(ys)
    dz = maximum(zs) - minimum(zs)
    return sqrt(dx*dx + dy*dy + dz*dz)
end


"""
    kde_peaks_print(pk; method="")

Print KDE peaks summary to console.
"""
function kde_peaks_print(pk; method::String="")
    if pk.peak1_prom == 0.0 && pk.peak2_prom == 0.0
        println("KDE Peaks: None found")
        return
    end

    r = x -> round(x, digits=1)
    p1 = pk.peak1_prom > 0 ? "P1=[$(r(pk.peak1_left)),$(r(pk.peak1_right))] prom=$(round(pk.peak1_prom, digits=2))" : "P1=none"
    p2 = pk.peak2_prom > 0 ? "P2=[$(r(pk.peak2_left)),$(r(pk.peak2_right))] prom=$(round(pk.peak2_prom, digits=2))" : "P2=none"

    title = isempty(method) ? "KDE Peaks" : "KDE Peaks ($method)"
    println("$title: $p1 | $p2")
end


# ============================================================
# Markdown display functions (for Pluto notebooks)
# ============================================================

using Markdown

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


"""
    diffusion_params_md(dfpars)

Return markdown summary of diffusion parameters for Pluto display.
"""
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


"""
    kde_peaks_md(pk)

Return markdown summary of KDE peaks for Pluto display.
"""
function kde_peaks_md(pk)
    r = x -> round(x, digits=1)

    if pk.peak1_prom == 0.0 && pk.peak2_prom == 0.0
        return Markdown.parse("**KDE Peaks:** None found")
    end

    p1 = pk.peak1_prom > 0 ? "P1=[$(r(pk.peak1_left)),$(r(pk.peak1_right))] prom=$(round(pk.peak1_prom, digits=2))" : "P1=none"
    p2 = pk.peak2_prom > 0 ? "P2=[$(r(pk.peak2_left)),$(r(pk.peak2_right))] prom=$(round(pk.peak2_prom, digits=2))" : "P2=none"

    return Markdown.parse("**KDE Peaks:** $p1 | $p2")
end


"""
    blobs_md(blobs)

Return markdown summary of blob analysis for Pluto display.
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

