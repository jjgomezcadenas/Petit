# Track Path Reconstruction Functions
# Functions for reconstructing and analyzing track paths

using LinearAlgebra
using NearestNeighbors
using DataFrames

"""
    reconstruct_central_path(track::Tracks,
                             path::Vector{Int};
                             filter_radius::Float64)

Given a track and its primary skeleton (as voxel indices `path`),
apply the spatial charge filter of Li & Zeng (2018):

For each skeleton voxel, replace it by the charge-weighted barycentre
of all voxels within `filter_radius`. Returns a DataFrame with the
smoothed central path.

Arguments
---------
track         :: Tracks           # must have track.voxels with x,y,z,energy
path          :: Vector{Int}      # indices into track.voxels, ordered skeleton
filter_radius :: Float64          # R_f, in same units as x,y,z (e.g. mm)

Returns
-------
central_path  :: DataFrame with columns:
    idx_skel   :: Int        # original position index along skeleton (1..length(path))
    x, y, z    :: Float64    # smoothed coordinates
    energy_sum :: Float64    # total energy in the neighbourhood (optional)
"""
function reconstruct_central_path(track::Tracks,
                                  path::Vector{Int};
                                  filter_radius::Float64)

    vox = track.voxels
    Nvox = nrow(vox)
    Ns   = length(path)

    Ns == 0 && return DataFrame(idx_skel = Int[], x = Float64[],
                                y = Float64[], z = Float64[],
                                energy_sum = Float64[])

    # Build coordinate matrix for ALL voxels: 3×Nvox
    coords = Matrix{Float64}(undef, 3, Nvox)
    @inbounds for i in 1:Nvox
        coords[1, i] = vox.x[i]
        coords[2, i] = vox.y[i]
        coords[3, i] = vox.z[i]
    end

    tree = KDTree(coords)

    xs = Vector{Float64}(undef, Ns)
    ys = Vector{Float64}(undef, Ns)
    zs = Vector{Float64}(undef, Ns)
    Es = Vector{Float64}(undef, Ns)

    R = filter_radius

    @inbounds for (k, idx) in enumerate(path)
        # position of skeleton voxel k
        px = coords[1, idx]
        py = coords[2, idx]
        pz = coords[3, idx]

        # Use view into coords matrix (zero allocation)
        q = @view coords[:, idx]

        # all voxels within radius R
        neigh = inrange(tree, q, R, true)

        if isempty(neigh)
            # fallback: keep original voxel if no neighbours (rare)
            xs[k] = px
            ys[k] = py
            zs[k] = pz
            Es[k] = vox.energy[idx]
            continue
        end

        Ex = 0.0
        Ey = 0.0
        Ez = 0.0
        Etot = 0.0

        for j in neigh
            e = vox.energy[j]
            Etot += e
            Ex   += e * coords[1, j]
            Ey   += e * coords[2, j]
            Ez   += e * coords[3, j]
        end

        xs[k] = Ex / Etot
        ys[k] = Ey / Etot
        zs[k] = Ez / Etot
        Es[k] = Etot
    end

    return DataFrame(
        idx_skel   = collect(1:Ns),
        x          = xs,
        y          = ys,
        z          = zs,
        energy_sum = Es,
    )
end



"""
    reconstruct_end_point_energy(central_path::DataFrame; n_points::Int=1)

Given a reconstructed central path (from `reconstruct_central_path`),
compute the blob energies at each endpoint. Blob 1 (Eb1) is the one
with higher energy, Blob 2 (Eb2) has lower energy.

# Arguments
- `central_path::DataFrame`: Output from `reconstruct_central_path` with columns
  idx_skel, x, y, z, energy_sum
- `n_points::Int`: Number of points from each end to sum for blob energy (default: 1)

# Returns
NamedTuple with:
- `blob1`: (x, y, z, energy, idx_start, idx_end) - higher energy blob (Eb1)
- `blob2`: (x, y, z, energy, idx_start, idx_end) - lower energy blob (Eb2)
- `Eb1`: Float64 - energy of blob 1
- `Eb2`: Float64 - energy of blob 2
- `ratio`: Eb1 / Eb2
"""
function reconstruct_end_point_energy(central_path::DataFrame; n_points::Int=1)

    Npath = nrow(central_path)
    Npath ≥ 2 || error("central_path must have at least 2 points")

    # Clamp n_points to avoid overlap
    n = min(n_points, Npath ÷ 2)
    n = max(1, n)

    # Sum energy from first n points (start endpoint)
    E_start = sum(central_path.energy_sum[1:n])
    x_start = central_path.x[1]
    y_start = central_path.y[1]
    z_start = central_path.z[1]

    # Sum energy from last n points (end endpoint)
    E_end = sum(central_path.energy_sum[(Npath-n+1):Npath])
    x_end = central_path.x[end]
    y_end = central_path.y[end]
    z_end = central_path.z[end]

    # Blob1 = higher energy, Blob2 = lower energy
    if E_start >= E_end
        blob1 = (x = x_start, y = y_start, z = z_start,
                 energy = E_start,
                 idx_start = 1, idx_end = n)
        blob2 = (x = x_end, y = y_end, z = z_end,
                 energy = E_end,
                 idx_start = Npath - n + 1, idx_end = Npath)
    else
        blob1 = (x = x_end, y = y_end, z = z_end,
                 energy = E_end,
                 idx_start = Npath - n + 1, idx_end = Npath)
        blob2 = (x = x_start, y = y_start, z = z_start,
                 energy = E_start,
                 idx_start = 1, idx_end = n)
    end

    Eb1 = blob1.energy
    Eb2 = blob2.energy
    ratio = Eb1 / max(Eb2, eps())

    return (blob1 = blob1,
            blob2 = blob2,
            Eb1 = Eb1,
            Eb2 = Eb2,
            asymmetry = (Eb1-Eb2)/(Eb1+Eb2))
end




"""
    get_raw_path(track, path_indices)

Extract the raw (unsmoothed) path from track voxels using path indices.
Returns a DataFrame with x, y, z, s (arc-length).
"""
function get_raw_path(track, path_indices::Vector{Int})
    voxels = track.voxels
    n = length(path_indices)

    xs = voxels.x[path_indices]
    ys = voxels.y[path_indices]
    zs = voxels.z[path_indices]

    # Compute cumulative arc-length
    ss = Vector{Float64}(undef, n)
    ss[1] = 0.0
    for i in 2:n
        dx = xs[i] - xs[i-1]
        dy = ys[i] - ys[i-1]
        dz = zs[i] - zs[i-1]
        ss[i] = ss[i-1] + sqrt(dx^2 + dy^2 + dz^2)
    end

    return DataFrame(x=xs, y=ys, z=zs, s=ss)
end

"""
    project_voxels_to_path(voxels, path)

Project all voxels onto the path, returning arc-length for each voxel.
Uses nearest-neighbor projection.
"""
function project_voxels_to_path(voxels::DataFrame, path::DataFrame)
    n_voxels = nrow(voxels)
    s_voxels = Vector{Float64}(undef, n_voxels)

    for i in 1:n_voxels
        vx, vy, vz = voxels.x[i], voxels.y[i], voxels.z[i]

        # Find closest path point
        min_dist = Inf
        closest_idx = 1
        for j in 1:nrow(path)
            dx = vx - path.x[j]
            dy = vy - path.y[j]
            dz = vz - path.z[j]
            dist = sqrt(dx^2 + dy^2 + dz^2)
            if dist < min_dist
                min_dist = dist
                closest_idx = j
            end
        end

        s_voxels[i] = path.s[closest_idx]
    end

    return s_voxels
end

"""
    compute_extreme_distances(reco_path::DataFrame, mc_path::DataFrame)

Compute euclidean distances between reco and MC path extremes.

Finds the optimal pairing (minimizing total distance) between:
- reco_path extremes: [1,:] and [end,:]
- mc_path extremes: [1,:] and [end,:]

Returns a NamedTuple with:
- `d1`: distance from reco extreme 1 to its matched MC extreme
- `d2`: distance from reco extreme 2 to its matched MC extreme
- `total`: d1 + d2
- `pairing`: :direct (reco1→mc1, reco2→mc2) or :crossed (reco1→mc2, reco2→mc1)
"""
function compute_extreme_distances(reco_path::DataFrame, mc_path::DataFrame)
    # Handle empty paths
    if nrow(reco_path) == 0 || nrow(mc_path) == 0
        return (d1=NaN, d2=NaN, total=NaN, pairing=:none)
    end

    # Reco extremes
    rx1, ry1, rz1 = reco_path.x[1], reco_path.y[1], reco_path.z[1]
    rx2, ry2, rz2 = reco_path.x[end], reco_path.y[end], reco_path.z[end]

    # MC extremes
    mx1, my1, mz1 = mc_path.x[1], mc_path.y[1], mc_path.z[1]
    mx2, my2, mz2 = mc_path.x[end], mc_path.y[end], mc_path.z[end]

    # Direct pairing: reco1→mc1, reco2→mc2
    d_r1_m1 = sqrt((rx1 - mx1)^2 + (ry1 - my1)^2 + (rz1 - mz1)^2)
    d_r2_m2 = sqrt((rx2 - mx2)^2 + (ry2 - my2)^2 + (rz2 - mz2)^2)
    direct_total = d_r1_m1 + d_r2_m2

    # Crossed pairing: reco1→mc2, reco2→mc1
    d_r1_m2 = sqrt((rx1 - mx2)^2 + (ry1 - my2)^2 + (rz1 - mz2)^2)
    d_r2_m1 = sqrt((rx2 - mx1)^2 + (ry2 - my1)^2 + (rz2 - mz1)^2)
    crossed_total = d_r1_m2 + d_r2_m1

    # Choose optimal pairing
    if direct_total <= crossed_total
        return (d1=d_r1_m1, d2=d_r2_m2, total=direct_total, pairing=:direct)
    else
        return (d1=d_r1_m2, d2=d_r2_m1, total=crossed_total, pairing=:crossed)
    end
end
