# Track Blob Analysis Functions
# Functions for analyzing energy blobs at track endpoints

using DataFrames
using NearestNeighbors
using Statistics

"""
    find_blob_energies(track::Tracks,
                       central_path::DataFrame;
                       radius::Float64)

Given a reconstructed central path and the full set of voxels, sum the
energy in a sphere of radius `radius` around each endpoint of the
central path. Blob1 has higher energy, Blob2 has lower energy.

# Arguments
- `track::Tracks`: Must have `track.voxels::DataFrame` with columns `x, y, z, energy`
- `central_path::DataFrame`: Smoothed path from `reconstruct_central_path` with x, y, z columns
- `radius::Float64`: Radius of the spherical neighbourhood around each endpoint (mm)

# Returns
NamedTuple with:
- `blob1`: (x, y, z, energy, voxels_idx, endpoint_index) - higher energy blob (Eb1)
- `blob2`: (x, y, z, energy, voxels_idx, endpoint_index) - lower energy blob (Eb2)
- `Eb1`: Float64 - energy of blob 1 (higher)
- `Eb2`: Float64 - energy of blob 2 (lower)
- `asymmetry`: (Eb1 - Eb2) / (Eb1 + Eb2)
"""
function find_blob_energies(track::Tracks,
                            central_path::DataFrame;
                            radius::Float64)

    vox  = track.voxels
    Nvox = nrow(vox)
    Nvox == 0 && error("Track has no voxels")

    Npath = nrow(central_path)
    Npath ≥ 2 || error("central_path must have at least 2 points")

    # Build KDTree on all voxels
    coords = Matrix{Float64}(undef, 3, Nvox)
    @inbounds for i in 1:Nvox
        coords[1, i] = vox.x[i]
        coords[2, i] = vox.y[i]
        coords[3, i] = vox.z[i]
    end
    tree = KDTree(coords)

    # Helper: sum energy in sphere
    function energy_in_sphere(xc, yc, zc)
        center = [xc, yc, zc]
        idxs = inrange(tree, center, radius, true)

        if isempty(idxs)
            # Fallback: use nearest voxel if none within radius
            idxs, _ = knn(tree, center, 1, true)
        end

        Etot = 0.0
        @inbounds for j in idxs
            Etot += vox.energy[j]
        end
        return Etot, idxs
    end

    # Smoothed endpoints
    x_start = central_path.x[1]
    y_start = central_path.y[1]
    z_start = central_path.z[1]

    x_end = central_path.x[end]
    y_end = central_path.y[end]
    z_end = central_path.z[end]

    E_start, idxs_start = energy_in_sphere(x_start, y_start, z_start)
    E_end, idxs_end = energy_in_sphere(x_end, y_end, z_end)

    # Blob1 = higher energy, Blob2 = lower energy
    if E_start >= E_end
        blob1 = (x = x_start, y = y_start, z = z_start,
                 energy = E_start,
                 voxels_idx = idxs_start,
                 endpoint_index = 1)
        blob2 = (x = x_end, y = y_end, z = z_end,
                 energy = E_end,
                 voxels_idx = idxs_end,
                 endpoint_index = Npath)
    else
        blob1 = (x = x_end, y = y_end, z = z_end,
                 energy = E_end,
                 voxels_idx = idxs_end,
                 endpoint_index = Npath)
        blob2 = (x = x_start, y = y_start, z = z_start,
                 energy = E_start,
                 voxels_idx = idxs_start,
                 endpoint_index = 1)
    end

    Eb1 = blob1.energy * 1e+3
    Eb2 = blob2.energy * 1e+3

    return (blob1 = blob1,
            blob2 = blob2,
            Eb1 = Eb1,
            Eb2 = Eb2,
            asymmetry = (Eb1 - Eb2) / (Eb1 + Eb2))
end

"""
    energy_in_spheres_around_extremes(track, walk_result, radius)

Calculate energy within spheres of given radius around track endpoints.

# Arguments
- `track::Tracks`: Track object
- `walk_result`: Result from `walk_track_from_extremes`
- `radius::Float64`: Sphere radius in mm

# Returns
NamedTuple with blob1 (higher energy) and blob2 (lower energy) results:
energy, voxel_count, center coordinates for each blob.
"""
function energy_in_spheres_around_extremes(track::Tracks, walk_result, radius::Float64)

    # Check if we have valid extremes
    if isnothing(walk_result.extremes[1])
        return (start_sphere_energy = 0.0,
                end_sphere_energy = 0.0,
                start_voxel_count = 0,
                end_voxel_count = 0,
                start_center = (NaN, NaN, NaN),
                end_center = (NaN, NaN, NaN))
    end

    start_voxel, end_voxel = walk_result.extremes

    # Get coordinates of extreme voxels
    start_center = (start_voxel.x, start_voxel.y, start_voxel.z)
    end_center = (end_voxel.x, end_voxel.y, end_voxel.z)

    # Initialize counters
    start_sphere_energy = 0.0
    end_sphere_energy = 0.0
    start_voxel_count = 0
    end_voxel_count = 0

    # Check each voxel in the track
    for i in 1:nrow(track.voxels)
        voxel_x = track.voxels.x[i]
        voxel_y = track.voxels.y[i]
        voxel_z = track.voxels.z[i]
        voxel_energy = track.voxels.energy[i]

        # Calculate distance to start voxel
        dist_to_start = euclidean_distance(voxel_x, voxel_y, voxel_z,
                                          start_center[1], start_center[2], start_center[3])

        # Calculate distance to end voxel
        dist_to_end = euclidean_distance(voxel_x, voxel_y, voxel_z,
                                        end_center[1], end_center[2], end_center[3])

        # Check if voxel is within start sphere
        if dist_to_start <= radius
            start_sphere_energy += voxel_energy
            start_voxel_count += 1
        end

        # Check if voxel is within end sphere
        if dist_to_end <= radius
            end_sphere_energy += voxel_energy
            end_voxel_count += 1
        end
    end

    # Determine which sphere has more energy
    # blob1 is the sphere with larger energy
    if start_sphere_energy >= end_sphere_energy
        blob1_energy = start_sphere_energy
        blob1_voxel_count = start_voxel_count
        blob1_center = start_center
        blob2_energy = end_sphere_energy
        blob2_voxel_count = end_voxel_count
        blob2_center = end_center
    else
        blob1_energy = end_sphere_energy
        blob1_voxel_count = end_voxel_count
        blob1_center = end_center
        blob2_energy = start_sphere_energy
        blob2_voxel_count = start_voxel_count
        blob2_center = start_center
    end

    # Return both old format (for compatibility) and new format
    return (start_sphere_energy = start_sphere_energy,
            end_sphere_energy = end_sphere_energy,
            start_voxel_count = start_voxel_count,
            end_voxel_count = end_voxel_count,
            start_center = start_center,
            end_center = end_center,
            # New format: blob1 has higher energy than blob2
            blob1_energy = blob1_energy,
            blob1_voxel_count = blob1_voxel_count,
            blob1_center = blob1_center,
            blob2_energy = blob2_energy,
            blob2_voxel_count = blob2_voxel_count,
            blob2_center = blob2_center)
end


"""
    energy_blobs_from_path(walk_result, n::Int; energy_col::Symbol=:energy)

Calculate blob energies using n voxels from each end of the walk path.

This function uses the ordered path_voxels from walk_result directly,
taking n voxels from the start and n voxels from the end of the path.

# Arguments
- `walk_result`: Result from walk_track_from_extremes containing path_voxels
- `n::Int`: Number of voxels to include from each end (n=1 means just the extreme voxel)
- `energy_col::Symbol`: Column to use for energy values (:energy or :electrons, default: :energy)

# Returns
NamedTuple with:
- `eb1::Float64`: Energy of blob with higher energy
- `eb2::Float64`: Energy of blob with lower energy
- `n1::Int`: Number of voxels in blob1
- `n2::Int`: Number of voxels in blob2
- `start_energy::Float64`: Energy at start of path
- `end_energy::Float64`: Energy at end of path
"""
function energy_blobs_from_path(walk_result, n::Int; energy_col::Symbol=:energy)
    # Check if we have valid path_voxels
    if isnothing(walk_result.path_voxels) || nrow(walk_result.path_voxels) == 0
        return (eb1 = 0.0, eb2 = 0.0, n1 = 0, n2 = 0,
                start_energy = 0.0, end_energy = 0.0)
    end

    path_df = walk_result.path_voxels
    nvoxels = nrow(path_df)

    # Check if requested column exists
    if !hasproperty(path_df, energy_col)
        error("Column :$energy_col not found in path_voxels. Available: $(names(path_df))")
    end

    # Get energy values from specified column
    energy_values = Float64.(path_df[!, energy_col])

    # Ensure n is at least 1 and doesn't exceed half the path
    n = max(1, n)
    n_start = min(n, nvoxels ÷ 2)  # Don't overlap with end
    n_end = min(n, nvoxels - n_start)  # Take remaining from end

    # Sum energy from first n_start voxels (start of path)
    start_energy = sum(energy_values[1:n_start])

    # Sum energy from last n_end voxels (end of path)
    end_energy = sum(energy_values[(nvoxels - n_end + 1):nvoxels])

    # Assign blob1 to higher energy, blob2 to lower
    if start_energy >= end_energy
        eb1 = start_energy
        eb2 = end_energy
        n1 = n_start
        n2 = n_end
    else
        eb1 = end_energy
        eb2 = start_energy
        n1 = n_end
        n2 = n_start
    end

    return (eb1 = eb1, eb2 = eb2, n1 = n1, n2 = n2,
            start_energy = start_energy, end_energy = end_energy)
end


"""
    blob_asymmetry(tracks::Vector{Tracks}; i0::Int=1, il::Int=10, r0::Int=5, rl::Int=15)

Compute blob energy asymmetry |eb1-eb2|/(eb1+eb2) over a range of events and sphere radii.

# Arguments
- `tracks::Vector{Tracks}`: Vector of track objects
- `i0::Int`: First event index (default: 1)
- `il::Int`: Last event index (default: 10)
- `r0::Int`: Minimum sphere radius in mm (default: 5)
- `rl::Int`: Maximum sphere radius in mm (default: 15)

# Returns
- `(mean_asymmetry, std_asymmetry)`: Vectors of mean and std of asymmetry for each radius
"""
function blob_asymmetry(tracks::Vector{Tracks}; i0::Int=1, il::Int=10, r0::Int=5, rl::Int=15)
    # Validate indices
    il = min(il, length(tracks))
    i0 = max(1, i0)

    if i0 > il
        error("i0 ($i0) must be <= il ($il)")
    end

    nevents = il - i0 + 1
    nradii = rl - r0 + 1
    DB = Matrix{Float64}(undef, nevents, nradii)

    for (ievent_idx, ievent) in enumerate(i0:il)
        trk = tracks[ievent]
        walk_result = walk_track_from_extremes(trk)

        for (ir, rr) in enumerate(r0:rl)
            sphere_energies = energy_in_spheres_around_extremes(trk, walk_result, Float64(rr))
            eb1 = sphere_energies.blob1_energy
            eb2 = sphere_energies.blob2_energy

            # Compute asymmetry, handle division by zero
            total = eb1 + eb2
            db = total > 0 ? abs(eb1 - eb2) / total : 0.0
            DB[ievent_idx, ir] = db
        end
    end

    return (vec(mean(DB, dims=1)), vec(std(DB, dims=1)))
end


"""
    blob_analysis_vs_radius(track, walk_result; r0::Int=5, rl::Int=15)

Compute blob energies, asymmetry and track length for a single track across a range of radii.

# Arguments
- `track`: A Tracks object
- `walk_result`: Result from walk_track_from_extremes
- `r0::Int`: Minimum sphere radius in mm (default: 5)
- `rl::Int`: Maximum sphere radius in mm (default: 15)

# Returns
NamedTuple with vectors indexed by radius:
- `EB1::Vector{Float64}`: Energy of higher energy blob (MeV)
- `EB2::Vector{Float64}`: Energy of lower energy blob (MeV)
- `DB::Vector{Float64}`: Blob asymmetry |eb1-eb2|/(eb1+eb2)
- `DL::Float64`: Drift/track length from walk_result (mm)
- `radii::Vector{Int}`: The radius values used
"""
function blob_analysis_vs_radius(track, walk_result; r0::Int=5, rl::Int=15)
    nradii = rl - r0 + 1
    EB1 = Vector{Float64}(undef, nradii)
    EB2 = Vector{Float64}(undef, nradii)
    DB = Vector{Float64}(undef, nradii)

    for (ir, rr) in enumerate(r0:rl)
        sphere_energies = energy_in_spheres_around_extremes(track, walk_result, Float64(rr))
        eb1 = sphere_energies.blob1_energy
        eb2 = sphere_energies.blob2_energy

        EB1[ir] = eb1
        EB2[ir] = eb2

        # Compute asymmetry
        total = eb1 + eb2
        DB[ir] = total > 0 ? abs(eb1 - eb2) / total : 0.0
    end

    DL = walk_result.total_length

    return (EB1=EB1, EB2=EB2, DB=DB, DL=DL, radii=collect(r0:rl))
end


"""
    blob_analysis_vs_radius(track; r0::Int=5, rl::Int=15)

Convenience method that computes walk_result internally.
"""
function blob_analysis_vs_radius(track; r0::Int=5, rl::Int=15)
    walk_result = walk_track_from_extremes(track)
    return blob_analysis_vs_radius(track, walk_result; r0=r0, rl=rl)
end


"""
    blob_asymmetry_from_path(tracks::Vector{Tracks}; i0::Int=1, il::Int=10,
                             n0::Int=1, nl::Int=10, energy_col::Symbol=:energy)

Compute blob energy asymmetry |eb1-eb2|/(eb1+eb2) using path voxels directly.

Uses `energy_blobs_from_path` which takes n voxels from each end of the walk path,
rather than spheres of a given radius around extremes.

# Arguments
- `tracks::Vector{Tracks}`: Vector of track objects
- `i0::Int`: First event index (default: 1)
- `il::Int`: Last event index (default: 10)
- `n0::Int`: Minimum number of voxels from each end (default: 1)
- `nl::Int`: Maximum number of voxels from each end (default: 10)
- `energy_col::Symbol`: Column to use for energy (:energy or :electrons, default: :energy)

# Returns
- `(mean_asymmetry, std_asymmetry)`: Vectors of mean and std of asymmetry for each n value
"""
function blob_asymmetry_from_path(tracks::Vector{Tracks}; i0::Int=1, il::Int=10,
                                   n0::Int=1, nl::Int=10, energy_col::Symbol=:energy)
    # Validate indices
    il = min(il, length(tracks))
    i0 = max(1, i0)

    if i0 > il
        error("i0 ($i0) must be <= il ($il)")
    end

    nevents = il - i0 + 1
    nn = nl - n0 + 1
    DB = Matrix{Float64}(undef, nevents, nn)

    for (ievent_idx, ievent) in enumerate(i0:il)
        trk = tracks[ievent]
        walk_result = walk_track_from_extremes(trk)

        for (in_idx, n) in enumerate(n0:nl)
            blobs = energy_blobs_from_path(walk_result, n; energy_col=energy_col)
            eb1 = blobs.eb1
            eb2 = blobs.eb2

            # Compute asymmetry, handle division by zero
            total = eb1 + eb2
            db = total > 0 ? abs(eb1 - eb2) / total : 0.0
            DB[ievent_idx, in_idx] = db
        end
    end

    return (vec(mean(DB, dims=1)), vec(std(DB, dims=1)))
end


"""
    energy_in_variable_spheres_around_extremes(track::Tracks, walk_result;
                                                seed_radius::Float64=3.0,
                                                step::Float64=1.0,
                                                max_radius::Float64=10.0,
                                                threshold::Float64=0.05)

Calculate energy in spheres around track extremes using adaptive radius expansion.

Starts with a seed radius and expands until energy variation becomes small or max radius is reached.

# Arguments
- `track::Tracks`: A Tracks object
- `walk_result`: Result from walk_track_from_extremes function
- `seed_radius::Float64`: Initial sphere radius in mm (default: 3.0)
- `step::Float64`: Radius increment in mm (default: 1.0)
- `max_radius::Float64`: Maximum sphere radius in mm (default: 10.0)
- `threshold::Float64`: Relative energy change threshold to stop expansion (default: 0.05)

# Returns
- NamedTuple with:
  - blob1_energy: Energy in sphere with higher energy
  - blob2_energy: Energy in sphere with lower energy
  - blob1_radius: Final radius used for blob1
  - blob2_radius: Final radius used for blob2
  - blob1_voxel_count: Number of voxels in blob1
  - blob2_voxel_count: Number of voxels in blob2
  - blob1_center: Coordinates of blob1 center (x, y, z)
  - blob2_center: Coordinates of blob2 center (x, y, z)
  - blob1_history: Vector of (radius, energy) pairs showing expansion history
  - blob2_history: Vector of (radius, energy) pairs showing expansion history

# Algorithm
For each extreme:
1. Start with seed_radius
2. Calculate energy within sphere
3. Expand radius by step
4. Calculate new energy
5. If relative change < threshold, stop
6. If radius >= max_radius, stop
7. Repeat from step 3
"""
function energy_in_variable_spheres_around_extremes(track::Tracks, walk_result;
                                                     seed_radius::Float64=3.0,
                                                     step::Float64=1.0,
                                                     max_radius::Float64=10.0,
                                                     threshold::Float64=0.05)
    # Check if we have valid extremes
    if isnothing(walk_result.extremes[1])
        return (blob1_energy = 0.0,
                blob2_energy = 0.0,
                blob1_radius = 0.0,
                blob2_radius = 0.0,
                blob1_voxel_count = 0,
                blob2_voxel_count = 0,
                blob1_center = (NaN, NaN, NaN),
                blob2_center = (NaN, NaN, NaN),
                blob1_history = Tuple{Float64,Float64}[],
                blob2_history = Tuple{Float64,Float64}[])
    end

    start_voxel, end_voxel = walk_result.extremes

    # Get coordinates of extreme voxels
    start_center = (start_voxel.x, start_voxel.y, start_voxel.z)
    end_center = (end_voxel.x, end_voxel.y, end_voxel.z)

    # Helper function to calculate energy in sphere at given radius
    function calculate_sphere_energy(center, radius)
        energy = 0.0
        count = 0
        for i in 1:nrow(track.voxels)
            voxel_x = track.voxels.x[i]
            voxel_y = track.voxels.y[i]
            voxel_z = track.voxels.z[i]
            voxel_energy = track.voxels.energy[i]

            dist = euclidean_distance(voxel_x, voxel_y, voxel_z,
                                     center[1], center[2], center[3])

            if dist <= radius
                energy += voxel_energy
                count += 1
            end
        end
        return energy, count
    end

    # Adaptive expansion for start sphere
    start_history = Tuple{Float64,Float64}[]
    current_radius = seed_radius
    prev_energy = 0.0
    start_final_energy = 0.0
    start_final_count = 0
    start_final_radius = seed_radius

    while current_radius <= max_radius
        energy, count = calculate_sphere_energy(start_center, current_radius)
        push!(start_history, (current_radius, energy))

        # Check convergence
        if energy > 0.0 && prev_energy > 0.0
            relative_change = abs(energy - prev_energy) / energy
            if relative_change < threshold
                start_final_energy = energy
                start_final_count = count
                start_final_radius = current_radius
                break
            end
        end

        start_final_energy = energy
        start_final_count = count
        start_final_radius = current_radius
        prev_energy = energy
        current_radius += step
    end

    # Adaptive expansion for end sphere
    end_history = Tuple{Float64,Float64}[]
    current_radius = seed_radius
    prev_energy = 0.0
    end_final_energy = 0.0
    end_final_count = 0
    end_final_radius = seed_radius

    while current_radius <= max_radius
        energy, count = calculate_sphere_energy(end_center, current_radius)
        push!(end_history, (current_radius, energy))

        # Check convergence
        if energy > 0.0 && prev_energy > 0.0
            relative_change = abs(energy - prev_energy) / energy
            if relative_change < threshold
                end_final_energy = energy
                end_final_count = count
                end_final_radius = current_radius
                break
            end
        end

        end_final_energy = energy
        end_final_count = count
        end_final_radius = current_radius
        prev_energy = energy
        current_radius += step
    end

    # Determine which is blob1 (higher energy) and blob2 (lower energy)
    if start_final_energy >= end_final_energy
        blob1_energy = start_final_energy
        blob1_radius = start_final_radius
        blob1_voxel_count = start_final_count
        blob1_center = start_center
        blob1_history = start_history
        blob2_energy = end_final_energy
        blob2_radius = end_final_radius
        blob2_voxel_count = end_final_count
        blob2_center = end_center
        blob2_history = end_history
    else
        blob1_energy = end_final_energy
        blob1_radius = end_final_radius
        blob1_voxel_count = end_final_count
        blob1_center = end_center
        blob1_history = end_history
        blob2_energy = start_final_energy
        blob2_radius = start_final_radius
        blob2_voxel_count = start_final_count
        blob2_center = start_center
        blob2_history = start_history
    end

    return (blob1_energy = blob1_energy,
            blob2_energy = blob2_energy,
            blob1_radius = blob1_radius,
            blob2_radius = blob2_radius,
            blob1_voxel_count = blob1_voxel_count,
            blob2_voxel_count = blob2_voxel_count,
            blob1_center = blob1_center,
            blob2_center = blob2_center,
            blob1_history = blob1_history,
            blob2_history = blob2_history)
end


"""
    energy_in_spheres_around_extremes(track, radius)

Convenience wrapper: computes walk_result internally.
"""
function energy_in_spheres_around_extremes(track::Tracks, radius::Float64)
    walk_result = walk_track_from_extremes(track)
    return energy_in_spheres_around_extremes(track, walk_result, radius)
end

