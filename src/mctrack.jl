"""
    MC Track Functions

Functions for processing Monte Carlo truth information from particle tracks.
Handles both single-electron and double-beta (two-electron) events.
"""

"""
    get_mc_extremes(event_df)

Get MC truth extremes from event hits DataFrame.

For single electron (particle_id=1 only): returns start and Bragg peak positions.
For double-beta (particle_id=1 and 2): returns the two Bragg peak positions.

Returns (pos1, pos2) tuples or (nothing, nothing) if no primary hits found.
"""
function get_mc_extremes(event_df::DataFrame)
    # Check for double-beta (two primary particles)
    particle_ids = unique(event_df.particle_id)
    has_two_particles = (1 in particle_ids) && (2 in particle_ids)

    if has_two_particles
        # Double-beta: return Bragg peaks of both electrons
        hits1 = filter(row -> row.particle_id == 1, event_df)
        hits2 = filter(row -> row.particle_id == 2, event_df)

        if nrow(hits1) == 0 || nrow(hits2) == 0
            return nothing, nothing
        end

        sort!(hits1, :time)
        sort!(hits2, :time)

        # Last hit of each particle = Bragg peak
        bragg1 = last(hits1)
        bragg2 = last(hits2)

        pos1 = (Float64(bragg1.x), Float64(bragg1.y), Float64(bragg1.z))
        pos2 = (Float64(bragg2.x), Float64(bragg2.y), Float64(bragg2.z))

        return pos1, pos2
    else
        # Single electron: particle_id = 1
        primary_hits = filter(row -> row.particle_id == 1, event_df)

        if nrow(primary_hits) == 0
            return nothing, nothing
        end

        sort!(primary_hits, :time)

        first_hit = first(primary_hits)
        last_hit = last(primary_hits)

        pos1 = (Float64(first_hit.x), Float64(first_hit.y), Float64(first_hit.z))
        pos2 = (Float64(last_hit.x), Float64(last_hit.y), Float64(last_hit.z))

        return pos1, pos2
    end
end

"""
    voxelize_particle_hits(hits_df, mcvox_size)

Helper function to voxelize hits from a single particle.
Returns vectors (xs, ys, zs, Es) of voxelized positions and energies.
"""
function voxelize_particle_hits(hits_df::DataFrame, mcvox_size::Float64)
    if nrow(hits_df) == 0
        return Float64[], Float64[], Float64[], Float64[]
    end

    hits = copy(hits_df)
    hits.vox_ix = floor.(Int, hits.x ./ mcvox_size)
    hits.vox_iy = floor.(Int, hits.y ./ mcvox_size)
    hits.vox_iz = floor.(Int, hits.z ./ mcvox_size)

    # Group consecutive hits in same voxel
    n_hits = nrow(hits)
    voxel_groups = Vector{UnitRange{Int}}()

    i = 1
    while i <= n_hits
        start_i = i
        current_vox = (hits.vox_ix[i], hits.vox_iy[i], hits.vox_iz[i])

        while i <= n_hits &&
              (hits.vox_ix[i], hits.vox_iy[i], hits.vox_iz[i]) == current_vox
            i += 1
        end

        push!(voxel_groups, start_i:(i-1))
    end

    # Compute voxel positions and energies
    n_voxels = length(voxel_groups)
    xs = Vector{Float64}(undef, n_voxels)
    ys = Vector{Float64}(undef, n_voxels)
    zs = Vector{Float64}(undef, n_voxels)
    Es = Vector{Float64}(undef, n_voxels)

    for (k, group) in enumerate(voxel_groups)
        hits_in_group = @view hits[group, :]
        total_E = sum(hits_in_group.energy)

        if total_E > 0
            xs[k] = sum(hits_in_group.x .* hits_in_group.energy) / total_E
            ys[k] = sum(hits_in_group.y .* hits_in_group.energy) / total_E
            zs[k] = sum(hits_in_group.z .* hits_in_group.energy) / total_E
        else
            xs[k] = mean(hits_in_group.x)
            ys[k] = mean(hits_in_group.y)
            zs[k] = mean(hits_in_group.z)
        end
        Es[k] = total_E
    end

    return xs, ys, zs, Es
end

"""
    compute_mc_path(event_df, mcvox_size)

Compute the Monte Carlo path by voxelizing primary particle hits.

For single electron: path follows particle_id=1 in time order.
For double-beta: path goes from Bragg_peak_1 → vertex → Bragg_peak_2
  (electron 1 reversed, then electron 2 forward)

Returns a DataFrame with columns: x, y, z, energy, s (arc-length), primary_electron
  - primary_electron: 1 for all voxels in single-electron events
                      1 or 2 for double-beta, indicating which electron
"""
function compute_mc_path(event_df::DataFrame, mcvox_size::Float64)
    # Check for double-beta (two primary particles)
    particle_ids = unique(event_df.particle_id)
    has_two_particles = (1 in particle_ids) && (2 in particle_ids)

    if has_two_particles
        # Double-beta decay: combine both electrons
        hits1 = filter(row -> row.particle_id == 1, event_df)
        hits2 = filter(row -> row.particle_id == 2, event_df)

        if nrow(hits1) == 0 || nrow(hits2) == 0
            return DataFrame(x=Float64[], y=Float64[], z=Float64[],
                            energy=Float64[], s=Float64[], primary_electron=Int[])
        end

        # Sort each by time
        sort!(hits1, :time)
        sort!(hits2, :time)

        # Voxelize each particle separately
        x1, y1, z1, E1 = voxelize_particle_hits(hits1, mcvox_size)
        x2, y2, z2, E2 = voxelize_particle_hits(hits2, mcvox_size)

        if isempty(x1) || isempty(x2)
            return DataFrame(x=Float64[], y=Float64[], z=Float64[],
                            energy=Float64[], s=Float64[], primary_electron=Int[])
        end

        # Create continuous path: electron1 reversed (Bragg→vertex) + electron2 (vertex→Bragg)
        # Reverse electron 1
        x1_rev = reverse(x1)
        y1_rev = reverse(y1)
        z1_rev = reverse(z1)
        E1_rev = reverse(E1)

        # Concatenate: electron1_reversed + electron2
        xs = vcat(x1_rev, x2)
        ys = vcat(y1_rev, y2)
        zs = vcat(z1_rev, z2)
        Es = vcat(E1_rev, E2)

        # Track which electron each voxel belongs to
        pe = vcat(fill(1, length(x1)), fill(2, length(x2)))
    else
        # Single electron: particle_id = 1
        primary_hits = filter(row -> row.particle_id == 1, event_df)

        if nrow(primary_hits) == 0
            return DataFrame(x=Float64[], y=Float64[], z=Float64[],
                            energy=Float64[], s=Float64[], primary_electron=Int[])
        end

        sort!(primary_hits, :time)
        xs, ys, zs, Es = voxelize_particle_hits(primary_hits, mcvox_size)

        if isempty(xs)
            return DataFrame(x=Float64[], y=Float64[], z=Float64[],
                            energy=Float64[], s=Float64[], primary_electron=Int[])
        end

        # All voxels belong to electron 1
        pe = fill(1, length(xs))
    end

    # Compute cumulative arc-length
    n_voxels = length(xs)
    ss = Vector{Float64}(undef, n_voxels)
    ss[1] = 0.0
    for k in 2:n_voxels
        dx = xs[k] - xs[k-1]
        dy = ys[k] - ys[k-1]
        dz = zs[k] - zs[k-1]
        ss[k] = ss[k-1] + sqrt(dx^2 + dy^2 + dz^2)
    end

    return DataFrame(x=xs, y=ys, z=zs, energy=Es, s=ss, primary_electron=pe)
end
