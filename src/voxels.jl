using NearestNeighbors

"""
    euclidean_distance(x1, y1, z1, x2, y2, z2)

Compute 3D Euclidean distance between two points.
"""
function euclidean_distance(x1, y1, z1, x2, y2, z2)
    return sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
end


"""
    voxel_distances(hitsdf, event_id; max_distance=Inf)

Compute all pairwise distances for voxels in a specific event.
"""
function voxel_distances(hitsdf::DataFrame, event_id::Int; max_distance=Inf)
    evtdf = get_event(hitsdf, event_id)
    voxel_distances(evtdf; max_distance=max_distance)
end


"""
    voxel_distances(df; max_distance=Inf)

Compute all pairwise distances between voxels. O(n²) complexity.

# Returns
- `Vector{Float64}`: All pairwise distances <= max_distance
"""
function voxel_distances(df::DataFrame; max_distance=Inf)
    all_distances = Float64[]
    n_voxels = nrow(df)

    if n_voxels < 2
        return Float64[]
    end

    for i in 1:(n_voxels-1)
        for j in (i+1):n_voxels
            dist = euclidean_distance(
                df.x[i], df.y[i], df.z[i],
                df.x[j], df.y[j], df.z[j]
            )
            if dist <= max_distance
                push!(all_distances, dist)
            end
        end
    end

    return all_distances
end


"""
    voxel_closest_distance(hitsdf, event_id; max_distance=Inf)

Compute closest neighbor distance for each voxel in a specific event.
"""
function voxel_closest_distance(hitsdf::DataFrame, event_id::Int; max_distance=Inf)
    evtdf = get_event(hitsdf, event_id)
    voxel_closest_distance(evtdf; max_distance=max_distance)
end


"""
    voxel_closest_distance(df; max_distance=Inf)

For each voxel, find distance to its nearest neighbor. O(n²) complexity.

# Returns
- `Vector{Float64}`: Closest neighbor distance for each voxel
"""
function voxel_closest_distance(df::DataFrame; max_distance=Inf)
    closest_distances = Float64[]
    n_voxels = nrow(df)

    if n_voxels < 2
        return Float64[]
    end

    for i in 1:n_voxels
        min_dist = Inf
        for j in 1:n_voxels
            if i != j
                dist = euclidean_distance(
                    df.x[i], df.y[i], df.z[i],
                    df.x[j], df.y[j], df.z[j]
                )
                if dist < min_dist
                    min_dist = dist
                end
            end
        end

        if min_dist <= max_distance && min_dist != Inf
            push!(closest_distances, min_dist)
        end
    end

    return closest_distances
end


"""
    closest_voxel_distances_fast(df)

Fast O(n log n) version using KDTree. Returns closest neighbor distance for each voxel.
"""
function closest_voxel_distances_fast(df::DataFrame)
    n = nrow(df)
    n < 2 && return Float64[]

    points = reduce(hcat, [[df.x[i], df.y[i], df.z[i]] for i in 1:n])
    tree = KDTree(points)
    _, dists = knn(tree, points, 2, true)

    return [d[2] for d in dists]
end


"""
    closest_voxel_distances_fast(hitsdf, event_id)

Fast version: compute closest neighbor distances for a specific event using KDTree.
"""
function closest_voxel_distances_fast(hitsdf::DataFrame, event_id::Int)
    evtdf = get_event(hitsdf, event_id)
    closest_voxel_distances_fast(evtdf)
end

"""
    voxel_energy(hitsdf, event_id; energy_column=:energy)

Extract energy values for a specific event.
"""
function voxel_energy(hitsdf::DataFrame, event_id::Int; energy_column::Symbol=:energy)
    evtdf = get_event(hitsdf, event_id)
    voxel_energy(evtdf; energy_column=energy_column)
end


"""
    voxel_energy(df; energy_column=:energy)

Return energy column as Float64 vector.
"""
function voxel_energy(df::DataFrame; energy_column::Symbol=:energy)
    return Float64.(df[!, energy_column])
end


"""
    voxelize_event(evtdf, voxel_size; ef=1e+5/2.5)

Voxelize a single event DataFrame.

# Arguments
- `evtdf::DataFrame`: Event data with x, y, z, energy columns
- `voxel_size::Float64`: Voxel size in mm
- `ef::Float64`: Energy to electrons conversion factor (default: 1e5/2.5)

# Returns
- `DataFrame`: Voxelized data with event_id, x, y, z, energy, electrons
"""
function voxelize_event(evtdf::DataFrame, voxel_size::Float64; ef::Float64=1e+5/2.5)

    # Get event_id from DataFrame (assume all rows have same event_id)
    event_id = nrow(evtdf) > 0 ? evtdf.event_id[1] : 0

    # Extract coordinates and energy
    x = Float64.(evtdf.x)
    y = Float64.(evtdf.y)
    z = Float64.(evtdf.z)
    e = Float64.(evtdf.energy)

    # Compute voxel indices
    ix = floor.(Int, x ./ voxel_size)
    iy = floor.(Int, y ./ voxel_size)
    iz = floor.(Int, z ./ voxel_size)

    # Aggregate energy per voxel
    voxel_energy = Dict{Tuple{Int, Int, Int}, Float64}()
    for (key, en) in zip(zip(ix, iy, iz), e)
        voxel_energy[key] = get(voxel_energy, key, 0.0) + en
    end

    # Build output rows
    voxel_rows = Vector{NamedTuple{(:event_id, :x, :y, :z, :energy, :electrons),
                                    Tuple{Int64, Float64, Float64, Float64, Float64, Int64}}}()

    for ((i, j, k), total_e) in voxel_energy
        push!(voxel_rows, (
            event_id = event_id,
            x = (i + 0.5) * voxel_size,
            y = (j + 0.5) * voxel_size,
            z = (k + 0.5) * voxel_size,
            energy = total_e,
            electrons = Int(round(total_e * ef))
        ))
    end

    return DataFrame(voxel_rows)
end


"""
    voxelize_event(hitsdf, event_number, voxel_size; ef=1e+5/2.5)

Voxelize a specific event from a hits DataFrame.
"""
function voxelize_event(hitsdf::DataFrame, event_number::Int, voxel_size::Float64;
                        ef::Float64=1e+5/2.5)
    evtdf = get_event(hitsdf, event_number)
    voxelize_event(evtdf, voxel_size; ef=ef)
end


"""
    voxelize_hits(hitsdf, voxel_size; ef=1e+5/2.5)

Voxelize all events in a hits DataFrame.

# Returns
- `DataFrame`: Combined voxelized data for all events
"""
function voxelize_hits(hitsdf::DataFrame, voxel_size::Float64; ef::Float64=1e+5/2.5)
    event_ids = unique(hitsdf.event_id)
    voxelized_events = [voxelize_event(hitsdf, eid, voxel_size; ef=ef) for eid in event_ids]
    return vcat(voxelized_events...)
end

