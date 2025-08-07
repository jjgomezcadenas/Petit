
function euclidean_distance(x1, y1, z1, x2, y2, z2)
    return sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
end


function voxel_distances(hitsdf::DataFrame, event_id::Int; max_distance=Inf)
    evtdf = get_event(hitsdf, event_id)
    voxel_distances(evtdf; max_distance=max_distance)
end
    
    
function voxel_distances(df::DataFrame; max_distance=Inf)
    all_distances = Float64[]

    n_voxels = nrow(df)
    
    # Return empty array for events with insufficient voxels
    if n_voxels < 2
        return Float64[]
    end
    
    # Compute all pairwise distances within this event
    for i in 1:(n_voxels-1)
        for j in (i+1):n_voxels
            dist = euclidean_distance(
                df.x[i], df.y[i], df.z[i],
                df.x[j], df.y[j], df.z[j]
            )
            
            # Only include distances within the specified maximum
            if dist <= max_distance
                push!(all_distances, dist)
            end
        end
    end

    return all_distances
end


function voxel_closest_distance(hitsdf::DataFrame, event_id::Int; max_distance=Inf)
    evtdf = get_event(hitsdf, event_id)
    voxel_closest_distance(evtdf; max_distance=max_distance)
end


function voxel_closest_distance(df::DataFrame; max_distance=Inf)

    closest_distances = Float64[]
    
    n_voxels = nrow(df)
    
    # Return empty array for events with insufficient voxels
    if n_voxels < 2
        return Float64[]
    end
        
    # For each voxel, find the closest neighbor
    for i in 1:n_voxels
        min_dist = Inf
        
        # Check distance to all other voxels in the same event
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
        
        # Only include distances within the specified maximum
        if min_dist <= max_distance && min_dist != Inf
            push!(closest_distances, min_dist)
        end
    end
    
    return closest_distances
end


function voxel_energy(hitsdf::DataFrame, event_id::Int)
    evtdf = get_event(hitsdf, event_id)
    voxel_energy(evtdf)
end


function voxel_energy(df::DataFrame)
    # Simply return the energy column as Float64 array
    return Float64.(df.energy)
end


function voxelize_hits(hitsdf::DataFrame, voxel_size::Float64)
    ghdf = groupby(hitsdf, :event_id)
    voxel_rows = Vector{NamedTuple{(:event_id, :x, :y, :z, :energy), Tuple{Int64, Float64, Float64, Float64, Float64}}}()

    for group in ghdf
        eid = first(group.event_id)
        
        x = Float64.(group.x)
        y = Float64.(group.y)
        z = Float64.(group.z)
        e = Float64.(group.energy)

        # Compute voxel indices (integers), then voxel center positions
        ix = floor.(Int, x ./ voxel_size)
        iy = floor.(Int, y ./ voxel_size)
        iz = floor.(Int, z ./ voxel_size)

        # Combine into a tuple key
        voxel_keys = zip(ix, iy, iz)

        # Aggregate energy per voxel
        voxel_energy = Dict{Tuple{Int, Int, Int}, Float64}()

        for (key, en) in zip(voxel_keys, e)
            voxel_energy[key] = get(voxel_energy, key, 0.0) + en
        end

        # Push one row per voxel
        for ((i, j, k), total_e) in voxel_energy
            push!(voxel_rows, (
                event_id = eid,
                x = (i + 0.5) * voxel_size,
                y = (j + 0.5) * voxel_size,
                z = (k + 0.5) * voxel_size,
                energy = total_e
            ))
        end
    end

    return DataFrame(voxel_rows)
end

