using LinearAlgebra
using Graphs
using NearestNeighbors
using SparseArrays


"""
    DiffusionParams

Diffusion and voxelization parameters for a track.

# Fields
- `ldrift::Float64`: Drift length (mm)
- `sigma_t::Float64`: Transverse diffusion (mm)
- `sigma_l::Float64`: Longitudinal diffusion (mm)
- `voxel_size::Float64`: Voxel size (mm)
- `max_distance::Float64`: Max distance for clustering (mm)
- `energy_threshold::Float64`: Energy threshold (keV)
- `nbins_df::Int`: Number of bins for diffusion histogram
- `nsigma_df::Float64`: Number of sigmas for histogram padding
"""
struct DiffusionParams
    ldrift::Float64
    sigma_t::Float64
    sigma_l::Float64
    voxel_size::Float64
    max_distance::Float64
    energy_threshold::Float64
    nbins_df::Int
    nsigma_df::Float64
end

"""
    DiffusionParams()

Create DiffusionParams with default values (for MC tracks with no diffusion).
"""
DiffusionParams() = DiffusionParams(0.0, 0.0, 0.0, 1.0, 1.5, 0.0, 100, 3.0)


"""
    Tracks

A track structure representing connected voxels in a particle track.

# Fields
- `voxels::DataFrame`: Voxel data (x, y, z, energy, electrons)
- `graph::SimpleGraph{Int}`: Graph connecting nearby voxels
- `components::Vector{Vector{Int}}`: Connected component indices
- `diffusion::DiffusionParams`: Diffusion and voxelization parameters
"""
struct Tracks
    voxels::DataFrame
    graph::SimpleGraph{Int}
    components::Vector{Vector{Int}}
    diffusion::DiffusionParams
end


"""
    build_tracks(hitsdf, event_id; max_distance=1.5, energy_threshold=0.0, diffusion=DiffusionParams())

Build tracks from hits DataFrame for a specific event.
Extracts event data and delegates to the DataFrame version.
"""
function build_tracks(hitsdf::DataFrame, event_id::Int;
                     max_distance::Float64=1.5,
                     energy_threshold::Float64=0.0,
                     diffusion::DiffusionParams=DiffusionParams(), method="KDT")
    event_data = get_event(hitsdf, event_id)

    if method == "Eucledian"
        build_tracks(event_data; max_distance=max_distance, energy_threshold=energy_threshold, diffusion=diffusion)
    elseif method == "KDT"
        build_tracks_kdtree(event_data; max_distance=max_distance, energy_threshold=energy_threshold, diffusion=diffusion)
    end
end


"""
    build_tracks(event_data; max_distance=1.5, energy_threshold=0.0, diffusion=DiffusionParams())

Build tracks from voxelized event data.

# Arguments
- `event_data::DataFrame`: Voxelized data with x, y, z, energy columns
- `max_distance::Float64`: Max distance (mm) to connect voxels (default: 1.5)
- `energy_threshold::Float64`: Min energy (MeV) to include voxel (default: 0.0)
- `diffusion::DiffusionParams`: Diffusion parameters for the tracks

# Returns
- `Vector{Tracks}`: List of tracks (connected components), sorted by energy
"""
function build_tracks(event_data::DataFrame;
                     max_distance::Float64=1.5,
                     energy_threshold::Float64=0.0,
                     diffusion::DiffusionParams=DiffusionParams())

    n_voxels = nrow(event_data)
    if n_voxels == 0
        return Tracks[]
    end

    # Filter voxels by energy threshold
    valid_indices = findall(e -> e >= energy_threshold, event_data.energy)
    if isempty(valid_indices)
        return Tracks[]
    end

    filtered_data = event_data[valid_indices, :]
    n_filtered = length(valid_indices)

    # Build graph: nodes are voxels, edges connect voxels within max_distance
    g = SimpleGraph(n_filtered)
    for i in 1:(n_filtered-1)
        for j in (i+1):n_filtered
            dist = euclidean_distance(
                filtered_data.x[i], filtered_data.y[i], filtered_data.z[i],
                filtered_data.x[j], filtered_data.y[j], filtered_data.z[j]
            )
            if dist <= max_distance
                add_edge!(g, i, j)
            end
        end
    end

    # Each connected component becomes a separate track
    components = connected_components(g)
    graphs = Tracks[]

    for component in components
        if !isempty(component)
            component_data = filtered_data[component, :]

            # Create subgraph with remapped vertex indices
            subgraph = SimpleGraph(length(component))
            vertex_map = Dict(old_v => new_v for (new_v, old_v) in enumerate(component))

            for edge in edges(g)
                src_old, dst_old = src(edge), dst(edge)
                if src_old in component && dst_old in component
                    add_edge!(subgraph, vertex_map[src_old], vertex_map[dst_old])
                end
            end

            push!(graphs, Tracks(component_data, subgraph, [component], diffusion))
        end
    end

    return graphs
end


### Using kdtree


# helper: coord matrix 3×N
function _coords_from_df(df::DataFrame)
    n = nrow(df)
    coords = Matrix{Float64}(undef, 3, n)
    @inbounds for i in 1:n
        coords[1, i] = df.x[i]
        coords[2, i] = df.y[i]
        coords[3, i] = df.z[i]
    end
    return coords
end

function build_tracks_kdtree(event_data::DataFrame;
                             max_distance::Float64 = 1.5,
                             energy_threshold::Float64 = 0.0,
                             diffusion::DiffusionParams = DiffusionParams())

    n_voxels = nrow(event_data)
    n_voxels == 0 && return Tracks[]

    # filter by energy
    valid_indices = findall(e -> e >= energy_threshold, event_data.energy)
    isempty(valid_indices) && return Tracks[]

    filtered_data = event_data[valid_indices, :]
    n_filtered = nrow(filtered_data)

    # trivial case
    n_filtered == 1 && return [Tracks(filtered_data, SimpleGraph(1), [collect(1:1)], diffusion)]

    # coords and KDTree
    coords = _coords_from_df(filtered_data)
    tree   = KDTree(coords)

    g = SimpleGraph(n_filtered)
    r  = max_distance

    # build graph using neighbors in range
    @inbounds for i in 1:n_filtered
        idxs = inrange(tree, coords[:, i], r, true)
        for j in idxs
            j == i && continue
            j <  i && continue   # avoid double edges; we'll see (i,j) once with i<j
            add_edge!(g, i, j)
        end
    end

    # components → tracks
    components = connected_components(g)
    tracks = Tracks[]

    for comp in components
        isempty(comp) && continue

        comp_data = filtered_data[comp, :]

        # induced subgraph for this component (returns tuple: subgraph, vertex_map)
        subg, _ = induced_subgraph(g, comp)

        # Tracks constructor: data, graph, list of components, diffusion params
        push!(tracks, Tracks(comp_data, subg, [comp], diffusion))
    end

    return tracks
end


"""
    build_tracks_knn(event_data::DataFrame; k=10, max_distance=Inf,
                     energy_threshold=0.0, diffusion=DiffusionParams(),
                     mutual=false)

Build tracks using k-Nearest Neighbor graph instead of radius graph.

# Arguments
- `event_data::DataFrame`: Voxelized event data with x, y, z, energy columns
- `k::Int`: Number of nearest neighbors per voxel (default: 10)
- `max_distance::Float64`: Maximum edge length allowed (default: Inf, no limit)
- `energy_threshold::Float64`: Minimum voxel energy in MeV (default: 0.0)
- `diffusion::DiffusionParams`: Diffusion parameters for the tracks
- `mutual::Bool`: If true, keep only mutual kNN edges (default: false)

# Returns
- `Vector{Tracks}`: Vector of tracks (connected components)

# Notes
kNN graph prevents shortcut edges across U-shaped bends by limiting each
voxel to its k closest neighbors, rather than all neighbors within a radius.
This fixes extreme-finding failures on curved tracks.
"""
function build_tracks_knn(event_data::DataFrame;
                          k::Int = 10,
                          max_distance::Float64 = Inf,
                          energy_threshold::Float64 = 0.0,
                          diffusion::DiffusionParams = DiffusionParams(),
                          mutual::Bool = false)

    n_voxels = nrow(event_data)
    n_voxels == 0 && return Tracks[]

    # Filter by energy
    valid_indices = findall(e -> e >= energy_threshold, event_data.energy)
    isempty(valid_indices) && return Tracks[]

    filtered_data = event_data[valid_indices, :]
    n_filtered = nrow(filtered_data)

    # Trivial case: single voxel
    n_filtered == 1 && return [Tracks(filtered_data, SimpleGraph(1), [collect(1:1)], diffusion)]

    # Adjust k if we have fewer points
    k_actual = min(k, n_filtered - 1)
    k_actual < 1 && return [Tracks(filtered_data, SimpleGraph(n_filtered), [collect(1:n_filtered)], diffusion)]

    # Build KDTree
    coords = _coords_from_df(filtered_data)
    tree = KDTree(coords)

    # Build kNN graph
    g = SimpleGraph(n_filtered)

    if mutual
        # Mutual kNN: keep edge only if both endpoints have each other in their kNN
        # First, collect all kNN sets
        knn_sets = Vector{Set{Int}}(undef, n_filtered)
        @inbounds for i in 1:n_filtered
            idxs, dists = knn(tree, coords[:, i], k_actual + 1, true)
            # Filter by max_distance and exclude self
            neighbors = Set{Int}()
            for (idx, d) in zip(idxs, dists)
                idx == i && continue
                if d <= max_distance
                    push!(neighbors, idx)
                end
            end
            knn_sets[i] = neighbors
        end

        # Add edge only if mutual
        @inbounds for i in 1:n_filtered
            for j in knn_sets[i]
                j <= i && continue  # Avoid duplicate edges
                if i in knn_sets[j]
                    add_edge!(g, i, j)
                end
            end
        end
    else
        # Standard kNN: add edge if j is in kNN(i)
        # Since we process all i, edges become symmetric
        @inbounds for i in 1:n_filtered
            idxs, dists = knn(tree, coords[:, i], k_actual + 1, true)
            for (idx, d) in zip(idxs, dists)
                idx == i && continue
                if d <= max_distance
                    add_edge!(g, i, idx)  # SimpleGraph handles duplicates
                end
            end
        end
    end

    # Connected components → tracks
    components = connected_components(g)
    tracks = Tracks[]

    for comp in components
        isempty(comp) && continue
        comp_data = filtered_data[comp, :]
        subg, _ = induced_subgraph(g, comp)
        push!(tracks, Tracks(comp_data, subg, [comp], diffusion))
    end

    return tracks
end


"""
    make_tracks(hitsdf, event_id; max_distance_mm=10.0, energy_threshold_kev=0.0,
                diffusion=DiffusionParams(), method="KDT", k=10)

Convenience wrapper: build tracks with parameters in mm and keV units.
"""
function make_tracks(hitsdf::DataFrame, event_id::Int;
                    max_distance_mm::Float64=10.0,
                    energy_threshold_kev::Float64=0.0,
                    diffusion::DiffusionParams=DiffusionParams(),
                    method::String="KDT",
                    k::Int=10)
    event_data = get_event(hitsdf, event_id)
    make_tracks(event_data;
                max_distance_mm=max_distance_mm,
                energy_threshold_kev=energy_threshold_kev,
                diffusion=diffusion,
                method=method,
                k=k)
end


"""
    make_tracks(event_data; max_distance_mm=1.0, energy_threshold_kev=0.0,
                diffusion=DiffusionParams(), method="KDT", k=10)

Build and sort tracks by total energy (descending).

# Arguments
- `event_data::DataFrame`: Voxelized event data
- `max_distance_mm::Float64`: Max voxel connection distance in mm
- `energy_threshold_kev::Float64`: Min voxel energy in keV
- `diffusion::DiffusionParams`: Diffusion parameters for the tracks
- `method::String`: Graph construction method:
  - "KDT": Radius graph (all neighbors within max_distance)
  - "kNN": k-Nearest Neighbor graph (k closest neighbors)
  - "kNN_mutual": Mutual kNN (edge only if both are in each other's kNN)
- `k::Int`: Number of neighbors for kNN methods (default: 10)

# Returns
- `Vector{Tracks}`: Tracks sorted by total energy (highest first)
"""
function make_tracks(event_data::DataFrame;
                    max_distance_mm::Float64=1.0,
                    energy_threshold_kev::Float64=0.0,
                    diffusion::DiffusionParams=DiffusionParams(),
                    method::String="KDT",
                    k::Int=10)

    energy_threshold = energy_threshold_kev * 1e-3  # Convert keV to MeV

    if method == "KDT"
        tracks = build_tracks_kdtree(event_data;
                                     max_distance=max_distance_mm,
                                     energy_threshold=energy_threshold,
                                     diffusion=diffusion)
    elseif method == "kNN"
        tracks = build_tracks_knn(event_data;
                                  k=k,
                                  max_distance=max_distance_mm,
                                  energy_threshold=energy_threshold,
                                  diffusion=diffusion,
                                  mutual=false)
    elseif method == "kNN_mutual"
        tracks = build_tracks_knn(event_data;
                                  k=k,
                                  max_distance=max_distance_mm,
                                  energy_threshold=energy_threshold,
                                  diffusion=diffusion,
                                  mutual=true)
    else
        # Fallback to old build_tracks method
        tracks = build_tracks(event_data;
                              max_distance=max_distance_mm,
                              energy_threshold=energy_threshold,
                              diffusion=diffusion)
    end

    # Sort by total energy (descending)
    if length(tracks) > 0
        track_energies = [sum(track.voxels.energy) for track in tracks]
        sorted_indices = sortperm(track_energies, rev=true)
        tracks = tracks[sorted_indices]
    end

    return tracks
end


"""
    find_track_extremes(track::Tracks; use_optimized=true)

Find track endpoints using combined topology/curvature analysis.

# Arguments
- `track::Tracks`: Track object
- `use_optimized::Bool`: If true, use optimized implementation (default: true)

# Returns
- `(extreme1_idx, extreme2_idx, path, confidence)`: Endpoint indices, path between them, confidence score (0-1)
"""
function find_track_extremes(track::Tracks; use_optimized::Bool=true)
    if use_optimized
        return find_track_extremes_opt(track)
    else
        return find_track_extremes_improved(track; method=:combined)
    end
end


"""
    find_path_bfs(g, start_vertex, end_vertex)

Find shortest path between vertices using breadth-first search.

# Returns
- `Vector{Int}`: Vertex indices along the path (empty if no path exists)
"""
function find_path_bfs(g::SimpleGraph, start_vertex::Int, end_vertex::Int)
    if start_vertex == end_vertex
        return [start_vertex]
    end

    # BFS to find shortest path
    queue = [(start_vertex, [start_vertex])]
    visited = Set{Int}()

    while !isempty(queue)
        current, path = popfirst!(queue)

        if current == end_vertex
            return path
        end

        if current in visited
            continue
        end

        push!(visited, current)

        for neighbor in neighbors(g, current)
            if neighbor ∉ visited
                push!(queue, (neighbor, vcat(path, neighbor)))
            end
        end
    end

    # No path found
    return Int[]
end


"""
    walk_track_from_extremes(track::Tracks)

Walk through a track from one endpoint to the other.

# Returns
NamedTuple with:
- `extremes`: (start_voxel, end_voxel) DataFrameRows
- `path_indices`: Vertex indices along the path
- `path_voxels`: DataFrame of voxels in path order
- `total_length`: Path length in mm
- `confidence`: Confidence score (0-1)
"""
function walk_track_from_extremes(track::Tracks)
    extreme1, extreme2, path, confidence = find_track_extremes(track)

    if isnothing(extreme1)
        return (extremes = (nothing, nothing),
                path_indices = Int[],
                path_voxels = DataFrame(),
                total_length = 0.0,
                confidence = 0.0)
    end

    # Get voxel data for extremes
    start_voxel = track.voxels[extreme1, :]
    end_voxel = track.voxels[extreme2, :]

    # Get voxels along the path
    path_voxels = track.voxels[path, :]

    # Calculate total path length
    total_length = 0.0
    for i in 1:length(path)-1
        v1, v2 = path[i], path[i+1]
        total_length += euclidean_distance(
            track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
            track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
        )
    end

    return (extremes = (start_voxel, end_voxel),
            path_indices = path,
            path_voxels = path_voxels,
            total_length = total_length,
            confidence = confidence)
end


"""
    find_track_extremes(trk; i=1)

Convenience function: find extremes for i-th track in a vector.
"""
function find_track_extremes(trk; i=1)
	return walk_track_from_extremes(trk[i])
end


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
    track_energies_keV(tracks::Tracks)

Compute total energy (keV) for each track in the collection.

# Returns
- `Vector{Float64}`: Energy in keV for each track
"""
function track_energies_keV(tracks::Tracks)
	E = Float64[]
	for i in 1:length(tracks)
		energy_kev = 1e+3 * sum(tracks[i].voxels.energy)
		push!(E, energy_kev)
	end
	E
end



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

# =============================================================================
# IMPROVED TRACK EXTREMES ALGORITHMS
# =============================================================================

"""
    find_track_extremes_improved(track; method=:combined)

Find track endpoints using specified algorithm.

# Arguments
- `track::Tracks`: Track object
- `method::Symbol`: Algorithm - `:topology`, `:curvature`, or `:combined`

# Returns
- `(extreme1_idx, extreme2_idx, path, confidence)`: Endpoints, path, and confidence score
"""
function find_track_extremes_improved(track::Tracks; method::Symbol=:combined)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0)
    elseif n_vertices == 1
        return (1, 1, [1], 1.0)
    end

    if method == :topology
        return find_extremes_topology_based(track)
    elseif method == :curvature
        return find_extremes_curvature_based(track)
    else  # :combined
        return find_extremes_combined(track)
    end
end


"""Find extremes using graph topology (degree-1 vertices preferred)."""
function find_extremes_topology_based(track::Tracks)
    g = track.graph
    n_vertices = nv(g)

    # Find vertices with degree 1 (true endpoints)
    endpoints = Int[]
    for v in vertices(g)
        if degree(g, v) == 1
            push!(endpoints, v)
        end
    end

    # If we have exactly 2 endpoints, use them (high confidence)
    if length(endpoints) == 2
        path = find_path_bfs(g, endpoints[1], endpoints[2])
        return (endpoints[1], endpoints[2], path, 0.9)
    end

    # If we have more than 2 endpoints, find the pair with maximum path length
    if length(endpoints) > 2
        max_path_length = 0.0
        best_pair = (endpoints[1], endpoints[2])
        best_path = Int[]

        for i in 1:length(endpoints)-1
            for j in i+1:length(endpoints)
                path = find_path_bfs(g, endpoints[i], endpoints[j])
                path_length = calculate_path_length(track, path)

                if path_length > max_path_length
                    max_path_length = path_length
                    best_pair = (endpoints[i], endpoints[j])
                    best_path = path
                end
            end
        end

        return (best_pair[1], best_pair[2], best_path, 0.8)
    end

    # Fallback: no clear endpoints, use distance-based method
    return find_extremes_distance_fallback(track, 0.6)
end


"""Find extremes using local curvature (low curvature = likely endpoint)."""
function find_extremes_curvature_based(track::Tracks)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices < 3
        return find_extremes_distance_fallback(track, 0.7)
    end

    # Calculate curvature at each vertex
    curvatures = calculate_vertex_curvatures(track)

    # Find vertices with lowest curvature (most "straight")
    # These are likely to be at track ends
    sorted_indices = sortperm(curvatures)

    # Try pairs of low-curvature vertices
    max_path_length = 0.0
    best_pair = (sorted_indices[1], sorted_indices[2])
    best_path = Int[]
    confidence = 0.0

    # Check first several low-curvature vertices
    n_candidates = min(6, n_vertices)

    for i in 1:n_candidates-1
        for j in i+1:n_candidates
            v1, v2 = sorted_indices[i], sorted_indices[j]

            # Skip if vertices are too close (likely same region)
            if are_vertices_too_close(track, v1, v2)
                continue
            end

            path = find_path_bfs(g, v1, v2)
            if !isempty(path)
                path_length = calculate_path_length(track, path)

                # Bonus for low curvature at endpoints
                curvature_bonus = 2.0 / (curvatures[v1] + curvatures[v2] + 1e-6)
                score = path_length + curvature_bonus

                if score > max_path_length
                    max_path_length = score
                    best_pair = (v1, v2)
                    best_path = path
                    confidence = 0.7 + 0.2 * curvature_bonus / 10.0  # Scale confidence
                end
            end
        end
    end

    return (best_pair[1], best_pair[2], best_path, min(confidence, 0.9))
end


"""Combined approach: try multiple methods and pick the best result."""
function find_extremes_combined(track::Tracks)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0)
    end

    # Check if this is a dense track (high connectivity)
    avg_degree = 2 * ne(g) / n_vertices
    is_dense_track = avg_degree > 6.0  # Threshold for dense tracks

    if is_dense_track
        # For dense tracks, use spatial-based method
        spatial_result = find_extremes_spatial_based(track)

        # Also try topology and curvature for comparison
        topo_result = find_extremes_topology_based(track)
        curv_result = find_extremes_curvature_based(track)

        # Compare path lengths and choose the best
        spatial_length = calculate_path_length(track, spatial_result[3])
        topo_length = calculate_path_length(track, topo_result[3])
        curv_length = calculate_path_length(track, curv_result[3])

        # Choose the method that gives the longest path

        if use_mst_fallback
            mst_result = find_extremes_mst_diameter(track, coords)
            mst_length = mst_result[5]
            if mst_length >= max(spatial_length, curv_length, topo_length)
                return (mst_result[1], mst_result[2], mst_result[3], 1.0)
            end
        end

        if spatial_length >= topo_length && spatial_length >= curv_length
            return spatial_result
        elseif topo_length >= curv_length
            return topo_result
        else
            return curv_result
        end
    else
        # For sparse tracks, use original logic
        # Try topology-based first
        topo_result = find_extremes_topology_based(track)

        # Try curvature-based
        curv_result = find_extremes_curvature_based(track)

        # If topology gives high confidence, use it
        if topo_result[4] >= 0.8
            return topo_result
        end

        # Otherwise, compare path lengths and choose longer path
        topo_length = calculate_path_length(track, topo_result[3])
        curv_length = calculate_path_length(track, curv_result[3])

        if curv_length > topo_length * 1.1  # 10% improvement threshold
            return curv_result
        else
            return topo_result
        end
    end
end


"""Find extremes using spatial analysis (for dense/highly-connected tracks)."""
function find_extremes_spatial_based(track::Tracks)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0)
    elseif n_vertices == 1
        return (1, 1, [1], 1.0)
    end

    # Get spatial extreme candidates
    spatial_candidates = [
        argmin(track.voxels.x),  # X min
        argmax(track.voxels.x),  # X max
        argmin(track.voxels.y),  # Y min
        argmax(track.voxels.y),  # Y max
        argmin(track.voxels.z),  # Z min
        argmax(track.voxels.z)   # Z max
    ]

    # Add minimum degree vertices (potential extremes even in dense tracks)
    min_degree = minimum([degree(g, v) for v in vertices(g)])
    min_degree_vertices = [v for v in vertices(g) if degree(g, v) == min_degree]

    # Combine and remove duplicates
    all_candidates = unique([spatial_candidates; min_degree_vertices])

    # Find pair with maximum path length
    max_path_length = 0.0
    best_pair = (all_candidates[1], all_candidates[min(2, end)])
    best_path = Int[]

    for i in 1:length(all_candidates)-1
        for j in i+1:length(all_candidates)
            v1, v2 = all_candidates[i], all_candidates[j]
            path = find_path_bfs(g, v1, v2)
            if !isempty(path)
                path_length = calculate_path_length(track, path)
                if path_length > max_path_length
                    max_path_length = path_length
                    best_pair = (v1, v2)
                    best_path = path
                end
            end
        end
    end

    # Calculate confidence based on path efficiency and characteristics
    if !isempty(best_path)
        pos1 = (track.voxels.x[best_pair[1]], track.voxels.y[best_pair[1]], track.voxels.z[best_pair[1]])
        pos2 = (track.voxels.x[best_pair[2]], track.voxels.y[best_pair[2]], track.voxels.z[best_pair[2]])
        straight_dist = euclidean_distance(pos1[1], pos1[2], pos1[3], pos2[1], pos2[2], pos2[3])

        if max_path_length > 0
            efficiency = straight_dist / max_path_length

            # Base confidence for spatial method
            confidence = 0.75

            # Bonus for good efficiency (straight tracks)
            if efficiency > 0.8
                confidence += 0.1
            end

            # Bonus for long paths
            confidence += min(0.1, max_path_length / 500.0)

            # Penalty if endpoints have very high degree (not true extremes)
            deg1 = degree(g, best_pair[1])
            deg2 = degree(g, best_pair[2])
            avg_endpoint_degree = (deg1 + deg2) / 2.0
            avg_track_degree = 2 * ne(g) / nv(g)

            if avg_endpoint_degree < avg_track_degree * 0.7  # Endpoints have lower connectivity
                confidence += 0.05
            end

            confidence = min(confidence, 0.95)  # Cap at 0.95
        else
            confidence = 0.5
        end
    else
        confidence = 0.3
    end

    return (best_pair[1], best_pair[2], best_path, confidence)
end


"""Calculate local curvature at each vertex (0=straight, 1=sharp turn)."""
function calculate_vertex_curvatures(track::Tracks)
    g = track.graph
    n_vertices = nv(g)
    curvatures = zeros(n_vertices)

    for v in vertices(g)
        neighbors = collect(Graphs.neighbors(g, v))

        if length(neighbors) < 2
            curvatures[v] = 0.0  # Endpoint - zero curvature
        elseif length(neighbors) == 2
            # Calculate angle between neighbors
            v1, v2 = neighbors[1], neighbors[2]

            # Vectors from current vertex to neighbors
            vec1 = [track.voxels.x[v1] - track.voxels.x[v],
                   track.voxels.y[v1] - track.voxels.y[v],
                   track.voxels.z[v1] - track.voxels.z[v]]

            vec2 = [track.voxels.x[v2] - track.voxels.x[v],
                   track.voxels.y[v2] - track.voxels.y[v],
                   track.voxels.z[v2] - track.voxels.z[v]]

            # Normalize vectors
            norm1 = sqrt(sum(vec1.^2))
            norm2 = sqrt(sum(vec2.^2))

            if norm1 > 1e-6 && norm2 > 1e-6
                vec1 ./= norm1
                vec2 ./= norm2

                # Dot product gives cos(angle)
                cos_angle = dot(vec1, vec2)
                cos_angle = clamp(cos_angle, -1.0, 1.0)

                # Curvature = 1 - |cos(angle)| (0 for straight, 1 for sharp turn)
                curvatures[v] = 1.0 - abs(cos_angle)
            else
                curvatures[v] = 0.0
            end
        else
            # Multiple neighbors - use average curvature
            total_curvature = 0.0
            count = 0

            for i in 1:length(neighbors)-1
                for j in i+1:length(neighbors)
                    v1, v2 = neighbors[i], neighbors[j]

                    vec1 = [track.voxels.x[v1] - track.voxels.x[v],
                           track.voxels.y[v1] - track.voxels.y[v],
                           track.voxels.z[v1] - track.voxels.z[v]]

                    vec2 = [track.voxels.x[v2] - track.voxels.x[v],
                           track.voxels.y[v2] - track.voxels.y[v],
                           track.voxels.z[v2] - track.voxels.z[v]]

                    norm1 = sqrt(sum(vec1.^2))
                    norm2 = sqrt(sum(vec2.^2))

                    if norm1 > 1e-6 && norm2 > 1e-6
                        vec1 ./= norm1
                        vec2 ./= norm2
                        cos_angle = clamp(dot(vec1, vec2), -1.0, 1.0)
                        total_curvature += 1.0 - abs(cos_angle)
                        count += 1
                    end
                end
            end

            curvatures[v] = count > 0 ? total_curvature / count : 1.0  # High curvature for junctions
        end
    end

    return curvatures
end


"""Calculate total Euclidean length along a path."""
function calculate_path_length(track::Tracks, path::Vector{Int})
    if length(path) < 2
        return 0.0
    end

    total_length = 0.0
    for i in 1:length(path)-1
        v1, v2 = path[i], path[i+1]
        dist = euclidean_distance(
            track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
            track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
        )
        total_length += dist
    end

    return total_length
end


"""Check if two vertices are too close to be distinct extremes."""
function are_vertices_too_close(track::Tracks, v1::Int, v2::Int, threshold::Float64=2.0)
    dist = euclidean_distance(
        track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
        track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
    )
    return dist < threshold
end


"""Fallback: find extremes by maximum Euclidean distance."""
function find_extremes_distance_fallback(track::Tracks, confidence::Float64)
    g = track.graph
    vertices_to_check = collect(vertices(g))

    max_dist = 0.0
    best_pair = (vertices_to_check[1], vertices_to_check[min(2, end)])

    for i in 1:length(vertices_to_check)-1
        for j in i+1:length(vertices_to_check)
            v1, v2 = vertices_to_check[i], vertices_to_check[j]
            dist = euclidean_distance(
                track.voxels.x[v1], track.voxels.y[v1], track.voxels.z[v1],
                track.voxels.x[v2], track.voxels.y[v2], track.voxels.z[v2]
            )
            if dist > max_dist
                max_dist = dist
                best_pair = (v1, v2)
            end
        end
    end

    path = find_path_bfs(g, best_pair[1], best_pair[2])
    return (best_pair[1], best_pair[2], path, confidence)
end


# =============================================================================
# OPTIMIZED TRACK EXTREMES ALGORITHMS
# =============================================================================

"""
    TrackCoords

Pre-extracted coordinates for faster access (avoids DataFrame column overhead).
"""
struct TrackCoords
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

"""Extract coordinates from track voxels once."""
function extract_coords(track::Tracks)
    return TrackCoords(
        collect(track.voxels.x),
        collect(track.voxels.y),
        collect(track.voxels.z)
    )
end

"""
    find_path_bfs_opt(g, start_vertex, end_vertex)

Optimized BFS using parent pointers instead of copying paths.
Returns (path, path_length) where path_length is the number of edges.
"""
function find_path_bfs_opt(g::SimpleGraph, start_vertex::Int, end_vertex::Int)
    if start_vertex == end_vertex
        return [start_vertex]
    end

    n = nv(g)
    parent = zeros(Int, n)
    visited = falses(n)

    # Use a simple queue (circular buffer would be faster but this is clearer)
    queue = Int[start_vertex]
    visited[start_vertex] = true

    while !isempty(queue)
        current = popfirst!(queue)

        for neighbor in neighbors(g, current)
            if !visited[neighbor]
                visited[neighbor] = true
                parent[neighbor] = current

                if neighbor == end_vertex
                    # Reconstruct path from parent pointers
                    path = Int[end_vertex]
                    node = end_vertex
                    while parent[node] != 0
                        node = parent[node]
                        pushfirst!(path, node)
                    end
                    return path
                end

                push!(queue, neighbor)
            end
        end
    end

    return Int[]  # No path found
end

"""
    calculate_path_length_opt(coords, path)

Calculate path length using pre-extracted coordinates.
"""
function calculate_path_length_opt(coords::TrackCoords, path::Vector{Int})
    length(path) < 2 && return 0.0

    total = 0.0
    @inbounds for i in 1:length(path)-1
        v1, v2 = path[i], path[i+1]
        dx = coords.x[v2] - coords.x[v1]
        dy = coords.y[v2] - coords.y[v1]
        dz = coords.z[v2] - coords.z[v1]
        total += sqrt(dx*dx + dy*dy + dz*dz)
    end
    return total
end

"""
    euclidean_distance_opt(coords, v1, v2)

Euclidean distance using pre-extracted coordinates.
"""
@inline function euclidean_distance_opt(coords::TrackCoords, v1::Int, v2::Int)
    dx = coords.x[v2] - coords.x[v1]
    dy = coords.y[v2] - coords.y[v1]
    dz = coords.z[v2] - coords.z[v1]
    return sqrt(dx*dx + dy*dy + dz*dz)
end

"""
    find_track_extremes_opt(track)

Optimized version of find_track_extremes using pre-extracted coordinates
and early-exit strategies.
"""
function find_track_extremes_opt(track::Tracks)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0)
    elseif n_vertices == 1
        return (1, 1, [1], 1.0)
    end

    # Pre-extract coordinates once
    coords = extract_coords(track)

    return find_extremes_combined_opt(track, coords)
end

"""
    find_extremes_topology_opt(track, coords)

Optimized topology-based extreme finding.
Now includes path coverage check to avoid false positives from adjacent degree-1 vertices.
"""
function find_extremes_topology_opt(track::Tracks, coords::TrackCoords)
    g = track.graph
    n_vertices = nv(g)

    # Compute track extent for coverage check
    x_range = maximum(coords.x) - minimum(coords.x)
    y_range = maximum(coords.y) - minimum(coords.y)
    z_range = maximum(coords.z) - minimum(coords.z)
    track_extent = sqrt(x_range^2 + y_range^2 + z_range^2)

    # Find degree-1 vertices efficiently
    endpoints = Int[]
    sizehint!(endpoints, 10)
    for v in 1:n_vertices
        if degree(g, v) == 1
            push!(endpoints, v)
        end
    end

    # Exactly 2 endpoints - check if they give good coverage
    if length(endpoints) == 2
        path = find_path_bfs_opt(g, endpoints[1], endpoints[2])
        path_length = calculate_path_length_opt(coords, path)
        coverage = path_length / max(track_extent, 1.0)

        # High confidence only if coverage is good
        if coverage >= 0.5
            confidence = 0.85 + 0.1 * min(1.0, coverage)  # 0.85-0.95
        else
            confidence = 0.5 + 0.3 * coverage  # 0.5-0.65 for poor coverage
        end

        return (endpoints[1], endpoints[2], path, confidence, path_length)
    end

    # More than 2 endpoints - find pair with maximum path length
    if length(endpoints) > 2
        max_path_length = 0.0
        best_pair = (endpoints[1], endpoints[2])
        best_path = Int[]

        for i in 1:length(endpoints)-1
            for j in i+1:length(endpoints)
                path = find_path_bfs_opt(g, endpoints[i], endpoints[j])
                path_length = calculate_path_length_opt(coords, path)

                if path_length > max_path_length
                    max_path_length = path_length
                    best_pair = (endpoints[i], endpoints[j])
                    best_path = path
                end
            end
        end

        coverage = max_path_length / max(track_extent, 1.0)
        if coverage >= 0.5
            confidence = 0.75 + 0.15 * min(1.0, coverage)  # 0.75-0.90
        else
            confidence = 0.4 + 0.3 * coverage  # 0.4-0.55 for poor coverage
        end

        return (best_pair[1], best_pair[2], best_path, confidence, max_path_length)
    end

    # Fallback: no clear endpoints
    return find_extremes_distance_fallback_opt(track, coords, 0.6)
end

"""
    calculate_vertex_curvatures_opt(track, coords)

Optimized curvature calculation avoiding temporary array allocations.
"""
function calculate_vertex_curvatures_opt(track::Tracks, coords::TrackCoords)
    g = track.graph
    n_vertices = nv(g)
    curvatures = zeros(n_vertices)

    @inbounds for v in 1:n_vertices
        nbrs = neighbors(g, v)
        n_nbrs = length(nbrs)

        if n_nbrs < 2
            curvatures[v] = 0.0  # Endpoint
        elseif n_nbrs == 2
            v1, v2 = nbrs[1], nbrs[2]

            # Compute vectors inline (no array allocation)
            dx1 = coords.x[v1] - coords.x[v]
            dy1 = coords.y[v1] - coords.y[v]
            dz1 = coords.z[v1] - coords.z[v]

            dx2 = coords.x[v2] - coords.x[v]
            dy2 = coords.y[v2] - coords.y[v]
            dz2 = coords.z[v2] - coords.z[v]

            norm1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
            norm2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2)

            if norm1 > 1e-6 && norm2 > 1e-6
                # Normalize and compute dot product
                inv_norm1 = 1.0 / norm1
                inv_norm2 = 1.0 / norm2
                cos_angle = (dx1*dx2 + dy1*dy2 + dz1*dz2) * inv_norm1 * inv_norm2
                cos_angle = clamp(cos_angle, -1.0, 1.0)
                curvatures[v] = 1.0 - abs(cos_angle)
            end
        else
            # Multiple neighbors - compute average curvature
            total_curv = 0.0
            count = 0

            for i in 1:n_nbrs-1
                for j in i+1:n_nbrs
                    vi, vj = nbrs[i], nbrs[j]

                    dx1 = coords.x[vi] - coords.x[v]
                    dy1 = coords.y[vi] - coords.y[v]
                    dz1 = coords.z[vi] - coords.z[v]

                    dx2 = coords.x[vj] - coords.x[v]
                    dy2 = coords.y[vj] - coords.y[v]
                    dz2 = coords.z[vj] - coords.z[v]

                    norm1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1)
                    norm2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2)

                    if norm1 > 1e-6 && norm2 > 1e-6
                        inv_norm1 = 1.0 / norm1
                        inv_norm2 = 1.0 / norm2
                        cos_angle = (dx1*dx2 + dy1*dy2 + dz1*dz2) * inv_norm1 * inv_norm2
                        cos_angle = clamp(cos_angle, -1.0, 1.0)
                        total_curv += 1.0 - abs(cos_angle)
                        count += 1
                    end
                end
            end

            curvatures[v] = count > 0 ? total_curv / count : 1.0
        end
    end

    return curvatures
end

"""
    find_extremes_curvature_opt(track, coords)

Optimized curvature-based extreme finding.
"""
function find_extremes_curvature_opt(track::Tracks, coords::TrackCoords)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices < 3
        return find_extremes_distance_fallback_opt(track, coords, 0.7)
    end

    curvatures = calculate_vertex_curvatures_opt(track, coords)
    sorted_indices = sortperm(curvatures)

    max_score = 0.0
    best_pair = (sorted_indices[1], sorted_indices[2])
    best_path = Int[]
    best_path_length = 0.0
    confidence = 0.0

    n_candidates = min(6, n_vertices)

    for i in 1:n_candidates-1
        for j in i+1:n_candidates
            v1, v2 = sorted_indices[i], sorted_indices[j]

            # Skip if too close
            if euclidean_distance_opt(coords, v1, v2) < 2.0
                continue
            end

            path = find_path_bfs_opt(g, v1, v2)
            if !isempty(path)
                path_length = calculate_path_length_opt(coords, path)
                curvature_bonus = 2.0 / (curvatures[v1] + curvatures[v2] + 1e-6)
                score = path_length + curvature_bonus

                if score > max_score
                    max_score = score
                    best_pair = (v1, v2)
                    best_path = path
                    best_path_length = path_length
                    confidence = 0.7 + 0.2 * curvature_bonus / 10.0
                end
            end
        end
    end

    return (best_pair[1], best_pair[2], best_path, min(confidence, 0.9), best_path_length)
end

"""
    find_extremes_spatial_opt(track, coords)

Optimized spatial-based extreme finding for dense tracks.
"""
function find_extremes_spatial_opt(track::Tracks, coords::TrackCoords)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0, 0.0)
    elseif n_vertices == 1
        return (1, 1, [1], 1.0, 0.0)
    end

    # Spatial extreme candidates
    spatial_candidates = [
        argmin(coords.x), argmax(coords.x),
        argmin(coords.y), argmax(coords.y),
        argmin(coords.z), argmax(coords.z)
    ]

    # Add minimum degree vertices
    min_deg = minimum(degree(g, v) for v in 1:n_vertices)
    for v in 1:n_vertices
        if degree(g, v) == min_deg
            push!(spatial_candidates, v)
        end
    end

    all_candidates = unique(spatial_candidates)

    max_path_length = 0.0
    best_pair = (all_candidates[1], all_candidates[min(2, length(all_candidates))])
    best_path = Int[]

    for i in 1:length(all_candidates)-1
        for j in i+1:length(all_candidates)
            v1, v2 = all_candidates[i], all_candidates[j]
            path = find_path_bfs_opt(g, v1, v2)
            if !isempty(path)
                path_length = calculate_path_length_opt(coords, path)
                if path_length > max_path_length
                    max_path_length = path_length
                    best_pair = (v1, v2)
                    best_path = path
                end
            end
        end
    end

    # Calculate confidence
    confidence = 0.75
    if !isempty(best_path) && max_path_length > 0
        straight_dist = euclidean_distance_opt(coords, best_pair[1], best_pair[2])
        efficiency = straight_dist / max_path_length

        if efficiency > 0.8
            confidence += 0.1
        end
        confidence += min(0.1, max_path_length / 500.0)

        deg1 = degree(g, best_pair[1])
        deg2 = degree(g, best_pair[2])
        avg_endpoint_degree = (deg1 + deg2) / 2.0
        avg_track_degree = 2 * ne(g) / n_vertices

        if avg_endpoint_degree < avg_track_degree * 0.7
            confidence += 0.05
        end
        confidence = min(confidence, 0.95)
    end

    return (best_pair[1], best_pair[2], best_path, confidence, max_path_length)
end

"""
    find_extremes_combined_opt(track, coords; use_energy_weighting=true)

Optimized combined approach with early-exit strategies.
Key optimizations:
1. Pre-extracted coordinates passed through
2. Path lengths returned with results (no recomputation)
3. Early exit when topology gives high confidence
4. Energy-weighted method for dense tracks (Bragg peak detection)
"""
function find_extremes_combined_opt(track::Tracks, coords::TrackCoords; use_energy_weighting::Bool=true, use_edge_energy_weighting::Bool=true, use_mst_fallback::Bool=false)
    g = track.graph
    n_vertices = nv(g)

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0)
    end

    # Check track density
    avg_degree = 2 * ne(g) / n_vertices
    is_dense = avg_degree > 6.0

    # Always try topology first - it's fast and often sufficient
    topo_result = find_extremes_topology_opt(track, coords)
    topo_confidence = topo_result[4]
    topo_length = topo_result[5]

    # Early exit for high-confidence topology result (2 degree-1 vertices found)
    if topo_confidence >= 0.85
        return (topo_result[1], topo_result[2], topo_result[3], topo_confidence)
    end

    if is_dense
        # Dense track: use energy-weighted method (best for Bragg peak detection)

        # Dense track: energy-aware traversal can suppress geometric shortcuts
        if use_edge_energy_weighting
            edge_result = find_extremes_edge_energy_weighted_opt(track, coords)
            edge_confidence = edge_result[4]
            edge_length = edge_result[5]
            if edge_confidence >= 0.85 || edge_length > topo_length * 1.1
                return (edge_result[1], edge_result[2], edge_result[3], edge_confidence)
            end
        end

        if use_energy_weighting
            energy_result = find_extremes_energy_weighted_opt(track, coords)
            energy_confidence = energy_result[4]
            energy_length = energy_result[5]

            if use_mst_fallback
                mst_result = find_extremes_mst_diameter(track, coords)
                mst_length = mst_result[5]
                if mst_length >= max(spatial_length, curv_length, topo_length, energy_length)
                    return (mst_result[1], mst_result[2], mst_result[3], 1.0)
                end
            end


            # Energy method is preferred for dense tracks
            if energy_confidence >= 0.85 || energy_length > topo_length * 1.1
                return (energy_result[1], energy_result[2], energy_result[3], energy_confidence)
            end
        end

        # Fallback to spatial method
        spatial_result = find_extremes_spatial_opt(track, coords)
        spatial_length = spatial_result[5]

        if spatial_length > topo_length * 1.2
            return (spatial_result[1], spatial_result[2], spatial_result[3], spatial_result[4])
        end

        # Try curvature as final fallback
        curv_result = find_extremes_curvature_opt(track, coords)
        curv_length = curv_result[5]

        # Pick the best by path length
        if use_energy_weighting
            energy_result = find_extremes_energy_weighted_opt(track, coords)
            energy_length = energy_result[5]

            if use_mst_fallback
                mst_result = find_extremes_mst_diameter(track, coords)
                mst_length = mst_result[5]
                if mst_length >= max(spatial_length, curv_length, topo_length, energy_length)
                    return (mst_result[1], mst_result[2], mst_result[3], 1.0)
                end
            end

            if energy_length >= spatial_length && energy_length >= curv_length && energy_length >= topo_length
                return (energy_result[1], energy_result[2], energy_result[3], energy_result[4])
            end
        end


        if use_mst_fallback
            mst_result = find_extremes_mst_diameter(track, coords)
            mst_length = mst_result[5]
            if mst_length >= max(spatial_length, curv_length, topo_length)
                return (mst_result[1], mst_result[2], mst_result[3], 1.0)
            end
        end

        if spatial_length >= topo_length && spatial_length >= curv_length
            return (spatial_result[1], spatial_result[2], spatial_result[3], spatial_result[4])
        elseif curv_length > topo_length
            return (curv_result[1], curv_result[2], curv_result[3], curv_result[4])
        else
            return (topo_result[1], topo_result[2], topo_result[3], topo_confidence)
        end
    else
        # Sparse track: topology + curvature comparison
        if topo_confidence >= 0.8
            return (topo_result[1], topo_result[2], topo_result[3], topo_confidence)
        end

        curv_result = find_extremes_curvature_opt(track, coords)
        curv_length = curv_result[5]

        if curv_length > topo_length * 1.1
            return (curv_result[1], curv_result[2], curv_result[3], curv_result[4])
        else
            return (topo_result[1], topo_result[2], topo_result[3], topo_confidence)
        end
    end
end

"""
    find_extremes_energy_weighted_opt(track, coords; min_coverage=0.6)

Energy-weighted extreme finding. Bragg peaks at track endpoints have high energy.
Combines spatial extent with energy to find true endpoints.

Confidence calibration based on:
1. Energy at endpoints (Bragg peak detection)
2. Path coverage (path_length / track_extent)
3. Endpoint separation quality
"""
function find_extremes_energy_weighted_opt(track::Tracks, coords::TrackCoords; min_coverage::Float64=0.6)
    g = track.graph
    n_vertices = nv(g)
    vox = track.voxels

    if n_vertices == 0
        return (nothing, nothing, Int[], 0.0, 0.0)
    elseif n_vertices == 1
        return (1, 1, [1], 1.0, 0.0)
    end

    # Get energy values
    energies = collect(vox.energy)
    mean_energy = mean(energies)
    max_energy = maximum(energies)

    # Track spatial extent (diagonal of bounding box)
    x_range = maximum(coords.x) - minimum(coords.x)
    y_range = maximum(coords.y) - minimum(coords.y)
    z_range = maximum(coords.z) - minimum(coords.z)
    track_extent = sqrt(x_range^2 + y_range^2 + z_range^2)

    # Find high-energy voxels (potential Bragg peaks)
    # Threshold: voxels with energy > mean + 0.5*(max-mean)
    energy_threshold = mean_energy + 0.5 * (max_energy - mean_energy)
    high_energy_indices = findall(e -> e >= energy_threshold, energies)

    # Also include spatial extremes
    spatial_candidates = [
        argmin(coords.x), argmax(coords.x),
        argmin(coords.y), argmax(coords.y),
        argmin(coords.z), argmax(coords.z)
    ]

    # Combine candidates: high energy + spatial extremes
    all_candidates = unique(vcat(high_energy_indices, spatial_candidates))

    # Limit candidates to avoid O(n²) explosion
    if length(all_candidates) > 20
        sorted_by_energy = sortperm(energies[high_energy_indices], rev=true)
        top_energy = high_energy_indices[sorted_by_energy[1:min(10, length(sorted_by_energy))]]
        all_candidates = unique(vcat(top_energy, spatial_candidates))
    end

    # Score each pair by: path_length * (1 + energy_factor) * coverage_bonus
    best_score = -Inf
    best_pair = (all_candidates[1], all_candidates[min(2, length(all_candidates))])
    best_path = Int[]
    best_path_length = 0.0
    best_coverage = 0.0

    for i in 1:length(all_candidates)-1
        for j in i+1:length(all_candidates)
            v1, v2 = all_candidates[i], all_candidates[j]

            # Skip if too close
            dist = euclidean_distance_opt(coords, v1, v2)
            if dist < 5.0  # mm, minimum separation
                continue
            end

            path = find_path_bfs_opt(g, v1, v2)
            if isempty(path)
                continue
            end

            path_length = calculate_path_length_opt(coords, path)

            # Coverage: how much of the track extent does this path cover?
            coverage = path_length / max(track_extent, 1.0)

            # Skip pairs with poor coverage
            if coverage < min_coverage
                continue
            end

            # Energy bonus: normalized energy at both endpoints
            e1_norm = energies[v1] / max_energy
            e2_norm = energies[v2] / max_energy
            energy_factor = 0.5 * (e1_norm + e2_norm)

            # Coverage bonus: reward paths that cover more of the track
            coverage_bonus = min(coverage, 1.5)  # cap at 1.5

            # Combined score
            score = path_length * (1.0 + energy_factor) * coverage_bonus

            if score > best_score
                best_score = score
                best_pair = (v1, v2)
                best_path = path
                best_path_length = path_length
                best_coverage = coverage
            end
        end
    end

    # If no valid pair found with min_coverage, relax and try again
    if best_score == -Inf
        return find_extremes_energy_weighted_relaxed(track, coords, energies, energy_threshold,
                                                      max_energy, track_extent, all_candidates)
    end

    # Calibrated confidence based on multiple factors
    e1 = energies[best_pair[1]]
    e2 = energies[best_pair[2]]
    both_high_energy = (e1 >= energy_threshold) && (e2 >= energy_threshold)
    one_high_energy = (e1 >= energy_threshold) || (e2 >= energy_threshold)

    # Base confidence from energy
    if both_high_energy
        conf_energy = 0.4
    elseif one_high_energy
        conf_energy = 0.25
    else
        conf_energy = 0.1
    end

    # Confidence from coverage (0.0 to 0.4)
    conf_coverage = min(0.4, 0.4 * best_coverage / 1.0)

    # Confidence from path quality (0.0 to 0.2)
    # Straight-line efficiency: how direct is the path?
    straight_dist = euclidean_distance_opt(coords, best_pair[1], best_pair[2])
    efficiency = straight_dist / max(best_path_length, 1.0)
    conf_efficiency = 0.2 * min(1.0, efficiency / 0.8)  # max at 80% efficiency

    # Total confidence (0.0 to 1.0)
    confidence = min(0.95, conf_energy + conf_coverage + conf_efficiency)

    return (best_pair[1], best_pair[2], best_path, confidence, best_path_length)
end


"""
Fallback when min_coverage cannot be satisfied.
Returns results with LOW confidence to flag unreliable extremes.
"""
function find_extremes_energy_weighted_relaxed(track::Tracks, coords::TrackCoords,
                                                energies::Vector{Float64}, energy_threshold::Float64,
                                                max_energy::Float64, track_extent::Float64,
                                                all_candidates::Vector{Int})
    g = track.graph

    best_score = -Inf
    best_pair = (all_candidates[1], all_candidates[min(2, length(all_candidates))])
    best_path = Int[]
    best_path_length = 0.0

    for i in 1:length(all_candidates)-1
        for j in i+1:length(all_candidates)
            v1, v2 = all_candidates[i], all_candidates[j]

            dist = euclidean_distance_opt(coords, v1, v2)
            if dist < 3.0
                continue
            end

            path = find_path_bfs_opt(g, v1, v2)
            if isempty(path)
                continue
            end

            path_length = calculate_path_length_opt(coords, path)
            e1_norm = energies[v1] / max_energy
            e2_norm = energies[v2] / max_energy
            energy_factor = 0.5 * (e1_norm + e2_norm)
            score = path_length * (1.0 + energy_factor)

            if score > best_score
                best_score = score
                best_pair = (v1, v2)
                best_path = path
                best_path_length = path_length
            end
        end
    end

    # LOW confidence for relaxed results - coverage-dependent
    # This flags events where we couldn't find good extremes
    coverage = best_path_length / max(track_extent, 1.0)
    confidence = min(0.5, 0.2 + 0.3 * coverage)  # max 0.5 for relaxed

    return (best_pair[1], best_pair[2], best_path, confidence, best_path_length)
end


"""
    find_extremes_distance_fallback_opt(track, coords, confidence)

Optimized fallback using maximum Euclidean distance.
"""
function find_extremes_distance_fallback_opt(track::Tracks, coords::TrackCoords, confidence::Float64)
    g = track.graph
    n_vertices = nv(g)

    max_dist = 0.0
    best_pair = (1, min(2, n_vertices))

    # For small tracks, check all pairs
    # For larger tracks, use spatial extremes as candidates
    if n_vertices <= 50
        for i in 1:n_vertices-1
            for j in i+1:n_vertices
                dist = euclidean_distance_opt(coords, i, j)
                if dist > max_dist
                    max_dist = dist
                    best_pair = (i, j)
                end
            end
        end
    else
        # Use spatial extremes for large tracks
        candidates = unique([
            argmin(coords.x), argmax(coords.x),
            argmin(coords.y), argmax(coords.y),
            argmin(coords.z), argmax(coords.z)
        ])

        for i in 1:length(candidates)-1
            for j in i+1:length(candidates)
                v1, v2 = candidates[i], candidates[j]
                dist = euclidean_distance_opt(coords, v1, v2)
                if dist > max_dist
                    max_dist = dist
                    best_pair = (v1, v2)
                end
            end
        end
    end

    path = find_path_bfs_opt(g, best_pair[1], best_pair[2])
    path_length = calculate_path_length_opt(coords, path)
    return (best_pair[1], best_pair[2], path, confidence, path_length)
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


############################
# Energy-weighted traversal #
############################

"""
    energy_weight_matrix(g, coords, energy; epsE=1e-6, α=1.0, β=1.0) -> SparseMatrixCSC{Float64,Int}

Build a sparse weight matrix for `dijkstra_shortest_paths` on an unweighted graph `g`.

Edge cost model (default):
    w_ij = (d_ij^α) / ((0.5*(E_i + E_j) + epsE)^β)

- `coords` must be 3×N (Float64) with columns aligned to graph vertices.
- `energy` is a length-N vector (Float64), typically `track.voxels.energy` or `track.voxels.electrons`.
- `epsE` prevents blow-up when energies are small.
- `α` controls geometric emphasis (α=1 is linear distance).
- `β` controls energy emphasis (β=1 is inverse-energy).
"""
function energy_weight_matrix(g::SimpleGraph,
                              coords::AbstractMatrix{<:Real},
                              energy::AbstractVector{<:Real};
                              epsE::Float64 = 1e-6,
                              α::Float64 = 1.0,
                              β::Float64 = 1.0)

    n = nv(g)
    I = Int[]
    J = Int[]
    V = Float64[]

    @inbounds for e in edges(g)
        i = src(e); j = dst(e)

        dx = float(coords[1,i] - coords[1,j])
        dy = float(coords[2,i] - coords[2,j])
        dz = float(coords[3,i] - coords[3,j])
        d  = sqrt(dx*dx + dy*dy + dz*dz)

        Ej = 0.5*(float(energy[i]) + float(energy[j])) + epsE
        w  = (d^α) / (Ej^β)

        # symmetric
        push!(I,i); push!(J,j); push!(V,w)
        push!(I,j); push!(J,i); push!(V,w)
    end

    return sparse(I, J, V, n, n)
end


"""
    dijkstra_path(g, src, dst, W) -> Vector{Int}

Reconstruct a shortest path from `src` to `dst` using Dijkstra parents.
Returns an empty vector if `dst` is unreachable.
"""
function dijkstra_path(g::SimpleGraph, srcv::Int, dstv::Int, W)
    sp = dijkstra_shortest_paths(g, srcv, W)
    parents = sp.parents

    dstv > length(parents) && return Int[]
    isinf(sp.dists[dstv]) && return Int[]

    path = Int[dstv]
    v = dstv
    while v != srcv
        p = parents[v]
        (p == 0 || p == v) && return Int[]  # unreachable / broken parent chain
        push!(path, p)
        v = p
    end
    reverse!(path)
    return path
end


"""
    find_extremes_edge_energy_weighted_opt(track, coords; epsE=1e-6, α=1.0, β=1.0)
        -> (extreme1, extreme2, path, confidence, path_length)

Energy-weighted *edge-cost* extreme finding:
1) Build sparse weight matrix from geometry + voxel energies
2) Double-sweep using weighted Dijkstra distances
3) Return endpoints and the weighted-shortest path between them

This is specifically designed to suppress "geometric shortcuts" that run through
low-support (low-energy) regions while preserving the connected manifold.
"""
function find_extremes_edge_energy_weighted_opt(track::Tracks, coords::TrackCoords;
                                               epsE::Float64 = 1e-6,
                                               α::Float64 = 1.0,
                                               β::Float64 = 1.0)

    g = track.graph
    n = nv(g)
    n == 0 && return (nothing, nothing, Int[], 0.0, 0.0)
    n == 1 && return (1, 1, [1], 1.0, 0.0)

    # 3×N coordinate matrix aligned to graph vertices
    coord_matrix = hcat(coords.x, coords.y, coords.z)'

    # energy vector aligned to vertices (use energy by default)
    E = Float64.(track.voxels.energy)

    # robust epsilon scaled to typical energy
    epsE_eff = max(epsE, 1e-6 * max(maximum(E), 1e-12))

    W = energy_weight_matrix(g, coord_matrix, E; epsE=epsE_eff, α=α, β=β)

    # Double sweep in weighted metric
    d1 = dijkstra_shortest_paths(g, 1, W).dists
    u = argmax(d1)
    d2 = dijkstra_shortest_paths(g, u, W).dists
    v = argmax(d2)

    path = dijkstra_path(g, u, v, W)
    isempty(path) && return (u, v, Int[], 0.0, 0.0)

    # Physical (Euclidean) path length along the returned vertex sequence
    path_length = calculate_path_length_from_coords(coord_matrix, path)

    # Confidence: combine coverage proxy + endpoint separation proxy
    # (keep it simple and diagnostic-driven; you can refine later)
    D_end = euclidean_distance(coords.x[u], coords.y[u], coords.z[u],
                               coords.x[v], coords.y[v], coords.z[v])
    η = (D_end > 0) ? (path_length / D_end) : Inf
    # η>1 indicates non-trivial curvature; extremely small η indicates shortcutting
    confidence = clamp(0.5 + 0.5 * tanh((η - 1.05) / 0.25), 0.0, 1.0)

    return (u, v, path, confidence, path_length)
end


#################
# MST + diameter #
#################

"""
    compute_mst_graph(g, coords) -> SimpleGraph

Compute an MST (minimum spanning tree) as an unweighted `SimpleGraph` using
Euclidean distances for edges already present in `g`.

Note: MST can be used as a fallback/topology-cleaner for pathological cases.
"""
function compute_mst_graph(g::SimpleGraph, coords::AbstractMatrix{<:Real})
    n = nv(g)
    (n <= 1 || ne(g) == 0) && return g

    # Build a dense distance matrix only for existing edges; others are ignored by `kruskal_mst`.
    distmx = fill(Inf, n, n)
    @inbounds for e in edges(g)
        i = src(e); j = dst(e)
        dx = float(coords[1,i] - coords[1,j])
        dy = float(coords[2,i] - coords[2,j])
        dz = float(coords[3,i] - coords[3,j])
        d  = sqrt(dx*dx + dy*dy + dz*dz)
        distmx[i,j] = d
        distmx[j,i] = d
    end

    mst_edges = kruskal_mst(g, distmx)  # returns an edge iterator/collection
    mst = SimpleGraph(n)
    for e in mst_edges
        add_edge!(mst, src(e), dst(e))
    end
    return mst
end


"""
    find_extremes_mst_diameter(track, coords) -> (extreme1, extreme2, path, confidence, path_length)

Fallback method:
1) Build MST from `track.graph`
2) Find diameter endpoints by double BFS (unweighted)
3) Return BFS path and Euclidean path length
"""
function find_extremes_mst_diameter(track::Tracks, coords::TrackCoords)
    g = track.graph
    n = nv(g)
    n == 0 && return (nothing, nothing, Int[], 0.0, 0.0)
    n == 1 && return (1, 1, [1], 1.0, 0.0)

    coord_matrix = hcat(coords.x, coords.y, coords.z)'
    mst = compute_mst_graph(g, coord_matrix)

    d1 = gdistances(mst, 1)
    u = argmax(d1)
    d2 = gdistances(mst, u)
    v = argmax(d2)

    path = find_path_bfs_opt(mst, u, v)
    path_length = calculate_path_length_from_coords(coord_matrix, path)

    return (u, v, path, 1.0, path_length)
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

    # Build KDTree over skeleton points
    Ns = length(path)
    skel = Matrix{Float64}(undef, 3, Ns)
    @inbounds for (k, idx) in enumerate(path)
        skel[1,k] = coords.x[idx]
        skel[2,k] = coords.y[idx]
        skel[3,k] = coords.z[idx]
    end
    tree = KDTree(skel)

    Ein = 0.0
    @inbounds for i in 1:n
        idxs, dists = knn(tree, [coords.x[i], coords.y[i], coords.z[i]], 1, true)
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

