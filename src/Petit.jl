module Petit

using DataFrames
using Graphs
using StatsBase
using Plots

include("histos.jl")



struct VGraph
    voxels::DataFrame
    graph::SimpleGraph{Int}
    components::Vector{Vector{Int}}
end


function voxelize_hits(ghdf::GroupedDataFrame, voxel_size::Float64)
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


function euclidean_distance(x1, y1, z1, x2, y2, z2)
    return sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
end


"""
    build_vgraph(gvdf::GroupedDataFrame, event_id::Int; max_distance::Float64=1.5, energy_threshold::Float64=0.0)

Build connected graphs of voxels from a grouped DataFrame using optimized graph algorithms.

# Algorithm Details

This function implements a voxel clustering algorithm based on spatial proximity and energy thresholding:

1. **Event Selection**: Extracts voxel data for the specified `event_id` from the grouped DataFrame
2. **Energy Filtering**: Filters out voxels with energy below `energy_threshold` to reduce noise
3. **Adjacency Calculation**: Computes Euclidean distances between all pairs of filtered voxels
4. **Graph Construction**: Creates edges between voxels whose distance ≤ `max_distance`
5. **Connected Components**: Uses Graphs.jl's optimized O(|V|) algorithm to find all connected components
6. **Subgraph Creation**: Generates individual subgraphs for each connected component with proper index remapping

# Arguments
- `gvdf::GroupedDataFrame`: DataFrame grouped by event_id, containing columns (event_id, x, y, z, energy)
- `event_id::Int`: Specific event to process
- `max_distance::Float64=1.5`: Maximum Euclidean distance for voxels to be considered adjacent
- `energy_threshold::Float64=0.0`: Minimum energy threshold for voxels to be included

# Returns
- `Vector{VGraph}`: Array of connected voxel graphs, each containing:
  - `voxels`: DataFrame with voxel data for the component
  - `graph`: Graphs.jl SimpleGraph object with adjacency information
  - `components`: Vector containing the component indices

# Performance
- Time complexity: O(n²) for distance calculations + O(|V|) for connected components
- Space complexity: O(n²) for adjacency matrix storage
- Optimized using Graphs.jl's high-performance connected components algorithm

# Example
```julia
# Build connected voxel graphs for event 42 with 1.5mm adjacency threshold
graphs = build_vgraph(grouped_voxel_df, 42; max_distance=1.5, energy_threshold=0.1)

# Access first component's data and graph structure
first_component = graphs[1]
voxel_data = first_component.voxels
graph_structure = first_component.graph
```
"""
function build_vgraph(gvdf::GroupedDataFrame, event_id::Int; 
                     max_distance::Float64=1.5, 
                     energy_threshold::Float64=0.0)
    
    # Get the specific event data
    event_data = nothing
    for group in gvdf
        if first(group.event_id) == event_id
            event_data = group
            break
        end
    end
    
    if event_data === nothing
        error("Event ID $event_id not found in grouped dataframe")
    end
    
    n_voxels = nrow(event_data)
    if n_voxels == 0
        return VGraph[]
    end
    
    # Filter voxels by energy threshold
    valid_indices = findall(e -> e >= energy_threshold, event_data.energy)
    if isempty(valid_indices)
        return VGraph[]
    end
    
    filtered_data = event_data[valid_indices, :]
    n_filtered = length(valid_indices)
    
    # Create graph with filtered voxels
    g = SimpleGraph(n_filtered)
    
    # Add edges based on distance
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
    
    # Find connected components using Graphs.jl optimized algorithm
    components = connected_components(g)
    
    # Create VGraph objects for each component
    graphs = VGraph[]
    
    for component in components
        if !isempty(component)
            component_data = filtered_data[component, :]
            
            # Create subgraph for this component
            subgraph_vertices = length(component)
            subgraph = SimpleGraph(subgraph_vertices)
            
            # Remap edges to new indices
            vertex_map = Dict(old_v => new_v for (new_v, old_v) in enumerate(component))
            
            for edge in edges(g)
                src_old, dst_old = src(edge), dst(edge)
                if src_old in component && dst_old in component
                    src_new = vertex_map[src_old]
                    dst_new = vertex_map[dst_old]
                    add_edge!(subgraph, src_new, dst_new)
                end
            end
            
            push!(graphs, VGraph(component_data, subgraph, [component]))
        end
    end
    
    return graphs
end

"""
    histogram_voxel_distances(gvdf::GroupedDataFrame; bins=50, max_distance=Inf)

Compute and histogram all pairwise distances between voxels for all events in the grouped DataFrame.

This function calculates the Euclidean distance between every pair of voxels within each event,
collecting all distances across all events and creating a histogram distribution.

# Arguments
- `gvdf::GroupedDataFrame`: DataFrame grouped by event_id, containing columns (event_id, x, y, z, energy)
- `bins::Int=50`: Number of histogram bins
- `max_distance::Float64=Inf`: Maximum distance to include in histogram (for filtering outliers)

# Returns
- `Histogram`: StatsBase histogram object containing the distance distribution

# Algorithm
1. Iterates through each event group in the grouped DataFrame
2. For each event, computes all pairwise Euclidean distances between voxels
3. Collects distances from all events into a single vector
4. Creates histogram with specified number of bins

# Performance
- Time complexity: O(∑ᵢ nᵢ²) where nᵢ is the number of voxels in event i
- Space complexity: O(∑ᵢ nᵢ²) for storing all distance pairs

# Example
```julia
# Histogram all voxel distances across all events
hist = histogram_voxel_distances(grouped_voxel_df; bins=100, max_distance=10.0)

# Access histogram data
distances = hist.edges[1]  # bin edges
counts = hist.weights      # bin counts
```
"""
function histogram_voxel_distances(gvdf::GroupedDataFrame; bins=50, maxevts=100, max_distance=Inf)
    all_distances = Float64[]
    
    for (i, group) in enumerate(gvdf)
		if i > maxevts 
			break
		end
        
        n_voxels = nrow(group)
        if n_voxels < 2
            continue  # Skip events with less than 2 voxels
        end
        
        # Compute all pairwise distances within this event
        for i in 1:(n_voxels-1)
            for j in (i+1):n_voxels
                dist = euclidean_distance(
                    group.x[i], group.y[i], group.z[i],
                    group.x[j], group.y[j], group.z[j]
                )
                
                # Only include distances within the specified maximum
                if dist <= max_distance
                    push!(all_distances, dist)
                end
            end
        end
    end
    
    # Create histogram using StatsBase
    if isempty(all_distances)
        # Return empty histogram if no distances found
        return fit(Histogram, Float64[], nbins=bins)
    else
        return histogram(all_distances, bins=bins, xlabel="distance", 
						 ylabel="Count", title="Distances", legend=false, color=:red, alpha=0.7)
    end
end

"""
    histogram_closest_distance(gvdf::GroupedDataFrame; bins=50, max_distance=Inf)

Compute and histogram the distance to the closest neighbor for every voxel across all events.

For each voxel in each event, this function finds the minimum distance to any other voxel 
within the same event, collecting all these minimum distances and creating a histogram.

# Arguments
- `gvdf::GroupedDataFrame`: DataFrame grouped by event_id, containing columns (event_id, x, y, z, energy)
- `bins::Int=50`: Number of histogram bins
- `max_distance::Float64=Inf`: Maximum distance to include in histogram (for filtering outliers)

# Returns
- `Plot`: Plots.jl histogram plot of the closest neighbor distance distribution

# Algorithm
1. For each event group in the grouped DataFrame:
   - For each voxel, compute distances to all other voxels in the same event
   - Find the minimum distance (closest neighbor)
   - Add this minimum distance to the collection
2. Create histogram of all closest neighbor distances

# Performance
- Time complexity: O(∑ᵢ nᵢ²) where nᵢ is the number of voxels in event i
- Space complexity: O(∑ᵢ nᵢ) for storing closest distances

# Example
```julia
# Histogram closest neighbor distances with 100 bins
hist = histogram_closest_distance(grouped_voxel_df; bins=100, max_distance=5.0)

# Access histogram data
distances = hist.edges[1]  # bin edges
counts = hist.weights      # bin counts
```
"""
function histogram_closest_distance(gvdf::GroupedDataFrame; bins=50, maxevts=100, max_distance=Inf)
    closest_distances = Float64[]
    
    for (i, group) in enumerate(gvdf)
		if i > maxevts 
			break
		end

        n_voxels = nrow(group)
        if n_voxels < 2
            continue  # Skip events with less than 2 voxels (no neighbors)
        end
        
        # For each voxel, find the closest neighbor
        for i in 1:n_voxels
            min_dist = Inf
            
            # Check distance to all other voxels in the same event
            for j in 1:n_voxels
                if i != j
                    dist = euclidean_distance(
                        group.x[i], group.y[i], group.z[i],
                        group.x[j], group.y[j], group.z[j]
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
    end
    
    # Create histogram using Plots
    if isempty(closest_distances)
        # Return empty histogram if no distances found
        return histogram(Float64[], bins=bins, xlabel="Distance", ylabel="Count", 
                        title="Closest Neighbor Distances", legend=false, color=:green, alpha=0.7)
    else
        return histogram(closest_distances, bins=bins, xlabel="Distance", ylabel="Count", 
                        title="Closest Neighbor Distances", legend=false, color=:green, alpha=0.7)
    end
end


function histogram_voxel_energy(gvdf::GroupedDataFrame; bins=50, title="Voxel Energy Distribution")
    all_energies = Float64[]
    
    for group in gvdf
        append!(all_energies, group.energy)
    end
    
    return histogram(all_energies, bins=bins, xlabel="Energy", ylabel="Count", 
                    title=title, legend=false, color=:blue, alpha=0.7)
end


export VGraph, voxelize_hits, build_vgraph, histogram_voxel_energy, histogram_voxel_distances, histogram_closest_distance

end
