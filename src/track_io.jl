"""
Functions for saving and loading Tracks objects to/from disk.
"""

using CSV
using JSON
using Serialization
using Glob

function save_track_csv(track::Tracks, filename::String)
    """
    Save a track to CSV files (voxels and graph structure).
    Creates two files: filename_voxels.csv and filename_graph.json

    Parameters:
    - track: Tracks object to save
    - filename: Base filename (without extension)
    """
    # Save voxels as CSV
    voxels_file = filename * "_voxels.csv"
    CSV.write(voxels_file, track.voxels)

    # Save graph structure as JSON
    graph_file = filename * "_graph.json"

    # Convert graph to edge list
    edge_list = []
    for edge in edges(track.graph)
        push!(edge_list, [src(edge), dst(edge)])
    end

    graph_data = Dict(
        "n_vertices" => nv(track.graph),
        "edges" => edge_list,
        "components" => track.components
    )

    open(graph_file, "w") do f
        JSON.print(f, graph_data, 2)
    end

    println("Track saved to:")
    println("  Voxels: $voxels_file")
    println("  Graph: $graph_file")
end

function load_track_csv(filename::String)
    """
    Load a track from CSV files.

    Parameters:
    - filename: Base filename (without extension)

    Returns:
    - Tracks object
    """
    # Load voxels
    voxels_file = filename * "_voxels.csv"
    voxels = CSV.read(voxels_file, DataFrame)

    # Load graph structure
    graph_file = filename * "_graph.json"
    graph_data = JSON.parsefile(graph_file)

    # Reconstruct graph
    n_vertices = graph_data["n_vertices"]
    g = SimpleGraph(n_vertices)

    for edge in graph_data["edges"]
        add_edge!(g, edge[1], edge[2])
    end

    components = graph_data["components"]

    return Tracks(voxels, g, components)
end

function save_track_binary(track::Tracks, filename::String)
    """
    Save a track using Julia's binary serialization (faster for large tracks).

    Parameters:
    - track: Tracks object to save
    - filename: Filename with .jls extension
    """
    if !endswith(filename, ".jls")
        filename = filename * ".jls"
    end

    serialize(filename, track)
    println("Track saved to: $filename")
end

function load_track_binary(filename::String)
    """
    Load a track from binary serialization.

    Parameters:
    - filename: Filename with .jls extension

    Returns:
    - Tracks object
    """
    if !endswith(filename, ".jls")
        filename = filename * ".jls"
    end

    return deserialize(filename)
end

function save_track_with_analysis(track::Tracks, filename::String; include_analysis::Bool=true)
    """
    Save a track with optional analysis results.

    Parameters:
    - track: Tracks object to save
    - filename: Base filename
    - include_analysis: Whether to include analysis results
    """
    # Save the track itself
    save_track_csv(track, filename)

    if include_analysis
        # Perform analysis
        println("Performing track analysis...")

        # Find extremes with improved algorithm
        extreme1, extreme2, path, confidence = find_track_extremes(track)

        # Walk the track
        walk_result = walk_track_from_extremes(track)

        # Calculate curvatures
        curvatures = calculate_vertex_curvatures(track)

        # Energy analysis (if we have extremes)
        sphere_results = if !isnothing(extreme1)
            energy_in_spheres_around_extremes(track, 2.0)  # 2mm radius
        else
            nothing
        end

        # Create analysis summary
        analysis = Dict(
            "track_info" => Dict(
                "n_voxels" => nrow(track.voxels),
                "n_vertices" => nv(track.graph),
                "n_edges" => ne(track.graph),
                "total_energy" => sum(track.voxels.energy)
            ),
            "extremes" => Dict(
                "extreme1" => extreme1,
                "extreme2" => extreme2,
                "confidence" => confidence,
                "path_length" => length(path)
            ),
            "walk_analysis" => Dict(
                "total_length_mm" => walk_result.total_length,
                "confidence" => walk_result.confidence
            ),
            "curvature_analysis" => Dict(
                "vertex_curvatures" => curvatures,
                "max_curvature" => maximum(curvatures),
                "mean_curvature" => mean(curvatures),
                "sharp_turn_vertices" => findall(c -> c > 0.5, curvatures)
            ),
            "sphere_analysis" => sphere_results !== nothing ? Dict(
                "radius_mm" => 2.0,
                "blob1_energy" => sphere_results.blob1_energy,
                "blob2_energy" => sphere_results.blob2_energy,
                "blob1_center" => sphere_results.blob1_center,
                "blob2_center" => sphere_results.blob2_center
            ) : nothing
        )

        # Save analysis
        analysis_file = filename * "_analysis.json"
        open(analysis_file, "w") do f
            JSON.print(f, analysis, 2)
        end

        println("Analysis saved to: $analysis_file")
    end
end

function print_track_summary(track::Tracks)
    """
    Print a summary of the track for inspection.
    """
    println("=" ^ 50)
    println("TRACK SUMMARY")
    println("=" ^ 50)

    println("Voxels: $(nrow(track.voxels))")
    println("Graph vertices: $(nv(track.graph))")
    println("Graph edges: $(ne(track.graph))")
    println("Components: $(length(track.components))")
    println("Total energy: $(round(sum(track.voxels.energy), digits=3))")

    # Show vertex degrees
    println("\nVertex degrees:")
    endpoints = Int[]
    for v in vertices(track.graph)
        deg = degree(track.graph, v)
        if deg == 1
            push!(endpoints, v)
        end
        println("  Vertex $v: degree $deg")
    end

    if !isempty(endpoints)
        println("Endpoints (degree-1): $endpoints")
    else
        println("No endpoints found (circular or complex topology)")
    end

    # Run analysis
    println("\n" * "-" ^ 30)
    println("TRACK ANALYSIS")
    println("-" ^ 30)

    extreme1, extreme2, path, confidence = find_track_extremes(track)
    println("Extremes: $extreme1 ↔ $extreme2 (confidence: $(round(confidence, digits=3)))")

    if !isnothing(extreme1)
        pos1 = (track.voxels.x[extreme1], track.voxels.y[extreme1], track.voxels.z[extreme1])
        pos2 = (track.voxels.x[extreme2], track.voxels.y[extreme2], track.voxels.z[extreme2])
        println("Positions: $pos1 ↔ $pos2")

        # Path analysis
        path_length = calculate_path_length(track, path)
        straight_dist = euclidean_distance(pos1[1], pos1[2], pos1[3], pos2[1], pos2[2], pos2[3])
        efficiency = straight_dist / path_length

        println("Path length: $(round(path_length, digits=2)) mm")
        println("Straight distance: $(round(straight_dist, digits=2)) mm")
        println("Path efficiency: $(round(efficiency, digits=3))")
    end

    # Curvature analysis
    curvatures = calculate_vertex_curvatures(track)
    if !isempty(curvatures)
        max_curv = maximum(curvatures)
        sharp_turns = findall(c -> c > 0.5, curvatures)

        println("Max curvature: $(round(max_curv, digits=3))")
        if !isempty(sharp_turns)
            println("Sharp turns at vertices: $sharp_turns")
        end
    end

    println("=" ^ 50)
end

"""
    save_tracks_to_hdf5(tracks::Vector{Tracks}, hdf5_file, batch_id::Int)

Save a batch of tracks to an HDF5 file.

# Arguments
- `tracks::Vector{Tracks}`: Vector of track objects to save
- `hdf5_file`: Open HDF5 file handle
- `batch_id::Int`: Batch identifier for organizing data

# Returns
- `Float64`: Elapsed time in seconds for saving
"""
function save_tracks_to_hdf5(tracks::Vector{Tracks}, hdf5_file, batch_id::Int)
    if isempty(tracks)
        return 0.0  # Return zero time
    end

    n_tracks = length(tracks)
    print("    Saving $n_tracks tracks to HDF5...")
    flush(stdout)

    start_time = time()
    for (idx, track) in enumerate(tracks)
        if idx % 10 == 0 || idx == n_tracks
            print("\r    Saving track $idx/$n_tracks to HDF5...")
            flush(stdout)
        end

        track_group = "batch_$(batch_id)/track_$(idx)"

        # Create group if it doesn't exist
        if !haskey(hdf5_file, track_group)
            g = create_group(hdf5_file, track_group)
        else
            g = hdf5_file[track_group]
        end

        # Save voxels data (this is the slow part for large tracks)
        voxels_data = Matrix(track.voxels)
        g["voxels"] = voxels_data
        g["voxel_columns"] = String.(names(track.voxels))

        # Save graph edges - optimize by preallocating
        n_edges = ne(track.graph)
        if n_edges > 0
            edge_matrix = zeros(Int, n_edges, 2)
            for (i, edge) in enumerate(edges(track.graph))
                edge_matrix[i, 1] = src(edge)
                edge_matrix[i, 2] = dst(edge)
            end
            g["graph_edges"] = edge_matrix
        else
            g["graph_edges"] = zeros(Int, 0, 2)
        end

        # Save graph metadata
        g["n_vertices"] = nv(track.graph)

        # Save components
        if !isempty(track.components)
            max_comp_len = maximum(length.(track.components))
            comp_matrix = zeros(Int, length(track.components), max_comp_len)
            for (i, comp) in enumerate(track.components)
                comp_matrix[i, 1:length(comp)] .= comp
            end
            g["components"] = comp_matrix
            g["component_lengths"] = length.(track.components)
        end

        # Save diffusion parameters
        dp = track.diffusion
        g["diffusion_ldrift"] = dp.ldrift
        g["diffusion_sigma_t"] = dp.sigma_t
        g["diffusion_sigma_l"] = dp.sigma_l
        g["diffusion_voxel_size"] = dp.voxel_size
        g["diffusion_max_distance"] = dp.max_distance
        g["diffusion_energy_threshold"] = dp.energy_threshold
        g["diffusion_nbins_df"] = dp.nbins_df
        g["diffusion_nsigma_df"] = dp.nsigma_df
    end
    elapsed_time = time() - start_time
    println("\r    Saved $n_tracks tracks to HDF5 in $(round(elapsed_time, digits=2))s        ")
    return elapsed_time
end

"""
    read_tracks_from_hdf5(output_file::String)

Read tracks from the output HDF5 file and return them as a vector.

# Arguments
- `output_file::String`: Path to the HDF5 output file

# Returns
- `Vector{Tracks}`: Vector of reconstructed track objects
- `Dict`: Metadata from the file
"""
function read_tracks_from_hdf5(output_file::String)
    println("Reading tracks from: $output_file")
    tracks = Tracks[]
    metadata = Dict{String, Any}()
    start_time = time()

    h5open(output_file, "r") do fid
        # Read metadata
        println("  Reading metadata...")
        for attr_name in keys(attrs(fid))
            metadata[attr_name] = read_attribute(fid, attr_name)
        end

        # Count total batches and tracks first
        batch_keys = filter(k -> startswith(k, "batch_"), collect(keys(fid)))
        total_batches = length(batch_keys)
        println("  Found $total_batches batches")

        # Read all batches
        for (batch_idx, batch_key) in enumerate(batch_keys)
            batch_group = fid[batch_key]
            track_keys = filter(k -> startswith(k, "track_"), collect(keys(batch_group)))
            n_tracks_in_batch = length(track_keys)

            print("  Reading batch $batch_idx/$total_batches ($n_tracks_in_batch tracks)...")
            flush(stdout)

            # Read all tracks in this batch
            for track_key in track_keys
                track_group = batch_group[track_key]

                # Reconstruct voxels DataFrame
                voxels_data = read(track_group["voxels"])
                voxels_columns = read(track_group["voxel_columns"])
                voxels = DataFrame(voxels_data, Symbol.(voxels_columns))

                # Reconstruct graph
                n_vertices = read(track_group["n_vertices"])
                g = SimpleGraph(n_vertices)

                graph_edges = read(track_group["graph_edges"])
                if size(graph_edges, 1) > 0
                    for i in 1:size(graph_edges, 1)
                        add_edge!(g, graph_edges[i, 1], graph_edges[i, 2])
                    end
                end

                # Reconstruct components
                comp_matrix = read(track_group["components"])
                comp_lengths = read(track_group["component_lengths"])
                components = Vector{Vector{Int}}(undef, length(comp_lengths))
                for i in 1:length(comp_lengths)
                    components[i] = Vector{Int}(comp_matrix[i, 1:comp_lengths[i]])
                end

                # Reconstruct diffusion parameters (with backwards compatibility)
                diffusion = if haskey(track_group, "diffusion_ldrift")
                    DiffusionParams(
                        read(track_group["diffusion_ldrift"]),
                        read(track_group["diffusion_sigma_t"]),
                        read(track_group["diffusion_sigma_l"]),
                        read(track_group["diffusion_voxel_size"]),
                        read(track_group["diffusion_max_distance"]),
                        read(track_group["diffusion_energy_threshold"]),
                        read(track_group["diffusion_nbins_df"]),
                        read(track_group["diffusion_nsigma_df"])
                    )
                else
                    DiffusionParams()  # Default for old files
                end

                # Create track object
                track = Tracks(voxels, g, components, diffusion)
                push!(tracks, track)
            end
            println(" done")
        end
    end

    elapsed_time = time() - start_time
    n_tracks = length(tracks)
    println("Loaded $n_tracks tracks from $output_file in $(round(elapsed_time, digits=2))s")
    if n_tracks > 0
        println("Average time per track (reading): $(round(1000*elapsed_time/n_tracks, digits=2))ms")
    end

    # Add timing info to metadata
    metadata["reading_time"] = elapsed_time
    metadata["avg_read_time_per_track_ms"] = n_tracks > 0 ? 1000*elapsed_time/n_tracks : 0.0

    return tracks, metadata
end

"""
    chain_track_files(file_paths::Vector{String})

Read and chain tracks from multiple HDF5 files.

# Arguments
- `file_paths::Vector{String}`: Vector of paths to HDF5 track files

# Returns
- `Vector{Tracks}`: Vector containing all tracks from all files
- `Vector{Dict}`: Vector of metadata dictionaries, one per file

# Example
```julia
files = ["file1.h5", "file2.h5", "file3.h5"]
all_tracks, metadatas = chain_track_files(files)
println("Total tracks: ", length(all_tracks))
```
"""
function chain_track_files(file_paths::Vector{String})
    all_tracks = Tracks[]
    all_metadata = Dict{String, Any}[]

    println("Chaining $(length(file_paths)) track files...")

    for (idx, file_path) in enumerate(file_paths)
        if !isfile(file_path)
            println("Warning: File not found, skipping: $file_path")
            continue
        end

        println("\n[$idx/$(length(file_paths))] Reading: $file_path")
        tracks, metadata = read_tracks_from_hdf5(file_path)

        append!(all_tracks, tracks)
        push!(all_metadata, metadata)

        println("  Added $(length(tracks)) tracks (total so far: $(length(all_tracks)))")
    end

    println("\n" * "="^60)
    println("CHAINING COMPLETE")
    println("="^60)
    println("Total files processed: $(length(all_metadata))")
    println("Total tracks: $(length(all_tracks))")

    return all_tracks, all_metadata
end

"""
    chain_track_files(file_paths::String...)

Read and chain tracks from multiple HDF5 files (variadic version).

# Arguments
- `file_paths::String...`: Paths to HDF5 track files as individual arguments

# Returns
- `Vector{Tracks}`: Vector containing all tracks from all files
- `Vector{Dict}`: Vector of metadata dictionaries, one per file

# Example
```julia
all_tracks, metadatas = chain_track_files("file1.h5", "file2.h5", "file3.h5")
```
"""
function chain_track_files(file_paths::String...)
    return chain_track_files(collect(file_paths))
end

"""
    merge_csv_files(directory::String; pattern::String="*.csv", sort_column::Union{Symbol,Nothing}=nothing)

Read and merge all CSV files from a directory into a single DataFrame.

# Arguments
- `directory::String`: Path to the directory containing CSV files
- `pattern::String`: Glob pattern for matching CSV files (default: "*.csv")
- `sort_column::Union{Symbol,Nothing}`: Optional column name to sort the merged DataFrame by (default: nothing)

# Returns
- `DataFrame`: Merged DataFrame containing all data from matching CSV files

# Example
```julia
# Merge all CSV files from multithreaded analysis
df = merge_csv_files("xev1mm")

# Merge with specific pattern and sort by confidence
df = merge_csv_files("xev1mm", pattern="xev1mm_th_*.csv", sort_column=:confidence)
```
"""
function merge_csv_files(directory::String; pattern::String="*.csv", sort_column::Union{Symbol,Nothing}=nothing)
    if !isdir(directory)
        error("Directory does not exist: $directory")
    end

    # Find all matching CSV files in the directory
    csv_files = Glob.glob(pattern, directory)

    if isempty(csv_files)
        error("No CSV files found matching pattern '$pattern' in directory: $directory")
    end

    println("Found $(length(csv_files)) CSV files to merge")

    # Read and merge all CSV files
    dfs = DataFrame[]

    for (idx, csv_file) in enumerate(csv_files)
        # Glob.glob returns full paths when given a directory, so use directly
        println("  [$idx/$(length(csv_files))] Reading: $(basename(csv_file))")
        df = CSV.read(csv_file, DataFrame)
        push!(dfs, df)
    end

    # Concatenate all DataFrames
    merged_df = vcat(dfs...)

    println("Merged $(length(csv_files)) files into DataFrame with $(nrow(merged_df)) rows")

    # Sort if requested
    if !isnothing(sort_column)
        if sort_column in names(merged_df)
            sort!(merged_df, sort_column)
            println("Sorted by column: $sort_column")
        else
            @warn "Sort column '$sort_column' not found in DataFrame. Available columns: $(names(merged_df))"
        end
    end

    return merged_df
end


# =============================================================================
# RECONSTRUCTION RESULTS I/O (tracks + central paths)
# =============================================================================

"""
    read_central_path_from_hdf5(group)

Read a central path DataFrame from an HDF5 group.
"""
function read_central_path_from_hdf5(group)
    cp_data = read(group["central_path_data"])
    cp_columns = read(group["central_path_columns"])

    if isempty(cp_data) || isempty(cp_columns)
        return DataFrame()
    end

    df = DataFrame(cp_data, Symbol.(cp_columns))
    return df
end

"""
    read_reco_result_from_hdf5(group)

Read a single reconstruction result from an HDF5 group.

Returns NamedTuple with:
- `track`: Tracks object
- `path`: DataFrame with raw path (x, y, z, s)
- `event_id`: Event ID
- `track_length`: Total track length in mm
- `confidence`: Confidence score
- `mc_path`: DataFrame with MC path (x, y, z, energy, s)
- `reco_s`: Vector of arc-lengths for each reco voxel
- `kde_s`: Vector of KDE evaluation points
- `reco_kde_f`: Vector of RECO energy density
- `mc_kde_f`: Vector of MC energy density
- `kde_bandwidth`: RECO KDE bandwidth used
- `mc_kde_bandwidth`: MC KDE bandwidth used
"""
function read_reco_result_from_hdf5(group; diffusion::DiffusionParams=DiffusionParams())
    # Read event metadata
    event_id = read(group["event_id"])
    track_length = read(group["track_length"])
    confidence = read(group["confidence"])

    # Read KDE bandwidths (with backward compatibility)
    kde_bandwidth = haskey(group, "kde_bandwidth") ? read(group["kde_bandwidth"]) : 5.0
    mc_kde_bandwidth = haskey(group, "mc_kde_bandwidth") ? read(group["mc_kde_bandwidth"]) : kde_bandwidth

    # Read voxels
    voxels_data = read(group["voxels"])
    voxel_columns = read(group["voxel_columns"])
    voxels_df = DataFrame(voxels_data, Symbol.(voxel_columns))

    # Read graph
    n_vertices = read(group["n_vertices"])
    edge_matrix = read(group["graph_edges"])

    graph = SimpleGraph(n_vertices)
    if size(edge_matrix, 1) > 0
        for i in 1:size(edge_matrix, 1)
            add_edge!(graph, edge_matrix[i, 1], edge_matrix[i, 2])
        end
    end

    # Read components
    comp_matrix = read(group["components"])
    components = Vector{Int}[]
    if size(comp_matrix, 1) > 0
        for i in 1:size(comp_matrix, 1)
            comp = filter(x -> x > 0, comp_matrix[i, :])
            push!(components, comp)
        end
    end

    # Create track with diffusion params
    track = Tracks(voxels_df, graph, components, diffusion)

    # Read path (new format: path_data, or old format: central_path_data)
    path = DataFrame()
    if haskey(group, "path_data")
        path_data = read(group["path_data"])
        path_columns = read(group["path_columns"])
        if !isempty(path_data) && !isempty(path_columns)
            path = DataFrame(path_data, Symbol.(path_columns))
        end
    elseif haskey(group, "central_path_data")
        # Backward compatibility with old files
        path = read_central_path_from_hdf5(group)
    end

    # Read MC path (new format)
    mc_path = DataFrame()
    if haskey(group, "mc_path_data")
        mc_path_data = read(group["mc_path_data"])
        mc_path_columns = read(group["mc_path_columns"])
        if !isempty(mc_path_data) && !isempty(mc_path_columns)
            mc_path = DataFrame(mc_path_data, Symbol.(mc_path_columns))
        end
    end

    # Read projected voxel arc-lengths
    reco_s = haskey(group, "reco_s") ? read(group["reco_s"]) : Float64[]

    # Read KDE results
    kde_s = haskey(group, "kde_s") ? read(group["kde_s"]) : Float64[]
    reco_kde_f = haskey(group, "reco_kde_f") ? read(group["reco_kde_f"]) : Float64[]
    mc_kde_s = haskey(group, "mc_kde_s") ? read(group["mc_kde_s"]) : kde_s  # backward compat: use kde_s
    mc_kde_f = haskey(group, "mc_kde_f") ? read(group["mc_kde_f"]) : Float64[]

    # Read extreme distances (with backward compatibility)
    d1 = haskey(group, "d1") ? read(group["d1"]) : NaN
    d2 = haskey(group, "d2") ? read(group["d2"]) : NaN

    return (track=track, path=path,
            event_id=event_id, track_length=track_length, confidence=confidence,
            mc_path=mc_path, reco_s=reco_s,
            kde_s=kde_s, reco_kde_f=reco_kde_f,
            mc_kde_s=mc_kde_s, mc_kde_f=mc_kde_f,
            kde_bandwidth=kde_bandwidth, mc_kde_bandwidth=mc_kde_bandwidth,
            d1=d1, d2=d2)
end

"""
    read_reco_results_from_hdf5(filepath)

Read all reconstruction results from an HDF5 file produced by track_reco_mt.jl.

# Returns
- `results`: Vector of NamedTuples with fields:
  - `track`: Tracks object
  - `path`: DataFrame (raw path with x, y, z, s)
  - `event_id`, `track_length`, `confidence`
  - `mc_path`: DataFrame (MC path with x, y, z, energy, s, primary_electron)
  - `reco_s`: arc-lengths of reco voxels
  - `kde_s`, `reco_kde_f`: RECO KDE evaluation grid and energy density
  - `mc_kde_s`, `mc_kde_f`: MC KDE evaluation grid and energy density (separate from RECO)
  - `kde_bandwidth`, `mc_kde_bandwidth`: bandwidths used
  - `d1`, `d2`: distances from reco extremes to matched MC extremes
- `metadata`: Dict with file metadata (sigma_t, sigma_l, voxel_size, kde_bandwidth, etc.)

# Example
```julia
results, metadata = read_reco_results_from_hdf5("test_reco_th_1.h5")
for r in results
    println("Event \$(r.event_id): track_length=\$(r.track_length) mm")
    # Access pre-computed KDE
    plot(r.kde_s, r.reco_kde_f, label="RECO")
    plot!(r.kde_s, r.mc_kde_f, label="MC")
end
```
"""
function read_reco_results_from_hdf5(filepath::String)
    results = []
    metadata = Dict{String, Any}()

    h5open(filepath, "r") do fid
        # Read metadata from attributes
        for key in keys(attrs(fid))
            metadata[key] = HDF5.read_attribute(fid, key)
        end

        # Reconstruct DiffusionParams from metadata
        dfpars = DiffusionParams(
            get(metadata, "ldrift_cm", 100.0),
            get(metadata, "sigma_t_mm", 3.0),
            get(metadata, "sigma_l_mm", 0.0),
            get(metadata, "voxel_size_mm", 3.0),
            get(metadata, "max_distance_mm", 4.5),
            get(metadata, "energy_threshold_kev", 0.25),
            get(metadata, "nbins", 100),
            get(metadata, "nsigma", 3.0)
        )

        # Find all track groups
        if !haskey(fid, "batch_1")
            return results, metadata
        end

        batch_group = fid["batch_1"]
        track_names = sort(keys(batch_group), by=x -> parse(Int, split(x, "_")[2]))

        for track_name in track_names
            track_group = batch_group[track_name]
            result = read_reco_result_from_hdf5(track_group; diffusion=dfpars)
            push!(results, result)
        end
    end

    return results, metadata
end

"""
    chain_reco_results(filepaths::Vector{String})

Read and chain multiple reconstruction HDF5 files.

# Returns
- `all_results`: Vector of all reconstruction results
- `all_metadata`: Vector of metadata dicts from each file

# Example
```julia
files = ["test_reco_th_1.h5", "test_reco_th_2.h5"]
results, metadata = chain_reco_results(files)
```
"""
function chain_reco_results(filepaths::Vector{String})
    all_results = []
    all_metadata = Dict{String, Any}[]

    for filepath in filepaths
        results, metadata = read_reco_results_from_hdf5(filepath)
        append!(all_results, results)
        push!(all_metadata, metadata)
    end

    return all_results, all_metadata
end

"""
    chain_reco_results(directory::String, pattern::String)

Read and chain reconstruction files matching a glob pattern.

# Example
```julia
results, metadata = chain_reco_results("/path/to/data", "test_reco_th_*.h5")
```
"""
function chain_reco_results(directory::String, pattern::String)
    filepaths = sort(Glob.glob(pattern, directory))
    if isempty(filepaths)
        @warn "No files found matching pattern '$pattern' in '$directory'"
        return [], Dict{String, Any}[]
    end
    return chain_reco_results(filepaths)
end