#!/usr/bin/env julia

"""
Simple ITACA track reconstruction script.

Features:
- Reads HDF5 files and loops over events
- Can read a directory of part files (*_part_*.h5) and chain them
- Optional 3D Gaussian diffusion before voxelization
- Two modes: interactive (with plots) or batch (statistics only)
- Builds tracks and computes statistics
- In interactive mode: plots raw hits, diffused hits (if enabled), voxelized hits, and tracks
- In batch mode: collects statistics and shows summary plots at the end
- Saves results (plots as PNG, statistics as JSON) to a results directory

Usage:
    julia itaca_simple_reco.jl <cmdir> <input_file> --voxel-size=X [options]

Options:
    --voxel-size=X      Voxel size in mm (required)
    --sigmat=X          Transverse diffusion sigma in mm (default: 0, no diffusion)
    --sigmal=X          Longitudinal diffusion sigma in mm (default: 0.1, only used if sigmat > 0)
    --method=X          Track building method: "kNN" (default) or "KDT" (radius graph)
    --k=N               Number of neighbors for kNN method (default: 10)
    --max-distance=X    Max distance in mm (default: Inf for kNN, 2.5*voxel_size for KDT)
    --ievt=N            First event to process (default: 1)
    --levt=N            Last event to process (default: -1, all)
    --energy-threshold=X Energy threshold in keV (default: 1.0)
    --grid-size=X       Grid size for occupancy histograms in mm (default: 160)
    --batch             Run in batch mode (no interactive plots, show statistics at end)
    --readdir           Read directory of part files (*_part_*.h5) and chain them
    --save-results=X    Save results to <filename>_results/ directory (default: true)
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using HDF5
using DataFrames
using CSV
using Statistics
using StatsBase
using Plots
using JSON
using Dates

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

#=============================================================================
# Helper Functions
=============================================================================#

"""
    load_hits_from_directory(dirpath::String)

Load and chain all HDF5 part files (*_part_*.h5) from a directory.
Returns a single DataFrame with all hits concatenated.
"""
function load_hits_from_directory(dirpath::String)
    if !isdir(dirpath)
        error("Directory not found: $dirpath")
    end

    # Find all part files matching pattern *_part_*.h5 or *.next_part_*.h5
    all_files = readdir(dirpath)
    part_files = filter(f -> occursin(r"[._]part_\d+\.h5$", f), all_files)

    if isempty(part_files)
        error("No part files (*_part_*.h5 or *.next_part_*.h5) found in directory: $dirpath")
    end

    # Sort by part number
    sort!(part_files, by = f -> begin
        m = match(r"[._]part_(\d+)\.h5$", f)
        m !== nothing ? parse(Int, m.captures[1]) : 0
    end)

    println("Found $(length(part_files)) part files:")
    for f in part_files
        println("  - $f")
    end

    # Load and concatenate all files
    all_hits = DataFrame[]
    for f in part_files
        filepath = joinpath(dirpath, f)
        dfs = Petit.get_dataset_dfs(filepath)
        push!(all_hits, dfs["hits"])
    end

    # Concatenate all DataFrames
    hitsdf = vcat(all_hits...)
    println("Total hits loaded: $(nrow(hitsdf))")
    println("Total events: $(length(unique(hitsdf.event_id)))")

    return hitsdf
end

"""
    wait_for_enter(message::String="Press Enter to continue...")

Wait for user to press Enter before continuing.
"""
function wait_for_enter(message::String="Press Enter to continue...")
    print(message)
    readline()
    println()
end

"""
    wait_for_input(message::String="Press Enter to continue, 's' to save plot...")

Wait for user input and return the choice.
Returns:
- :continue if Enter was pressed
- :save if 's' or 'S' was pressed
- :quit if 'q' or 'Q' was pressed
"""
function wait_for_input(message::String="Press Enter to continue, 's' to save plot...")
    print(message)
    response = strip(readline())
    println()
    if lowercase(response) == "s"
        return :save
    elseif lowercase(response) == "q"
        return :quit
    else
        return :continue
    end
end

"""
    save_interactive_plot(p, results_dir::String, plot_name::String)

Save a plot to the interactive results directory.
"""
function save_interactive_plot(p, results_dir::String, plot_name::String)
    if !isdir(results_dir)
        mkpath(results_dir)
        println("Created interactive results directory: $results_dir")
    end
    filepath = joinpath(results_dir, "$(plot_name).png")
    savefig(p, filepath)
    println("Saved plot to: $filepath")
end

"""
    compute_track_z_extension(track::Petit.Tracks)

Compute the longitudinal (z) extension of a track.

Returns a NamedTuple with:
- z_min: Minimum z coordinate (mm)
- z_max: Maximum z coordinate (mm)
- z_length: Longitudinal extension z_max - z_min (mm)
"""
function compute_track_z_extension(track::Petit.Tracks)
    voxels = track.voxels

    if nrow(voxels) == 0
        return (z_min=NaN, z_max=NaN, z_length=0.0)
    end

    z_min = minimum(voxels.z)
    z_max = maximum(voxels.z)
    z_length = z_max - z_min

    return (z_min=z_min, z_max=z_max, z_length=z_length)
end

"""
    compute_track_bbox(track::Petit.Tracks)

Compute the bounding box of a track in the XY plane.

Returns a NamedTuple with:
- x_extent: Width in X direction (max(x) - min(x)) in mm
- y_extent: Width in Y direction (max(y) - min(y)) in mm
- x_min, x_max: X coordinate range
- y_min, y_max: Y coordinate range
"""
function compute_track_bbox(track::Petit.Tracks)
    voxels = track.voxels

    if nrow(voxels) == 0
        return (x_extent=0.0, y_extent=0.0, x_min=NaN, x_max=NaN, y_min=NaN, y_max=NaN)
    end

    x_min = minimum(voxels.x)
    x_max = maximum(voxels.x)
    y_min = minimum(voxels.y)
    y_max = maximum(voxels.y)

    x_extent = x_max - x_min
    y_extent = y_max - y_min

    return (x_extent=x_extent, y_extent=y_extent, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
end

# Occupancy analysis constants
const TIME_SLOT_MM = 2.0      # 2 mm time slots (1 ms = 1 mm drift)
const W_VALUE_EV = 22.0       # W-value for xenon: ~22 eV per ion pair
# Note: PIXEL_SIZE_MM is now passed as a parameter (= max_distance from track building)

# Ion threshold for pixel statistics (pixels with signal above this are counted)
const ION_THRESHOLD = 10.0

"""
    compute_occupancy_analysis(track::Petit.Tracks; pixel_size_mm::Float64, nbins::Int=50)

Compute detector pixel occupancy per time slot for a single track.

The analysis:
1. Centers an nbins x nbins pixel grid on the track barycenter
2. Slices the track into time slots (2 mm bins in z)
3. For each time slot, fills a 2D histogram counting IONS (not voxels)
4. Ion count = energy(MeV) * 1e6 / W_value (W = 22 eV for xenon)
5. Computes summary statistics across all time slots

Arguments:
- track: The track to analyze
- pixel_size_mm: Pixel size in mm (should match max_distance from track building)
- nbins: Number of bins per side (default: 50, giving 50x50 grid)

Returns a NamedTuple with:
- n_time_slots: Number of time slots
- z_length: Track z-extension (mm)
- barycenter: (x, y, z) barycenter coordinates
- avg_counts: nbins x nbins matrix of average ion counts per bin
- max_counts: nbins x nbins matrix of maximum ion counts per bin
- min_counts: nbins x nbins matrix of minimum ion counts per bin (NaN for bins never hit)
- total_counts: nbins x nbins matrix of total ion counts (sum across all time slots)
- occupancy_per_slot: Vector of occupancy fractions per time slot
- avg_occupancy: Average occupancy across all time slots
- nbins: Number of bins per side (for reference)
- pixel_size_mm: Pixel size used (for reference)
- avg_n_pixels_above_threshold: Average number of pixels with signal > ION_THRESHOLD per time slot
- avg_max_signal: Average of the maximum signal per time slot
- avg_min_signal: Average of the minimum signal (above ION_THRESHOLD) per time slot
- avg_mean_signal: Average of the mean signal (for pixels above ION_THRESHOLD) per time slot
"""
function compute_occupancy_analysis(track::Petit.Tracks; pixel_size_mm::Float64, nbins::Int=50)
    voxels = track.voxels
    pixel_grid_size = nbins * pixel_size_mm  # Total grid size in mm

    if nrow(voxels) == 0
        return (n_time_slots=0, z_length=0.0, barycenter=(x=NaN, y=NaN, z=NaN),
                avg_counts=zeros(nbins, nbins),
                max_counts=zeros(nbins, nbins),
                min_counts=fill(NaN, nbins, nbins),
                total_counts=zeros(nbins, nbins),
                occupancy_per_slot=Float64[], avg_occupancy=0.0, nbins=nbins,
                pixel_size_mm=pixel_size_mm,
                avg_n_pixels_above_threshold=0.0, avg_max_signal=0.0,
                avg_min_signal=0.0, avg_mean_signal=0.0)
    end

    # Compute barycenter and z-extension
    bary = compute_barycenter(track)
    z_ext = compute_track_z_extension(track)

    # Number of time slots (ceil to include all voxels)
    n_time_slots = max(1, ceil(Int, z_ext.z_length / TIME_SLOT_MM))

    # Initialize accumulator matrices for statistics
    sum_counts = zeros(nbins, nbins)
    max_counts = fill(-Inf, nbins, nbins)
    min_counts = fill(Inf, nbins, nbins)
    occupancy_per_slot = Float64[]

    # Per-slot statistics for new metrics
    n_pixels_above_threshold_per_slot = Float64[]
    max_signal_per_slot = Float64[]
    min_signal_per_slot = Float64[]
    mean_signal_per_slot = Float64[]

    # Define pixel grid edges centered on barycenter
    half_grid = pixel_grid_size / 2.0
    x_edges = range(bary.x - half_grid, bary.x + half_grid, length=nbins + 1)
    y_edges = range(bary.y - half_grid, bary.y + half_grid, length=nbins + 1)

    # Process each time slot
    for slot_idx in 1:n_time_slots
        z_slot_min = z_ext.z_min + (slot_idx - 1) * TIME_SLOT_MM
        z_slot_max = z_slot_min + TIME_SLOT_MM

        # Filter voxels in this time slot
        mask = (voxels.z .>= z_slot_min) .& (voxels.z .< z_slot_max)
        slot_voxels = voxels[mask, :]

        # Create 2D histogram for this slot (counting ions, not voxels)
        slot_counts = zeros(Float64, nbins, nbins)

        for i in 1:nrow(slot_voxels)
            x = slot_voxels.x[i]
            y = slot_voxels.y[i]
            energy_MeV = slot_voxels.energy[i]

            # Convert energy to number of ions: E(MeV) * 1e6 / W_value(eV)
            n_ions = energy_MeV * 1e6 / W_VALUE_EV

            # Find bin indices
            xi = searchsortedlast(collect(x_edges), x)
            yi = searchsortedlast(collect(y_edges), y)

            # Check bounds (1 to nbins)
            if 1 <= xi <= nbins && 1 <= yi <= nbins
                slot_counts[xi, yi] += n_ions
            end
        end

        # Update statistics
        sum_counts .+= slot_counts
        max_counts .= max.(max_counts, slot_counts)

        # Only update min for bins that have signal (count > 0) in this slot
        for i in 1:nbins, j in 1:nbins
            if slot_counts[i, j] > 0
                min_counts[i, j] = min(min_counts[i, j], slot_counts[i, j])
            end
        end

        # Compute occupancy for this slot (fraction of bins with signal)
        n_occupied = count(slot_counts .> 0)
        occupancy = n_occupied / (nbins * nbins)
        push!(occupancy_per_slot, occupancy)

        # Compute new pixel statistics for this slot (pixels above ION_THRESHOLD)
        pixels_above = slot_counts[slot_counts .> ION_THRESHOLD]
        n_above = length(pixels_above)
        push!(n_pixels_above_threshold_per_slot, Float64(n_above))

        if n_above > 0
            push!(max_signal_per_slot, maximum(pixels_above))
            push!(min_signal_per_slot, minimum(pixels_above))
            push!(mean_signal_per_slot, mean(pixels_above))
        else
            # No pixels above threshold in this slot - use 0 for max/mean, NaN for min
            push!(max_signal_per_slot, 0.0)
            push!(min_signal_per_slot, NaN)
            push!(mean_signal_per_slot, 0.0)
        end
    end

    # Compute average counts
    avg_counts = sum_counts ./ n_time_slots

    # Set bins that were never hit to NaN (for min_counts) or 0 (for max_counts)
    min_counts[min_counts .== Inf] .= NaN
    max_counts[max_counts .== -Inf] .= 0.0

    # Average occupancy
    avg_occupancy = isempty(occupancy_per_slot) ? 0.0 : mean(occupancy_per_slot)

    # Compute averages of new metrics
    avg_n_pixels_above_threshold = isempty(n_pixels_above_threshold_per_slot) ? 0.0 : mean(n_pixels_above_threshold_per_slot)
    avg_max_signal = isempty(max_signal_per_slot) ? 0.0 : mean(max_signal_per_slot)
    # For min signal, only average over slots that had pixels above threshold
    valid_min_signals = filter(!isnan, min_signal_per_slot)
    avg_min_signal = isempty(valid_min_signals) ? 0.0 : mean(valid_min_signals)
    avg_mean_signal = isempty(mean_signal_per_slot) ? 0.0 : mean(mean_signal_per_slot)

    return (n_time_slots=n_time_slots, z_length=z_ext.z_length,
            barycenter=(x=bary.x, y=bary.y, z=bary.z),
            avg_counts=avg_counts, max_counts=max_counts, min_counts=min_counts,
            total_counts=sum_counts,
            occupancy_per_slot=occupancy_per_slot, avg_occupancy=avg_occupancy,
            nbins=nbins, pixel_size_mm=pixel_size_mm,
            avg_n_pixels_above_threshold=avg_n_pixels_above_threshold,
            avg_max_signal=avg_max_signal, avg_min_signal=avg_min_signal,
            avg_mean_signal=avg_mean_signal)
end

"""
    plot_occupancy_heatmaps(occ; event_id::Int=0)

Plot four 2D heatmaps showing pixel occupancy analysis:
1. Average counts per bin (across all time slots)
2. Maximum counts per bin
3. Minimum counts per bin (only bins with at least 1 count shown)
4. Total counts (projected histogram - sum across all time slots)

Arguments:
- occ: NamedTuple from compute_occupancy_analysis()
- event_id: Event ID for title (optional)

Returns a combined plot with 4 heatmaps in a 2x2 layout.
"""
function plot_occupancy_heatmaps(occ; event_id::Int=0)
    # Get nbins and pixel_size from occupancy result
    nbins = occ.nbins
    pixel_size_mm = occ.pixel_size_mm
    pixel_grid_size = nbins * pixel_size_mm

    # Create axis labels (relative to barycenter)
    half_grid = pixel_grid_size / 2.0

    # Compute tick positions and labels dynamically based on nbins
    tick_positions = [1, round(Int, nbins/4), round(Int, nbins/2), round(Int, 3*nbins/4), nbins]
    tick_labels = [string(round(Int, -half_grid)),
                   string(round(Int, -half_grid/2)),
                   "0",
                   string(round(Int, half_grid/2)),
                   string(round(Int, half_grid))]

    # Common plot settings
    common_opts = (
        xlabel="X (mm)", ylabel="Y (mm)",
        xticks=(tick_positions, tick_labels),
        yticks=(tick_positions, tick_labels),
        aspect_ratio=:equal,
        colorbar=true,
        titlefontsize=9
    )

    # Heatmap 1: Average ion counts
    p1 = heatmap(occ.avg_counts', c=:viridis, title="Avg Ions/Bin"; common_opts...)

    # Heatmap 2: Maximum ion counts
    p2 = heatmap(occ.max_counts', c=:hot, title="Max Ions/Bin"; common_opts...)

    # Heatmap 3: Minimum ion counts (already has NaN for bins never hit)
    p3 = heatmap(occ.min_counts', c=:blues, title="Min Ions/Bin (>0)"; common_opts...)

    # Heatmap 4: Total ion counts (projected)
    p4 = heatmap(occ.total_counts', c=:inferno, title="Total Ions (Projected)"; common_opts...)

    # Combine into 2x2 layout
    pixel_str = isinteger(pixel_size_mm) ? "$(Int(pixel_size_mm))" : "$pixel_size_mm"
    title_str = event_id > 0 ? "Ion Occupancy - Event $event_id ($(occ.n_time_slots) slots, $(nbins)x$(nbins), $(pixel_str)mm pixels)" :
                               "Ion Occupancy ($(occ.n_time_slots) slots, $(nbins)x$(nbins), $(pixel_str)mm pixels)"
    plt = plot(p1, p2, p3, p4, layout=(2, 2), size=(800, 750),
               plot_title=title_str, plot_titlefontsize=10)

    return plt
end

"""
    compute_barycenter(track::Petit.Tracks)

Compute the energy-weighted barycenter of a track.

Returns a NamedTuple with:
- x, y, z: Barycenter coordinates (mm)
- total_energy: Total energy in the track (keV)
"""
function compute_barycenter(track::Petit.Tracks)
    voxels = track.voxels

    if nrow(voxels) == 0
        return (x=NaN, y=NaN, z=NaN, total_energy=0.0)
    end

    # Energy in keV (voxels.energy is in MeV)
    energies_keV = 1e+3 .* voxels.energy
    total_energy = sum(energies_keV)

    if total_energy == 0.0
        # If no energy, use simple centroid
        return (x=mean(voxels.x), y=mean(voxels.y), z=mean(voxels.z),
                total_energy=0.0)
    end

    # Energy-weighted barycenter
    x_bary = sum(voxels.x .* energies_keV) / total_energy
    y_bary = sum(voxels.y .* energies_keV) / total_energy
    z_bary = sum(voxels.z .* energies_keV) / total_energy

    return (x=x_bary, y=y_bary, z=z_bary, total_energy=total_energy)
end

"""
    plot_all_tracks(tracks::Vector{Petit.Tracks}, event_id::Int)

Plot all tracks in a single figure:
- Track 1 (main): viridis colormap based on energy
- Track 2: solid blue
- Track 3: solid green
- Track 4+: solid black

Returns a plot with XY, XZ, YZ projections and 3D scatter.
"""
function plot_all_tracks(tracks::Vector{Petit.Tracks}, event_id::Int)
    n_tracks = length(tracks)
    if n_tracks == 0
        return plot(title="No tracks")
    end

    # Collect all voxels for axis limits
    all_x, all_y, all_z = Float64[], Float64[], Float64[]
    for t in tracks
        append!(all_x, t.voxels.x)
        append!(all_y, t.voxels.y)
        append!(all_z, t.voxels.z)
    end

    # Compute padded limits
    xmid, xrange = mean((minimum(all_x), maximum(all_x))), maximum(all_x) - minimum(all_x)
    ymid, yrange = mean((minimum(all_y), maximum(all_y))), maximum(all_y) - minimum(all_y)
    zmid, zrange = mean((minimum(all_z), maximum(all_z))), maximum(all_z) - minimum(all_z)

    xlim = (xmid - 0.65 * max(xrange, 1.0), xmid + 0.65 * max(xrange, 1.0))
    ylim = (ymid - 0.65 * max(yrange, 1.0), ymid + 0.65 * max(yrange, 1.0))
    zlim = (zmid - 0.65 * max(zrange, 1.0), zmid + 0.65 * max(zrange, 1.0))

    # Track colors: track 2 = blue, track 3 = green, track 4+ = black
    satellite_colors = [:blue, :green, :black]

    # Main track data
    main = tracks[1].voxels
    x1, y1, z1 = Float64.(main.x), Float64.(main.y), Float64.(main.z)
    e1 = Float64.(main.energy)

    cmap = cgrad(:viridis)
    ms_main = 4.0
    ms_sat = 6.0  # Slightly larger for satellites to stand out

    # === XY Projection ===
    p1 = scatter(x1, y1, marker_z=e1, ms=ms_main, markerstrokewidth=0,
                 xlabel="X (mm)", ylabel="Y (mm)", title="XY Projection",
                 xlims=xlim, ylims=ylim, color=cmap, colorbar=false, label="Track 1")

    # Overlay satellites
    for i in 2:n_tracks
        sat = tracks[i].voxels
        col = i <= 3 ? satellite_colors[i-1] : :black
        label = i <= 3 ? "Track $i" : (i == 4 ? "Track 4+" : "")
        scatter!(p1, Float64.(sat.x), Float64.(sat.y), ms=ms_sat, color=col,
                 markerstrokewidth=1, markerstrokecolor=:white, label=label)
    end

    # === XZ Projection ===
    p2 = scatter(x1, z1, marker_z=e1, ms=ms_main, markerstrokewidth=0,
                 xlabel="X (mm)", ylabel="Z (mm)", title="XZ Projection",
                 xlims=xlim, ylims=zlim, color=cmap, colorbar=false, label="")

    for i in 2:n_tracks
        sat = tracks[i].voxels
        col = i <= 3 ? satellite_colors[i-1] : :black
        scatter!(p2, Float64.(sat.x), Float64.(sat.z), ms=ms_sat, color=col,
                 markerstrokewidth=1, markerstrokecolor=:white, label="")
    end

    # === YZ Projection ===
    p3 = scatter(y1, z1, marker_z=e1, ms=ms_main, markerstrokewidth=0,
                 xlabel="Y (mm)", ylabel="Z (mm)", title="YZ Projection",
                 xlims=ylim, ylims=zlim, color=cmap, colorbar=false, label="")

    for i in 2:n_tracks
        sat = tracks[i].voxels
        col = i <= 3 ? satellite_colors[i-1] : :black
        scatter!(p3, Float64.(sat.y), Float64.(sat.z), ms=ms_sat, color=col,
                 markerstrokewidth=1, markerstrokecolor=:white, label="")
    end

    # === 3D Scatter ===
    p4 = scatter(x1, y1, z1, marker_z=e1, ms=ms_main, markerstrokewidth=0,
                 xlabel="X", ylabel="Y", zlabel="Z", title="3D View",
                 color=cmap, colorbar_title="Energy (MeV)", label="", camera=(30, 30))

    for i in 2:n_tracks
        sat = tracks[i].voxels
        col = i <= 3 ? satellite_colors[i-1] : :black
        scatter!(p4, Float64.(sat.x), Float64.(sat.y), Float64.(sat.z),
                 ms=ms_sat, color=col, markerstrokewidth=1, markerstrokecolor=:white, label="")
    end

    # Combine into layout
    l = @layout [a b; c d]
    plt = plot(p1, p2, p3, p4, layout=l, size=(1000, 800),
               plot_title="Event $event_id - $n_tracks tracks (T1=viridis, T2=blue, T3=green, T4+=black)")

    return plt
end

"""
    closest_voxel_distance(track1::Petit.Tracks, track2::Petit.Tracks)

Compute the minimum distance between any voxel in track1 and any voxel in track2.
Returns the minimum distance and the indices of the closest voxel pair.
"""
function closest_voxel_distance(track1::Petit.Tracks, track2::Petit.Tracks)
    vox1 = track1.voxels
    vox2 = track2.voxels

    min_dist = Inf
    closest_i, closest_j = 0, 0

    for i in 1:nrow(vox1)
        for j in 1:nrow(vox2)
            dist = Petit.euclidean_distance(
                vox1.x[i], vox1.y[i], vox1.z[i],
                vox2.x[j], vox2.y[j], vox2.z[j]
            )
            if dist < min_dist
                min_dist = dist
                closest_i, closest_j = i, j
            end
        end
    end

    return (distance=min_dist, voxel_idx_track1=closest_i, voxel_idx_track2=closest_j)
end

"""
    compute_track_distances(tracks::Vector{Petit.Tracks})

Compute distances from the main track (highest energy) to each satellite track.
Uses closest voxel pair distance (minimum distance between any two voxels).

Returns a NamedTuple with:
- main_track_idx: Index of the main track (1-based)
- satellite_distances: Vector of closest distances from main track to each satellite (mm)
- satellite_indices: Vector of satellite track indices
"""
function compute_track_distances(tracks::Vector{Petit.Tracks})
    n_tracks = length(tracks)

    if n_tracks == 0
        return (main_track_idx=0, satellite_distances=Float64[], satellite_indices=Int[])
    end

    if n_tracks == 1
        return (main_track_idx=1, satellite_distances=Float64[], satellite_indices=Int[])
    end

    # Main track is index 1 (tracks are already sorted by energy)
    main_idx = 1
    main_track = tracks[main_idx]

    # Compute closest voxel distances to satellites
    satellite_distances = Float64[]
    satellite_indices = Int[]

    for i in 2:n_tracks
        result = closest_voxel_distance(main_track, tracks[i])
        push!(satellite_distances, result.distance)
        push!(satellite_indices, i)
    end

    return (main_track_idx=main_idx,
            satellite_distances=satellite_distances,
            satellite_indices=satellite_indices)
end

"""
    print_track_summary(tracks::Vector{Petit.Tracks}, event_id::Int)

Print detailed summary of all tracks for an event.
Waits for user to press Enter after displaying.
"""
function print_track_summary(tracks::Vector{Petit.Tracks}, event_id::Int)
    n_tracks = length(tracks)

    println()
    println("═" ^ 70)
    println("          TRACK SUMMARY - Event $event_id")
    println("═" ^ 70)
    println()
    println("Number of tracks: $n_tracks")
    println()

    if n_tracks == 0
        println("No tracks found for this event.")
        println("═" ^ 70)
        wait_for_enter()
        return
    end

    # Get track statistics
    stats = Petit.track_stats(tracks)

    # Compute barycenters
    barycenters = [compute_barycenter(t) for t in tracks]

    # Print info for each track
    for i in 1:n_tracks
        track_type = i == 1 ? "MAIN TRACK" : "Satellite $i"
        println("─" ^ 70)
        println("  Track $i ($track_type)")
        println("─" ^ 70)
        println("    Number of voxels:  $(stats.n_voxels[i])")
        println("    Total energy:      $(round(stats.energies_keV[i], digits=2)) keV")

        bary = barycenters[i]
        println("    Barycenter (x,y,z): ($(round(bary.x, digits=2)), $(round(bary.y, digits=2)), $(round(bary.z, digits=2))) mm")
        println()
    end

    # Compute and print distances between main track and satellites
    if n_tracks > 1
        distances = compute_track_distances(tracks)

        println("─" ^ 70)
        println("  DISTANCES FROM MAIN TRACK TO SATELLITES")
        println("─" ^ 70)

        for (i, (sat_idx, dist)) in enumerate(zip(distances.satellite_indices, distances.satellite_distances))
            println("    Main → Satellite $sat_idx: $(round(dist, digits=2)) mm")
        end
        println()
    end

    println("═" ^ 70)
    wait_for_enter()
end

#=============================================================================
# Statistics Structure
=============================================================================#

"""
    EventStats

Structure to hold statistics for a single event.
"""
struct EventStats
    event_id::Int
    n_hits::Int
    n_voxels::Int
    n_tracks::Int
    total_energy_keV::Float64
    main_track_energy_keV::Float64
    main_track_voxels::Int
    track2_energy_keV::Float64          # Energy of track 2 (0 if doesn't exist)
    track2_voxels::Int                  # Number of voxels in track 2 (0 if doesn't exist)
    track3plus_energy_keV::Float64      # Combined energy of tracks 3, 4, 5, ... (0 if none)
    satellite_distances::Vector{Float64}
    voxel_min_distances::Vector{Float64}  # Closest neighbor distance for each voxel
    # Single-track analysis (only computed for 1-track events)
    z_length::Float64                   # Longitudinal extension of main track (mm)
    x_extent::Float64                   # Track bounding box X extent (mm)
    y_extent::Float64                   # Track bounding box Y extent (mm)
    avg_occupancy::Float64              # Average pixel occupancy per time slot
    n_time_slots::Int                   # Number of 2mm time slots
    occupancy_per_slot::Vector{Float64} # Occupancy per time slot
    # Pixel signal statistics (single-track, averaged over time slots)
    avg_n_pixels_above_threshold::Float64  # Avg number of pixels with signal > 10 ions
    avg_max_signal::Float64                # Avg of max signal per time slot
    avg_min_signal::Float64                # Avg of min signal (above threshold) per time slot
    avg_mean_signal::Float64               # Avg of mean signal (above threshold) per time slot
end

"""
    plot_statistics(stats::Vector{EventStats}; grid_size::Float64=160.0, pixel_size_mm::Float64=2.0)

Plot summary statistics from batch processing.
Returns three plot frames:
- Frame 1 (2x2): Track multiplicity, Total energy, Pixel occupancy %, Track 2 energy
- Frame 2 (1x2): Z-length histogram, Track bounding box 2D histogram
- Frame 3 (2x2): Ion occupancy heatmaps (sample event)
"""
function plot_statistics(stats::Vector{EventStats}; grid_size::Float64=160.0, pixel_size_mm::Float64=2.0)
    n_events = length(stats)
    if n_events == 0
        println("No events to plot.")
        return nothing, nothing, nothing
    end

    # Compute nbins from grid_size / pixel_size
    nbins = round(Int, grid_size / pixel_size_mm)

    # Extract data
    n_tracks_vec = [s.n_tracks for s in stats]
    total_energies = [s.total_energy_keV for s in stats]

    # Count track multiplicities
    max_tracks = maximum(n_tracks_vec)
    track_counts = [count(==(i), n_tracks_vec) for i in 0:max_tracks]

    # =========================================================================
    # FRAME 1 (2x2): Track multiplicity, Total energy, Pixel occupancy, Track 2 energy
    # =========================================================================

    # === Plot 1: Number of tracks histogram ===
    p1 = bar(0:max_tracks, track_counts,
             xlabel="N Tracks", ylabel="Events",
             title="Track Multiplicity", titlefontsize=9,
             legend=false, color=:steelblue, bar_width=0.8)

    # Add percentage labels (offset to the right, only first 3 bars)
    max_count = maximum(track_counts)
    y_offset = 0.01 * max_count  # 1% of max height (minimal offset)
    for (i, cnt) in enumerate(track_counts)
        if i > 3
            break  # Only label first 3 bars (0, 1, 2 tracks)
        end
        pct = round(100 * cnt / n_events, digits=1)
        if cnt > 0
            annotate!(p1, i - 1 + 0.45, cnt + y_offset, text("$pct%", 7, :left))
        end
    end

    # === Plot 2: Total energy distribution ===
    p2 = histogram(total_energies, bins=50,
                   xlabel="Energy (keV)", ylabel="Events",
                   title="Total Energy", titlefontsize=9,
                   legend=false, color=:red)

    # === Plot 3: Pixel occupancy % (for single-track events only) ===
    avg_occupancies = [s.avg_occupancy for s in stats if s.n_tracks == 1 && s.avg_occupancy > 0]
    if !isempty(avg_occupancies)
        p3 = histogram(100.0 .* avg_occupancies, bins=30,
                       xlabel="Avg Occupancy (%)", ylabel="Events",
                       title="Pixel Occupancy (1-trk)", titlefontsize=9,
                       legend=false, color=:coral)
    else
        p3 = plot(title="Pixel Occupancy\n(no 1-trk data)", titlefontsize=9, legend=false)
    end

    # === Plot 4: Track 2 energy (for events with ≥2 tracks) ===
    track2_energies = [s.track2_energy_keV for s in stats if s.n_tracks >= 2]
    if !isempty(track2_energies)
        p4 = histogram(track2_energies, bins=50,
                       xlabel="Energy (keV)", ylabel="Events",
                       title="Track 2 Energy", titlefontsize=9,
                       legend=false, color=:blue)
    else
        p4 = plot(title="Track 2 Energy\n(no data)", titlefontsize=9, legend=false)
    end

    # Combine into layout (2x2) - Frame 1
    plt1 = plot(p1, p2, p3, p4, layout=(2, 2), size=(900, 700),
                plot_title="Statistics (N=$n_events)", plot_titlefontsize=10)

    # =========================================================================
    # FRAME 2 (1x2): Z-length histogram, Track bounding box 2D histogram
    # =========================================================================

    # Z-length histogram (for single-track events only)
    z_lengths = [s.z_length for s in stats if s.n_tracks == 1 && s.z_length > 0]
    if !isempty(z_lengths)
        p5 = histogram(z_lengths, bins=50,
                       xlabel="Z-length (mm)", ylabel="Events",
                       title="Track Z-Extension (1-trk)", titlefontsize=9,
                       legend=false, color=:teal)
    else
        p5 = plot(title="Track Z-Extension\n(no 1-trk data)", titlefontsize=9, legend=false)
    end

    # Track bounding box 2D histogram (for single-track events only)
    x_extents = [s.x_extent for s in stats if s.n_tracks == 1 && s.x_extent > 0]
    y_extents = [s.y_extent for s in stats if s.n_tracks == 1 && s.y_extent > 0]

    if !isempty(x_extents) && !isempty(y_extents)
        # Range: 0 to grid_size, bins = grid_size / pixel_size
        edges = range(0, grid_size, length=nbins + 1)
        p6 = histogram2d(x_extents, y_extents, bins=(edges, edges),
                         xlabel="X Extent (mm)", ylabel="Y Extent (mm)",
                         title="Track Bounding Box (1-trk)", titlefontsize=9,
                         c=:viridis, colorbar=true, aspect_ratio=:equal)
    else
        p6 = plot(title="Track Bounding Box\n(no 1-trk data)", titlefontsize=9, legend=false)
    end

    # Frame 2: Z-length and bounding box
    n_single_track = count(s -> s.n_tracks == 1, stats)
    plt2 = plot(p5, p6, layout=(1, 2), size=(1000, 450),
                plot_title="Single-Track Analysis (N=$n_single_track)", plot_titlefontsize=10)

    # Frame 3 is handled separately (occupancy heatmaps from sample event)
    # Return nothing for plt3 - it will be generated from sample_occ in batch mode
    return plt1, plt2, nothing
end

"""
    print_statistics_summary(stats::Vector{EventStats})

Print a text summary of batch statistics.
"""
function print_statistics_summary(stats::Vector{EventStats})
    n_events = length(stats)

    println("\n" * "═" ^ 70)
    println("                    BATCH PROCESSING SUMMARY")
    println("═" ^ 70)
    println("Total events processed: $n_events")
    println()

    # Track multiplicity
    n_tracks_vec = [s.n_tracks for s in stats]
    max_tracks = maximum(n_tracks_vec)

    println("Track Multiplicity:")
    for i in 0:min(max_tracks, 5)
        count = count_equal(n_tracks_vec, i)
        pct = round(100 * count / n_events, digits=1)
        label = i == 0 ? "0 tracks" : (i == 1 ? "1 track " : "$i tracks")
        println("  $label: $count events ($pct%)")
    end
    if max_tracks > 5
        count = sum(n_tracks_vec .> 5)
        pct = round(100 * count / n_events, digits=1)
        println("  >5 tracks: $count events ($pct%)")
    end

    println()
    println("Energy Statistics:")
    total_energies = [s.total_energy_keV for s in stats]
    println("  Total energy:  mean=$(round(mean(total_energies), digits=1)) keV, std=$(round(std(total_energies), digits=1)) keV")

    main_energies = [s.main_track_energy_keV for s in stats if s.n_tracks > 0]
    if !isempty(main_energies)
        println("  Main track:    mean=$(round(mean(main_energies), digits=1)) keV, std=$(round(std(main_energies), digits=1)) keV")
    end

    track2_energies = [s.track2_energy_keV for s in stats if s.n_tracks >= 2]
    if !isempty(track2_energies)
        println("  Track 2:       mean=$(round(mean(track2_energies), digits=1)) keV, std=$(round(std(track2_energies), digits=1)) keV (N=$(length(track2_energies)) events)")
    end

    track3plus_energies = [s.track3plus_energy_keV for s in stats if s.n_tracks >= 3]
    if !isempty(track3plus_energies)
        println("  Track 3+:      mean=$(round(mean(track3plus_energies), digits=1)) keV, std=$(round(std(track3plus_energies), digits=1)) keV (N=$(length(track3plus_energies)) events)")
    end

    println()
    println("Voxel Statistics:")
    n_voxels_vec = [s.n_voxels for s in stats]
    println("  Voxels/event:  mean=$(round(mean(n_voxels_vec), digits=1)), std=$(round(std(n_voxels_vec), digits=1))")

    # Voxel distance statistics
    all_voxel_dists = Float64[]
    for s in stats
        append!(all_voxel_dists, s.voxel_min_distances)
    end
    if !isempty(all_voxel_dists)
        println("  Voxel nearest neighbor dist: mean=$(round(mean(all_voxel_dists), digits=2)) mm, std=$(round(std(all_voxel_dists), digits=2)) mm")
    end

    # Single-track analysis
    single_track_stats = filter(s -> s.n_tracks == 1, stats)
    n_single = length(single_track_stats)
    if n_single > 0
        println()
        println("Single-Track Analysis (N=$n_single events):")

        z_lengths = [s.z_length for s in single_track_stats if s.z_length > 0]
        if !isempty(z_lengths)
            println("  Z-extension: mean=$(round(mean(z_lengths), digits=1)) mm, std=$(round(std(z_lengths), digits=1)) mm")
        end

        x_extents = [s.x_extent for s in single_track_stats if s.x_extent > 0]
        if !isempty(x_extents)
            println("  X-extent: mean=$(round(mean(x_extents), digits=1)) mm, std=$(round(std(x_extents), digits=1)) mm")
        end

        y_extents = [s.y_extent for s in single_track_stats if s.y_extent > 0]
        if !isempty(y_extents)
            println("  Y-extent: mean=$(round(mean(y_extents), digits=1)) mm, std=$(round(std(y_extents), digits=1)) mm")
        end

        avg_occupancies = [100 * s.avg_occupancy for s in single_track_stats if s.avg_occupancy > 0]
        if !isempty(avg_occupancies)
            println("  Avg occupancy: mean=$(round(mean(avg_occupancies), digits=1))%, std=$(round(std(avg_occupancies), digits=1))%")
        end
    end

    println("═" ^ 70)
end

# Helper function for counting
count_equal(vec, val) = count(==(val), vec)

#=============================================================================
# Results Saving Functions
=============================================================================#

"""
    create_results_dir(input_file::String, voxel_size::Float64, max_distance::Float64;
                       sigmat::Float64=0.0, sigmal::Float64=0.1,
                       energy_threshold::Float64=0.0,
                       ievt::Int=1, levt::Int=-1)

Create a results directory based on the input filename and parameters.
The directory is created under the script directory (itacaScripts/).
Directory name format:
- Without diffusion: <basename>_voxel_<X>mm_maxd_<Y>mm_ethr_<E>keV_evt_<ievt>_<levt>
- With diffusion: <basename>_sigmat_<S>mm_sigmal_<L>mm_voxel_<X>mm_maxd_<Y>mm_ethr_<E>keV_evt_<ievt>_<levt>
Returns the path to the created directory.
"""
function create_results_dir(input_file::String, voxel_size::Float64, max_distance::Float64;
                            sigmat::Float64=0.0, sigmal::Float64=0.1,
                            energy_threshold::Float64=0.0,
                            ievt::Int=1, levt::Int=-1)
    # Get the directory where this script lives
    script_dir = dirname(@__FILE__)

    # Remove extension from filename
    base_name = splitext(input_file)[1]

    # Format voxel size (remove trailing zeros)
    voxel_str = replace(string(round(voxel_size, digits=2)), r"\.0$" => "")

    # Format max distance
    if isinf(max_distance)
        maxd_str = "Inf"
    else
        maxd_str = replace(string(round(max_distance, digits=2)), r"\.0$" => "")
    end

    # Format event range
    levt_str = levt < 0 ? "end" : string(levt)
    evt_str = "evt_$(ievt)_$(levt_str)"

    # Format energy threshold
    ethr_str = replace(string(round(energy_threshold, digits=1)), r"\.0$" => "")

    # Build directory name
    if sigmat > 0.0
        # Include diffusion parameters
        sigmat_str = replace(string(round(sigmat, digits=2)), r"\.0$" => "")
        sigmal_str = replace(string(round(sigmal, digits=2)), r"\.0$" => "")
        results_dir = joinpath(script_dir, "$(base_name)_sigmat_$(sigmat_str)mm_sigmal_$(sigmal_str)mm_voxel_$(voxel_str)mm_maxd_$(maxd_str)mm_ethr_$(ethr_str)keV_$(evt_str)")
    else
        results_dir = joinpath(script_dir, "$(base_name)_voxel_$(voxel_str)mm_maxd_$(maxd_str)mm_ethr_$(ethr_str)keV_$(evt_str)")
    end

    if !isdir(results_dir)
        mkpath(results_dir)
        println("Created results directory: $results_dir")
    else
        println("Results directory exists: $results_dir")
    end

    return results_dir
end

"""
    save_statistics_json(stats::Vector{EventStats}, results_dir::String, params::Dict;
                         nbins::Int=50, pixel_size_mm::Float64=2.0, grid_size::Float64=160.0)

Save statistics to a JSON file.
"""
function save_statistics_json(stats::Vector{EventStats}, results_dir::String, params::Dict;
                              nbins::Int=50, pixel_size_mm::Float64=2.0, grid_size::Float64=160.0)
    n_events = length(stats)

    # Extract summary statistics
    n_tracks_vec = [s.n_tracks for s in stats]
    total_energies = [s.total_energy_keV for s in stats]
    main_energies = [s.main_track_energy_keV for s in stats if s.n_tracks > 0]
    track2_energies = [s.track2_energy_keV for s in stats if s.n_tracks >= 2]
    track3plus_energies = [s.track3plus_energy_keV for s in stats if s.n_tracks >= 3]
    track1_voxels = [s.main_track_voxels for s in stats if s.n_tracks > 0]
    track2_voxels = [s.track2_voxels for s in stats if s.n_tracks >= 2]
    n_voxels_vec = [s.n_voxels for s in stats]

    all_voxel_dists = Float64[]
    for s in stats
        append!(all_voxel_dists, s.voxel_min_distances)
    end

    # Single-track analysis data
    single_track_stats = filter(s -> s.n_tracks == 1, stats)
    z_lengths = [s.z_length for s in single_track_stats if s.z_length > 0]
    x_extents = [s.x_extent for s in single_track_stats if s.x_extent > 0]
    y_extents = [s.y_extent for s in single_track_stats if s.y_extent > 0]
    avg_occupancies = [s.avg_occupancy for s in single_track_stats if s.avg_occupancy > 0]

    # Track multiplicity counts
    max_tracks = maximum(n_tracks_vec)
    track_mult = Dict{String, Int}()
    for i in 0:max_tracks
        track_mult["$i"] = count_equal(n_tracks_vec, i)
    end

    # Build JSON structure
    results = Dict(
        "metadata" => Dict(
            "timestamp" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
            "n_events" => n_events,
            "parameters" => params
        ),
        "track_multiplicity" => Dict(
            "counts" => track_mult,
            "percentages" => Dict(k => round(100 * v / n_events, digits=2) for (k, v) in track_mult)
        ),
        "energy_statistics" => Dict(
            "total_energy_keV" => Dict(
                "mean" => round(mean(total_energies), digits=2),
                "std" => round(std(total_energies), digits=2),
                "min" => round(minimum(total_energies), digits=2),
                "max" => round(maximum(total_energies), digits=2)
            ),
            "main_track_energy_keV" => isempty(main_energies) ? nothing : Dict(
                "mean" => round(mean(main_energies), digits=2),
                "std" => round(std(main_energies), digits=2),
                "n_events" => length(main_energies)
            ),
            "track2_energy_keV" => isempty(track2_energies) ? nothing : Dict(
                "mean" => round(mean(track2_energies), digits=2),
                "std" => round(std(track2_energies), digits=2),
                "n_events" => length(track2_energies)
            ),
            "track3plus_energy_keV" => isempty(track3plus_energies) ? nothing : Dict(
                "mean" => round(mean(track3plus_energies), digits=2),
                "std" => round(std(track3plus_energies), digits=2),
                "n_events" => length(track3plus_energies)
            )
        ),
        "voxel_statistics" => Dict(
            "voxels_per_event" => Dict(
                "mean" => round(mean(n_voxels_vec), digits=2),
                "std" => round(std(n_voxels_vec), digits=2)
            ),
            "track1_voxels" => isempty(track1_voxels) ? nothing : Dict(
                "mean" => round(mean(track1_voxels), digits=2),
                "std" => round(std(track1_voxels), digits=2)
            ),
            "track2_voxels" => isempty(track2_voxels) ? nothing : Dict(
                "mean" => round(mean(track2_voxels), digits=2),
                "std" => round(std(track2_voxels), digits=2)
            ),
            "nearest_neighbor_dist_mm" => isempty(all_voxel_dists) ? nothing : Dict(
                "mean" => round(mean(all_voxel_dists), digits=3),
                "std" => round(std(all_voxel_dists), digits=3),
                "min" => round(minimum(all_voxel_dists), digits=3),
                "max" => round(maximum(all_voxel_dists), digits=3)
            )
        ),
        "single_track_analysis" => Dict(
            "n_single_track_events" => length(single_track_stats),
            "z_length_mm" => isempty(z_lengths) ? nothing : Dict(
                "mean" => round(mean(z_lengths), digits=2),
                "std" => round(std(z_lengths), digits=2),
                "min" => round(minimum(z_lengths), digits=2),
                "max" => round(maximum(z_lengths), digits=2)
            ),
            "x_extent_mm" => isempty(x_extents) ? nothing : Dict(
                "mean" => round(mean(x_extents), digits=2),
                "std" => round(std(x_extents), digits=2),
                "min" => round(minimum(x_extents), digits=2),
                "max" => round(maximum(x_extents), digits=2)
            ),
            "y_extent_mm" => isempty(y_extents) ? nothing : Dict(
                "mean" => round(mean(y_extents), digits=2),
                "std" => round(std(y_extents), digits=2),
                "min" => round(minimum(y_extents), digits=2),
                "max" => round(maximum(y_extents), digits=2)
            ),
            "avg_occupancy_pct" => isempty(avg_occupancies) ? nothing : Dict(
                "mean" => round(100 * mean(avg_occupancies), digits=2),
                "std" => round(100 * std(avg_occupancies), digits=2),
                "min" => round(100 * minimum(avg_occupancies), digits=2),
                "max" => round(100 * maximum(avg_occupancies), digits=2)
            ),
            "pixel_grid" => Dict(
                "nbins" => nbins,
                "pixel_size_mm" => pixel_size_mm,
                "grid_size_mm" => grid_size,
                "time_slot_mm" => TIME_SLOT_MM
            )
        ),
        "per_event_data" => [
            Dict(
                "event_id" => s.event_id,
                "n_hits" => s.n_hits,
                "n_voxels" => s.n_voxels,
                "n_tracks" => s.n_tracks,
                "total_energy_keV" => round(s.total_energy_keV, digits=2),
                "main_track_energy_keV" => round(s.main_track_energy_keV, digits=2),
                "main_track_voxels" => s.main_track_voxels,
                "track2_energy_keV" => round(s.track2_energy_keV, digits=2),
                "track2_voxels" => s.track2_voxels,
                "track3plus_energy_keV" => round(s.track3plus_energy_keV, digits=2),
                "z_length_mm" => s.n_tracks == 1 ? round(s.z_length, digits=2) : nothing,
                "x_extent_mm" => s.n_tracks == 1 ? round(s.x_extent, digits=2) : nothing,
                "y_extent_mm" => s.n_tracks == 1 ? round(s.y_extent, digits=2) : nothing,
                "avg_occupancy_pct" => s.n_tracks == 1 ? round(100 * s.avg_occupancy, digits=2) : nothing,
                "n_time_slots" => s.n_tracks == 1 ? s.n_time_slots : nothing
            ) for s in stats
        ]
    )

    # Save to file
    json_path = joinpath(results_dir, "statistics.json")
    open(json_path, "w") do f
        JSON.print(f, results, 2)  # indent=2 for readability
    end
    println("Saved statistics to: $json_path")

    return json_path
end

"""
    save_statistics_csv(stats::Vector{EventStats}, results_dir::String)

Save statistics to two CSV files for easy merging of independent batch jobs:

1. `event_stats.csv` - All events with general track statistics:
   - event_id, n_hits, n_voxels, n_tracks, total_energy_keV
   - main_track_energy_keV, main_track_voxels
   - track2_energy_keV, track2_voxels, track3plus_energy_keV
   - track1_track2_distance_mm (closest voxel distance between track 1 and 2, 0 if no track 2)

2. `single_track_stats.csv` - Only events with n_tracks == 1:
   - event_id, total_energy_keV, n_voxels
   - z_length_mm, x_extent_mm, y_extent_mm
   - avg_occupancy, n_time_slots
   - avg_n_pixels_above_threshold (avg number of pixels with signal > 10 ions per time slot)
   - avg_max_signal (avg of max signal per time slot)
   - avg_min_signal (avg of min signal above threshold per time slot)
   - avg_mean_signal (avg of mean signal above threshold per time slot)

Returns tuple of (event_stats_path, single_track_stats_path).
"""
function save_statistics_csv(stats::Vector{EventStats}, results_dir::String)
    # Build event_stats DataFrame (all events)
    event_stats_df = DataFrame(
        event_id = [s.event_id for s in stats],
        n_hits = [s.n_hits for s in stats],
        n_voxels = [s.n_voxels for s in stats],
        n_tracks = [s.n_tracks for s in stats],
        total_energy_keV = [round(s.total_energy_keV, digits=2) for s in stats],
        main_track_energy_keV = [round(s.main_track_energy_keV, digits=2) for s in stats],
        main_track_voxels = [s.main_track_voxels for s in stats],
        track2_energy_keV = [round(s.track2_energy_keV, digits=2) for s in stats],
        track2_voxels = [s.track2_voxels for s in stats],
        track3plus_energy_keV = [round(s.track3plus_energy_keV, digits=2) for s in stats],
        track1_track2_distance_mm = [
            isempty(s.satellite_distances) ? 0.0 : round(s.satellite_distances[1], digits=2)
            for s in stats
        ]
    )

    # Save event_stats.csv
    event_stats_path = joinpath(results_dir, "event_stats.csv")
    CSV.write(event_stats_path, event_stats_df)
    println("Saved event statistics to: $event_stats_path")

    # Build single_track_stats DataFrame (only n_tracks == 1)
    single_track_events = filter(s -> s.n_tracks == 1, stats)

    single_track_df = DataFrame(
        event_id = [s.event_id for s in single_track_events],
        total_energy_keV = [round(s.total_energy_keV, digits=2) for s in single_track_events],
        n_voxels = [s.n_voxels for s in single_track_events],
        z_length_mm = [round(s.z_length, digits=2) for s in single_track_events],
        x_extent_mm = [round(s.x_extent, digits=2) for s in single_track_events],
        y_extent_mm = [round(s.y_extent, digits=2) for s in single_track_events],
        avg_occupancy = [round(s.avg_occupancy, digits=4) for s in single_track_events],
        n_time_slots = [s.n_time_slots for s in single_track_events],
        avg_n_pixels_above_threshold = [round(s.avg_n_pixels_above_threshold, digits=2) for s in single_track_events],
        avg_max_signal = [round(s.avg_max_signal, digits=1) for s in single_track_events],
        avg_min_signal = [round(s.avg_min_signal, digits=1) for s in single_track_events],
        avg_mean_signal = [round(s.avg_mean_signal, digits=1) for s in single_track_events]
    )

    # Save single_track_stats.csv
    single_track_path = joinpath(results_dir, "single_track_stats.csv")
    CSV.write(single_track_path, single_track_df)
    println("Saved single-track statistics to: $single_track_path ($(nrow(single_track_df)) events)")

    return (event_stats_path, single_track_path)
end

"""
    save_plots(plt1, plt2, plt3, results_dir::String)

Save plots as PNG files:
- statistics_summary.png: Track multiplicity, total energy, pixel occupancy, track 2 energy (Frame 1)
- single_track_analysis.png: Z-length and bounding box 2D histogram (Frame 2)
- sample_occupancy_heatmaps.png: Ion occupancy heatmaps from sample event (Frame 3, if available)
"""
function save_plots(plt1, plt2, plt3, results_dir::String)
    # Save main statistics plot (Frame 1)
    plot1_path = joinpath(results_dir, "statistics_summary.png")
    savefig(plt1, plot1_path)
    println("Saved plot to: $plot1_path")

    # Save single-track analysis plot (Frame 2: z-length and bounding box)
    plot2_path = joinpath(results_dir, "single_track_analysis.png")
    savefig(plt2, plot2_path)
    println("Saved plot to: $plot2_path")

    # Save sample occupancy heatmaps (Frame 3, if available)
    plot3_path = nothing
    if plt3 !== nothing
        plot3_path = joinpath(results_dir, "sample_occupancy_heatmaps.png")
        savefig(plt3, plot3_path)
        println("Saved plot to: $plot3_path")
    end

    return plot1_path, plot2_path, plot3_path
end

#=============================================================================
# Main Event Loop
=============================================================================#

"""
    process_events(cmdir, input_file; kwargs...)

Main function to process events.
- batch=false: Interactive mode with plots and user prompts
- batch=true: Batch mode, collects statistics and shows summary at end
- save_results=true: Save plots and statistics to results directory
- readdir=false: If true, read directory of part files and chain them
- method="kNN": Track building method ("kNN" or "KDT")
- k=10: Number of neighbors for kNN method
- max_distance: Max distance (default: Inf for kNN, 2.5*voxel_size for KDT)
- sigmat=0.0: Transverse diffusion sigma in mm (0 = no diffusion)
- sigmal=0.1: Longitudinal diffusion sigma in mm (only used if sigmat > 0)
- grid_size=160.0: Grid size for occupancy histograms in mm (nbins = grid_size / pixel_size)
"""
function process_events(cmdir::String, input_file::String;
                        voxel_size::Float64,
                        method::String="kNN",
                        max_distance::Float64=0.0,
                        k::Int=10,
                        ievt::Int=1,
                        levt::Int=-1,
                        energy_threshold::Float64=0.0,
                        batch::Bool=false,
                        readdir::Bool=false,
                        save_results::Bool=true,
                        sigmat::Float64=0.0,
                        sigmal::Float64=0.1,
                        grid_size::Float64=160.0)

    # Set default max_distance if not provided
    # For kNN: default is Inf (no cutoff), for KDT: default is 2.5 * voxel_size
    if max_distance <= 0.0
        max_distance = method == "kNN" ? Inf : 2.5 * voxel_size
    end

    # Create results directory if saving results (both batch and interactive modes)
    results_dir = nothing
    if save_results
        results_dir = create_results_dir(input_file, voxel_size, max_distance;
                                         sigmat=sigmat, sigmal=sigmal,
                                         energy_threshold=energy_threshold,
                                         ievt=ievt, levt=levt)
    end

    # Load data
    mode_str = batch ? "Batch Mode" : "Interactive Mode"
    println("\n" * "=" ^ 70)
    println("ITACA SIMPLE RECONSTRUCTION - $mode_str")
    println("=" ^ 70)

    if readdir
        # Read directory of part files and chain them
        dirpath = joinpath(cmdir, input_file)
        println("Loading data from directory: $dirpath")
        hitsdf = load_hits_from_directory(dirpath)
    else
        # Read single file
        println("Loading data from: $cmdir/$input_file")
        filepath = joinpath(cmdir, input_file)
        if !isfile(filepath)
            error("File not found: $filepath")
        end
        dfs = Petit.get_dataset_dfs(filepath)
        hitsdf = dfs["hits"]
    end

    # Get unique event IDs (handle both 'event' and 'event_id' column names)
    event_col = :event in names(hitsdf) ? :event : :event_id
    all_event_ids = sort(unique(hitsdf[!, event_col]))
    n_total_events = length(all_event_ids)

    println("Total events in file: $n_total_events")

    # Determine event range
    first_event = max(1, ievt)
    last_event = levt < 0 ? n_total_events : min(levt, n_total_events)

    events_to_process = all_event_ids[first_event:last_event]
    n_events = length(events_to_process)

    println("Processing events $first_event to $last_event ($n_events events)")
    println("\nParameters:")
    if sigmat > 0.0
        println("  Diffusion:         ENABLED")
        println("    sigma_t:         $(round(sigmat, digits=3)) mm")
        println("    sigma_l:         $(round(sigmal, digits=3)) mm")
    else
        println("  Diffusion:         DISABLED")
    end
    println("  Voxel size:        $(round(voxel_size, digits=3)) mm")
    println("  Track method:      $method")
    if method == "kNN"
        println("  k (neighbors):     $k")
    end
    if isinf(max_distance)
        println("  Max distance:      Inf (no cutoff)")
    else
        println("  Max distance:      $(round(max_distance, digits=3)) mm")
    end
    println("  Energy threshold:  $(round(energy_threshold, digits=3)) keV")
    println("  Mode:              $mode_str")
    println("=" ^ 70)

    # Statistics collection for batch mode
    all_stats = EventStats[]
    sample_occ = nothing  # Store first single-track occupancy for batch mode heatmaps
    sample_occ_event_id = 0

    if !batch
        wait_for_enter("\nPress Enter to start processing...")
    else
        println("\nProcessing...")
    end

    # Process each event
    for (idx, event_id) in enumerate(events_to_process)

        if !batch
            println("\n" * "=" ^ 70)
            println("EVENT $event_id  ($idx of $n_events)")
            println("=" ^ 70)
        elseif idx % 100 == 0 || idx == n_events
            print("\r  Processing event $idx / $n_events...")
        end

        # Get event hits
        event_df = Petit.get_event(hitsdf, event_id)
        n_hits = nrow(event_df)

        if n_hits == 0
            if !batch
                println("No hits for event $event_id, skipping...")
            end
            continue
        end

        total_energy = 1e+3 * sum(event_df.energy)  # Convert to keV

        if !batch
            println("Number of MC hits: $n_hits")
            println("Total energy: $(round(total_energy, digits=2)) keV")

            # Plot MC hits
            println("\nPlotting raw MC hits...")
            p1 = Petit.plot_hits(event_df; nbins=100)
            title!(p1[1], "Raw MC Hits - Event $event_id")
            display(p1)
            choice = wait_for_input("[Enter]=continue, [s]=save, [q]=quit: ")
            if choice == :save && results_dir !== nothing
                save_interactive_plot(p1, results_dir, "raw_hits_evt$(event_id)")
            elseif choice == :quit
                println("Exiting...")
                return all_stats
            end
        end

        # Apply diffusion if sigmat > 0
        if sigmat > 0.0
            # Diffuse hits using Monte Carlo method
            diffused_df = Petit.diffuse_xyz_image_mc(event_df;
                                                      sigma_t_mm=sigmat,
                                                      sigma_l_mm=sigmal,
                                                      nbins=300,
                                                      nsigma=3.0)
            n_diffused = nrow(diffused_df)

            if !batch
                println("\nApplying diffusion (σt=$(round(sigmat, digits=3)) mm, σl=$(round(sigmal, digits=3)) mm)...")
                println("Number of diffused voxels: $n_diffused")

                # Plot diffused hits
                println("\nPlotting diffused hits...")
                p_diff = Petit.plot_hits(diffused_df; nbins=100)
                title!(p_diff[1], "Diffused Hits - Event $event_id (σt=$(round(sigmat, digits=2)), σl=$(round(sigmal, digits=2)) mm)")
                display(p_diff)
                choice = wait_for_input("[Enter]=continue, [s]=save, [q]=quit: ")
                if choice == :save && results_dir !== nothing
                    save_interactive_plot(p_diff, results_dir, "diffused_hits_evt$(event_id)")
                elseif choice == :quit
                    println("Exiting...")
                    return all_stats
                end
            end

            # Voxelize diffused hits
            voxels = Petit.voxelize_event(diffused_df, voxel_size)
        else
            # No diffusion - voxelize raw hits directly
            voxels = Petit.voxelize_event(event_df, voxel_size)
        end
        n_voxels = nrow(voxels)

        # Compute voxel nearest neighbor distances (before building tracks)
        voxel_min_dists = Float64[]
        if n_voxels > 1
            voxel_min_dists = Petit.closest_voxel_distances_fast(voxels)
        end

        if !batch
            println("\nVoxelizing event with voxel size = $(round(voxel_size, digits=3)) mm...")
            println("Number of voxels: $n_voxels")

            # Plot voxelized event
            println("\nPlotting voxelized event...")
            p2 = Petit.plot_hits(voxels; nbins=100)
            vox_title = sigmat > 0.0 ? "Voxelized Diffused Hits" : "Voxelized Raw Hits"
            title!(p2[1], "$vox_title - Event $event_id (voxel=$(round(voxel_size, digits=2)) mm)")
            display(p2)
            choice = wait_for_input("[Enter]=continue, [s]=save, [q]=quit: ")
            if choice == :save && results_dir !== nothing
                save_interactive_plot(p2, results_dir, "voxelized_evt$(event_id)")
            elseif choice == :quit
                println("Exiting...")
                return all_stats
            end

            # Plot voxel distance distribution
            if !isempty(voxel_min_dists)
                println("\nPlotting voxel nearest neighbor distances...")
                p_dist = histogram(voxel_min_dists, bins=30,
                                   xlabel="Closest Neighbor Distance (mm)", ylabel="Count",
                                   title="Voxel NN Distances - Event $event_id",
                                   legend=false, color=:purple)
                display(p_dist)
                choice = wait_for_input("[Enter]=continue, [s]=save, [q]=quit: ")
                if choice == :save && results_dir !== nothing
                    save_interactive_plot(p_dist, results_dir, "voxel_distances_evt$(event_id)")
                elseif choice == :quit
                    println("Exiting...")
                    return all_stats
                end
            end
        end

        # Build tracks
        if method == "kNN"
            tracks = Petit.make_tracks(voxels;
                                       method="kNN",
                                       k=k,
                                       max_distance_mm=max_distance,
                                       energy_threshold_kev=energy_threshold)
        else  # KDT (radius graph)
            tracks = Petit.make_tracks(voxels;
                                       method="KDT",
                                       max_distance_mm=max_distance,
                                       energy_threshold_kev=energy_threshold)
        end

        n_tracks = length(tracks)

        # Compute statistics
        main_track_energy = n_tracks > 0 ? 1e+3 * sum(tracks[1].voxels.energy) : 0.0
        main_track_voxels = n_tracks > 0 ? nrow(tracks[1].voxels) : 0

        # Track 2 energy and voxels
        track2_energy = n_tracks >= 2 ? 1e+3 * sum(tracks[2].voxels.energy) : 0.0
        track2_voxels = n_tracks >= 2 ? nrow(tracks[2].voxels) : 0

        # Track 3+ combined energy
        track3plus_energy = 0.0
        if n_tracks >= 3
            for i in 3:n_tracks
                track3plus_energy += 1e+3 * sum(tracks[i].voxels.energy)
            end
        end

        distances = compute_track_distances(tracks)

        # Compute z-length, bounding box, and occupancy analysis for single-track events
        z_length = 0.0
        x_extent = 0.0
        y_extent = 0.0
        avg_occupancy = 0.0
        n_time_slots = 0
        occupancy_per_slot = Float64[]
        occ = nothing  # Store full occupancy result for plotting
        # New pixel signal statistics
        avg_n_pixels_above_threshold = 0.0
        avg_max_signal = 0.0
        avg_min_signal = 0.0
        avg_mean_signal = 0.0

        if n_tracks == 1
            z_ext = compute_track_z_extension(tracks[1])
            z_length = z_ext.z_length

            # Compute bounding box
            bbox = compute_track_bbox(tracks[1])
            x_extent = bbox.x_extent
            y_extent = bbox.y_extent

            # Use max_distance as pixel size (they should match for consistent analysis)
            pixel_size_for_occ = isinf(max_distance) ? voxel_size : max_distance
            # Compute nbins from grid_size / pixel_size
            nbins = round(Int, grid_size / pixel_size_for_occ)
            occ = compute_occupancy_analysis(tracks[1]; pixel_size_mm=pixel_size_for_occ, nbins=nbins)
            avg_occupancy = occ.avg_occupancy
            n_time_slots = occ.n_time_slots
            occupancy_per_slot = occ.occupancy_per_slot
            # Extract new pixel signal statistics
            avg_n_pixels_above_threshold = occ.avg_n_pixels_above_threshold
            avg_max_signal = occ.avg_max_signal
            avg_min_signal = occ.avg_min_signal
            avg_mean_signal = occ.avg_mean_signal

            # Store first single-track occupancy for batch mode heatmaps
            if sample_occ === nothing
                sample_occ = occ
                sample_occ_event_id = event_id
            end
        end

        # Store stats
        push!(all_stats, EventStats(
            event_id, n_hits, n_voxels, n_tracks,
            total_energy, main_track_energy, main_track_voxels,
            track2_energy, track2_voxels, track3plus_energy,
            distances.satellite_distances,
            voxel_min_dists,
            z_length, x_extent, y_extent, avg_occupancy, n_time_slots, occupancy_per_slot,
            avg_n_pixels_above_threshold, avg_max_signal, avg_min_signal, avg_mean_signal
        ))

        if !batch
            println("\nBuilding tracks...")
            println("  Method: $method")
            if method == "kNN"
                println("  k (neighbors): $k")
            end
            if isinf(max_distance)
                println("  Max distance: Inf (no cutoff)")
            else
                println("  Max distance: $(round(max_distance, digits=3)) mm")
            end
            println("  Energy threshold: $(round(energy_threshold, digits=3)) keV")
            println("Number of tracks found: $n_tracks")

            # For single-track events, print z-length and occupancy analysis
            if n_tracks == 1
                println()
                println("─" ^ 50)
                println("SINGLE-TRACK ANALYSIS")
                println("─" ^ 50)
                println("  Z-extension (longitudinal): $(round(z_length, digits=2)) mm")
                println("  Number of time slots (2mm): $n_time_slots")
                println("  Average pixel occupancy:    $(round(100 * avg_occupancy, digits=2))%")
                if !isempty(occupancy_per_slot) && length(occupancy_per_slot) >= 3
                    println("  Occupancy range: $(round(100 * minimum(occupancy_per_slot), digits=2))% - $(round(100 * maximum(occupancy_per_slot), digits=2))%")
                end
                println("─" ^ 50)

                # Plot occupancy heatmaps
                if occ !== nothing
                    println("\nPlotting pixel occupancy heatmaps...")
                    p_occ = plot_occupancy_heatmaps(occ; event_id=event_id)
                    display(p_occ)
                    choice = wait_for_input("[Enter]=continue, [s]=save, [q]=quit: ")
                    if choice == :save && results_dir !== nothing
                        save_interactive_plot(p_occ, results_dir, "occupancy_heatmaps_evt$(event_id)")
                    elseif choice == :quit
                        println("Exiting...")
                        break
                    end
                end
            end

            # Print track summary (waits for Enter internally)
            print_track_summary(tracks, event_id)

            # Plot all tracks in a single figure
            if n_tracks > 0
                println("\nPlotting all tracks (T1=viridis, T2=blue, T3=green, T4+=black)...")
                p3 = plot_all_tracks(tracks, event_id)
                display(p3)
                choice = wait_for_input("[Enter]=continue, [s]=save, [q]=quit: ")
                if choice == :save && results_dir !== nothing
                    save_interactive_plot(p3, results_dir, "tracks_evt$(event_id)")
                elseif choice == :quit
                    println("Exiting...")
                    break
                end

                # Ask if user wants to continue to next event
                if idx < n_events
                    print("\nContinue to next event? (Enter=yes, q=quit): ")
                    response = readline()
                    if lowercase(strip(response)) == "q"
                        println("\nExiting...")
                        break
                    end
                end
            end
        end
    end

    if batch
        println("\n")  # Clear the progress line

        # Print statistics summary
        print_statistics_summary(all_stats)

        # Compute pixel size for occupancy analysis
        pixel_size_for_occ = isinf(max_distance) ? voxel_size : max_distance

        # Plot statistics (three frames)
        println("\nGenerating summary plots...")
        plt1, plt2, _ = plot_statistics(all_stats; grid_size=grid_size, pixel_size_mm=pixel_size_for_occ)

        # Generate sample occupancy heatmaps (Frame 3) if we have single-track data
        plt3 = nothing
        if sample_occ !== nothing
            println("Generating sample occupancy heatmaps (Event $sample_occ_event_id)...")
            plt3 = plot_occupancy_heatmaps(sample_occ; event_id=sample_occ_event_id)
        end

        # Save results if requested
        if save_results && results_dir !== nothing
            println("\nSaving results...")
            params = Dict(
                "voxel_size_mm" => voxel_size,
                "sigmat_mm" => sigmat,
                "sigmal_mm" => sigmal,
                "diffusion_enabled" => sigmat > 0.0,
                "method" => method,
                "k" => method == "kNN" ? k : nothing,
                "max_distance_mm" => isinf(max_distance) ? "Inf" : max_distance,
                "energy_threshold_keV" => energy_threshold,
                "grid_size_mm" => grid_size,
                "ievt" => ievt,
                "levt" => levt,
                "input_file" => input_file,
                "cmdir" => cmdir
            )
            nbins_for_save = round(Int, grid_size / pixel_size_for_occ)
            save_statistics_json(all_stats, results_dir, params; nbins=nbins_for_save, pixel_size_mm=pixel_size_for_occ, grid_size=grid_size)
            save_statistics_csv(all_stats, results_dir)
            save_plots(plt1, plt2, plt3, results_dir)
        end

        # Display plots interactively (skip in batch mode)
        if !batch
            display(plt1)
            wait_for_enter("Press Enter for single-track analysis...")

            display(plt2)

            if plt3 !== nothing
                wait_for_enter("Press Enter for sample occupancy heatmaps...")
                display(plt3)
            end

            wait_for_enter("Press Enter to exit...")
        end
    end

    println("\n" * "=" ^ 70)
    println("PROCESSING COMPLETE")
    println("=" ^ 70)

    return all_stats
end

#=============================================================================
# Command Line Interface
=============================================================================#

function main()
    # Check minimum required arguments
    if length(ARGS) < 2
        println("Error: Missing required arguments")
        println()
        println("Usage: julia itaca_simple_reco.jl <cmdir> <input_file> --voxel-size=X [options]")
        println()
        println("Required arguments:")
        println("  cmdir           Directory containing the input file")
        println("  input_file      Name of the HDF5 input file")
        println("  --voxel-size=X  Voxel size in mm (required)")
        println()
        println("Optional arguments:")
        println("  --sigmat=X            Transverse diffusion sigma in mm (default: 0, no diffusion)")
        println("  --sigmal=X            Longitudinal diffusion sigma in mm (default: 0.1, only if sigmat > 0)")
        println("  --method=X            Track building method: 'kNN' (default) or 'KDT'")
        println("  --k=N                 Number of neighbors for kNN method (default: 10)")
        println("  --max-distance=X      Max distance in mm (default: Inf for kNN, 2.5*voxel_size for KDT)")
        println("  --ievt=N              First event to process (default: 1)")
        println("  --levt=N              Last event to process (default: -1, all)")
        println("  --energy-threshold=X  Energy threshold in keV (default: 1.0)")
        println("  --grid-size=X         Grid size for occupancy histograms in mm (default: 160)")
        println("  --batch               Run in batch mode (no interactive plots)")
        println("  --readdir             Read directory of part files (*_part_*.h5) and chain them")
        println("  --save-results=X      Save results to <filename>_results/ (default: true)")
        println()
        println("Examples:")
        println("  # Interactive mode with kNN (default method):")
        println("  julia itaca_simple_reco.jl /data/itaca/ bb0nu.h5 --voxel-size=2.0 --ievt=1 --levt=10")
        println()
        println("  # Batch mode with kNN (default):")
        println("  julia itaca_simple_reco.jl /data/itaca/ bb0nu.h5 --voxel-size=2.0 --batch")
        println()
        println("  # Batch mode with KDT (radius graph):")
        println("  julia itaca_simple_reco.jl /data/itaca/ bb0nu.h5 --voxel-size=2.0 --method=KDT --batch")
        println()
        println("  # kNN with custom k:")
        println("  julia itaca_simple_reco.jl /data/itaca/ bb0nu.h5 --voxel-size=2.0 --k=15 --batch")
        println()
        println("  # kNN with distance cutoff:")
        println("  julia itaca_simple_reco.jl /data/itaca/ bb0nu.h5 --voxel-size=2.0 --k=10 --max-distance=5.0 --batch")
        println()
        println("  # Read from directory of part files:")
        println("  julia itaca_simple_reco.jl /data/itaca/ filtered_ecut_2400_2500 --voxel-size=2.0 --readdir --batch")
        println()
        println("  # With diffusion (σt=1.0 mm, σl=0.5 mm):")
        println("  julia itaca_simple_reco.jl /data/itaca/ bb0nu.h5 --voxel-size=2.0 --sigmat=1.0 --sigmal=0.5 --batch")
        exit(1)
    end

    # Parse required arguments
    cmdir = ARGS[1]
    input_file = ARGS[2]

    # Default values
    voxel_size = nothing  # Required
    method = "kNN"        # Default: kNN (no max_distance dependency)
    k = 10                # Default k for kNN
    max_distance = 0.0    # For KDT method
    ievt = 1
    levt = -1
    energy_threshold = 0.0
    batch = false
    readdir = false       # Read directory of part files
    save_results = true   # Save results by default
    sigmat = 0.0          # Transverse diffusion (0 = no diffusion)
    sigmal = 0.1          # Longitudinal diffusion (only used if sigmat > 0)
    grid_size = 160.0     # Grid size for occupancy histograms in mm

    # Parse optional arguments
    for arg in ARGS[3:end]
        if arg == "--batch"
            batch = true
        elseif arg == "--readdir"
            readdir = true
        elseif startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) != 2
                println("Warning: Ignoring malformed argument: $arg")
                continue
            end
            key, value = parts

            try
                if key == "voxel-size"
                    voxel_size = parse(Float64, value)
                elseif key == "method"
                    method = String(value)
                    if !(method in ("kNN", "KDT"))
                        println("Error: --method must be 'kNN' or 'KDT'")
                        exit(1)
                    end
                elseif key == "k"
                    k = parse(Int, value)
                elseif key == "max-distance"
                    max_distance = parse(Float64, value)
                elseif key == "ievt"
                    ievt = parse(Int, value)
                elseif key == "levt"
                    levt = parse(Int, value)
                elseif key == "energy-threshold"
                    energy_threshold = parse(Float64, value)
                elseif key == "save-results"
                    save_results = lowercase(value) in ("true", "1", "yes")
                elseif key == "readdir"
                    readdir = lowercase(value) in ("true", "1", "yes")
                elseif key == "sigmat"
                    sigmat = parse(Float64, value)
                elseif key == "sigmal"
                    sigmal = parse(Float64, value)
                elseif key == "grid-size"
                    grid_size = parse(Float64, value)
                else
                    println("Warning: Unknown argument: --$key")
                end
            catch e
                println("Error parsing argument --$key=$value: $e")
                exit(1)
            end
        end
    end

    # Check required voxel-size
    if isnothing(voxel_size)
        println("Error: --voxel-size is required")
        exit(1)
    end

    # Run
    process_events(cmdir, input_file;
                   voxel_size=voxel_size,
                   method=method,
                   k=k,
                   max_distance=max_distance,
                   ievt=ievt,
                   levt=levt,
                   energy_threshold=energy_threshold,
                   batch=batch,
                   readdir=readdir,
                   save_results=save_results,
                   sigmat=sigmat,
                   sigmal=sigmal,
                   grid_size=grid_size)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
