"""
Kernel Density Estimation functions for track analysis.

Provides energy-weighted KDE along the longitudinal coordinate of a track
to identify Bragg peaks and estimate dE/dx profiles.
"""

using DataFrames
using Statistics
using LinearAlgebra
using KernelDensity
using Peaks

export compute_arc_length, project_voxel_to_path, project_voxels_to_path
export energy_weighted_kde, find_maxima, find_peaks, plot_kde_peaks

"""
    compute_arc_length(central_path::DataFrame)

Compute cumulative arc length s along the central path.

# Arguments
- `central_path::DataFrame`: Path with x, y, z columns

# Returns
- `Vector{Float64}`: Cumulative arc length at each path point, starting at 0
"""
function compute_arc_length(central_path::DataFrame)
    n = nrow(central_path)
    s = zeros(n)
    for i in 2:n
        dx = central_path.x[i] - central_path.x[i-1]
        dy = central_path.y[i] - central_path.y[i-1]
        dz = central_path.z[i] - central_path.z[i-1]
        s[i] = s[i-1] + sqrt(dx^2 + dy^2 + dz^2)
    end
    return s
end

"""
    project_voxel_to_path(voxel_pos, central_path, s_values)

Project a single voxel position onto the central path.

# Arguments
- `voxel_pos::Tuple{Float64,Float64,Float64}`: (x, y, z) position of voxel
- `central_path::DataFrame`: Path with x, y, z columns
- `s_values::Vector{Float64}`: Arc length values from `compute_arc_length`

# Returns
- `Float64`: Arc-length coordinate s of the closest point on the path
"""
function project_voxel_to_path(voxel_pos::Tuple{Float64,Float64,Float64},
                                central_path::DataFrame,
                                s_values::Vector{Float64})
    min_dist = Inf
    closest_idx = 1

    for i in 1:nrow(central_path)
        dx = voxel_pos[1] - central_path.x[i]
        dy = voxel_pos[2] - central_path.y[i]
        dz = voxel_pos[3] - central_path.z[i]
        dist = sqrt(dx^2 + dy^2 + dz^2)
        if dist < min_dist
            min_dist = dist
            closest_idx = i
        end
    end

    return s_values[closest_idx]
end

"""
    project_voxels_to_path(voxels::DataFrame, central_path::DataFrame, s_values::Vector{Float64})

Project all voxels onto the central path.

# Arguments
- `voxels::DataFrame`: Voxels with x, y, z, energy columns
- `central_path::DataFrame`: Path with x, y, z columns
- `s_values::Vector{Float64}`: Arc length values from `compute_arc_length`

# Returns
- `s_voxels::Vector{Float64}`: Arc-length coordinate for each voxel
- `energies::Vector{Float64}`: Energy of each voxel
"""
function project_voxels_to_path(voxels::DataFrame,
                                 central_path::DataFrame,
                                 s_values::Vector{Float64})
    n_voxels = nrow(voxels)
    s_voxels = Vector{Float64}(undef, n_voxels)
    energies = Vector{Float64}(undef, n_voxels)

    for i in 1:n_voxels
        pos = (voxels.x[i], voxels.y[i], voxels.z[i])
        s_voxels[i] = project_voxel_to_path(pos, central_path, s_values)
        energies[i] = voxels.energy[i]
    end

    return s_voxels, energies
end

"""
    energy_weighted_kde(s_samples, energies, s_eval; bandwidth=5.0)

Compute energy-weighted Kernel Density Estimation using KernelDensity.jl.

Estimates the longitudinal energy density:
    f(s) ∝ Σ E_i K((s - s_i)/h)

where K is a Gaussian kernel and h is the bandwidth.

# Arguments
- `s_samples::Vector{Float64}`: Arc-length positions of voxels
- `energies::Vector{Float64}`: Energy of each voxel
- `s_eval::AbstractVector`: Points at which to evaluate the KDE
- `bandwidth::Float64=5.0`: Kernel bandwidth in mm (default: 5.0)

# Returns
- `f::Vector{Float64}`: Energy density at each evaluation point
- `bandwidth::Float64`: The bandwidth used
"""
function energy_weighted_kde(s_samples::Vector{Float64},
                             energies::Vector{Float64},
                             s_eval::AbstractVector;
                             bandwidth::Float64=5.0)
    # Use KernelDensity with energy weights
    kde_obj = kde(s_samples; weights=energies, bandwidth=bandwidth)

    # Evaluate at requested points
    f = pdf(kde_obj, collect(s_eval))

    return f, bandwidth
end

"""
    find_maxima(x, y; min_prominence=0.1, boundary_width=5)

Find local maxima in y(x) with minimum prominence filtering.
Also checks boundary regions (important for Bragg peaks at track ends).

# Arguments
- `x::AbstractVector`: x coordinates
- `y::Vector{Float64}`: y values
- `min_prominence::Float64=0.1`: Minimum prominence as fraction of y range
- `boundary_width::Int=5`: Number of points to consider in boundary regions

# Returns
- `maxima_idx::Vector{Int}`: Indices of maxima, sorted by y value (descending)
- `prominences::Vector{Float64}`: Prominence of each maximum (same order)
"""
function find_maxima(x::AbstractVector, y::Vector{Float64};
                     min_prominence::Float64=0.1,
                     boundary_width::Int=5)
    n = length(y)
    maxima_idx = Int[]
    prominences = Float64[]

    y_range = maximum(y) - minimum(y)

    # Use small boundary region (just a few points near endpoints)
    boundary_region = min(boundary_width, 3)  # max 3 points

    y_max = maximum(y)
    y_min = minimum(y)

    # Check for maximum in left boundary region
    # Only accept if it's a significant fraction of the global max
    if boundary_region >= 1
        left_max_idx = argmax(y[1:boundary_region])
        left_max_val = y[left_max_idx]
        prominence = left_max_val - y_min
        # Must be prominent AND significant (>50% of global max height)
        if prominence > min_prominence * y_range && left_max_val > 0.5 * y_max
            push!(maxima_idx, left_max_idx)
            push!(prominences, prominence)
        end
    end

    # Check interior points (outside boundary regions)
    for i in (boundary_region + 1):(n - boundary_region)
        if y[i] > y[i-1] && y[i] > y[i+1]
            # Check prominence
            left_min = minimum(y[1:i])
            right_min = minimum(y[i:n])
            prominence = y[i] - max(left_min, right_min)

            if prominence > min_prominence * y_range
                push!(maxima_idx, i)
                push!(prominences, prominence)
            end
        end
    end

    # Check for maximum in right boundary region
    # Only accept if it's a significant fraction of the global max
    if boundary_region >= 1
        right_start = n - boundary_region + 1
        right_max_local = argmax(y[right_start:n])
        right_max_idx = right_start + right_max_local - 1
        right_max_val = y[right_max_idx]
        # Don't add if already found
        if !(right_max_idx in maxima_idx)
            prominence = right_max_val - y_min
            # Must be prominent AND significant (>50% of global max height)
            if prominence > min_prominence * y_range && right_max_val > 0.5 * y_max
                push!(maxima_idx, right_max_idx)
                push!(prominences, prominence)
            end
        end
    end

    # Sort by y value (descending)
    sort_order = sortperm(maxima_idx, by=i -> y[i], rev=true)
    maxima_idx = maxima_idx[sort_order]
    prominences = prominences[sort_order]

    return maxima_idx, prominences
end

"""
    analyze_track_energy_profile(track, central_path; bandwidth=5.0, n_eval=200, min_prominence=0.1)

Analyze the longitudinal energy profile of a track using KDE.

# Arguments
- `track`: Track object with voxels DataFrame
- `central_path::DataFrame`: Smoothed central path
- `bandwidth::Float64=5.0`: KDE bandwidth in mm
- `n_eval::Int=200`: Number of evaluation points
- `min_prominence::Float64=0.1`: Minimum prominence for peak detection (fraction of y range)

# Returns
Named tuple with:
- `s_eval`: Evaluation points (arc length)
- `f_kde`: Energy density values
- `maxima_idx`: Indices of detected maxima
- `maxima_s`: Arc-length positions of maxima
- `maxima_f`: Energy density values at maxima
- `prominences`: Prominence of each maximum
- `bandwidth`: Bandwidth used
- `track_length`: Total arc length
"""
function analyze_track_energy_profile(track, central_path::DataFrame;
                                      bandwidth::Float64=5.0,
                                      n_eval::Int=200,
                                      min_prominence::Float64=0.1)
    # Compute arc length along path
    s_path = compute_arc_length(central_path)
    track_length = s_path[end]

    # Project voxels onto path
    s_voxels, energies = project_voxels_to_path(track.voxels, central_path, s_path)

    # Compute KDE
    s_eval = range(0, track_length, length=n_eval)
    f_kde, h = energy_weighted_kde(s_voxels, energies, collect(s_eval); bandwidth=bandwidth)

    # Find maxima with prominences
    maxima_idx, prominences = find_maxima(collect(s_eval), f_kde; min_prominence=min_prominence)

    # Extract maxima positions and values
    maxima_s = [s_eval[i] for i in maxima_idx]
    maxima_f = [f_kde[i] for i in maxima_idx]

    return (
        s_eval = collect(s_eval),
        f_kde = f_kde,
        maxima_idx = maxima_idx,
        maxima_s = maxima_s,
        maxima_f = maxima_f,
        prominences = prominences,
        bandwidth = h,
        track_length = track_length,
        s_voxels = s_voxels,
        energies = energies
    )
end

"""
    get_reco_kde(track, path; bandwidth=5.0, n_eval=200)

Compute KDE for reconstructed track voxels projected onto a path.

# Arguments
- `track`: Track object with voxels DataFrame (must have x, y, z, energy columns)
- `path::DataFrame`: Central path with x, y, z columns
- `bandwidth::Float64=5.0`: KDE bandwidth in mm
- `n_eval::Int=200`: Number of evaluation points

# Returns
NamedTuple with:
- `track_length`: Total arc length of path
- `s`: Arc-length positions of voxels
- `E`: Energy of each voxel
- `kde_s`: Evaluation grid points
- `kde_f`: KDE values at evaluation points
"""
function get_reco_kde(track, path::DataFrame; bandwidth::Float64=5.0, n_eval::Int=200)
    # Compute arc length along path
    s_path = compute_arc_length(path)
    track_length = s_path[end]

    # Project RECO voxels onto path
    reco_s, reco_E = project_voxels_to_path(track.voxels, path, s_path)

    # Create evaluation grid for KDE
    kde_s = collect(range(0.0, track_length, length=n_eval))

    # Compute RECO KDE
    kde_f, _ = energy_weighted_kde(reco_s, reco_E, kde_s; bandwidth=bandwidth)

    (track_length = track_length, s = reco_s, E = reco_E, kde_s = kde_s, kde_f = kde_f)
end

"""
    get_mc_kde(mc_path; bandwidth=5.0, n_eval=200)

Compute KDE for MC path (which already has s and energy columns).

# Arguments
- `mc_path::DataFrame`: MC path with s and energy columns (from compute_mc_path)
- `bandwidth::Float64=5.0`: KDE bandwidth in mm
- `n_eval::Int=200`: Number of evaluation points

# Returns
NamedTuple with:
- `track_length`: Total arc length of MC path
- `s`: Arc-length positions from mc_path
- `E`: Energy values from mc_path
- `kde_s`: Evaluation grid points
- `kde_f`: KDE values at evaluation points
"""
function get_mc_kde(mc_path::DataFrame; bandwidth::Float64=5.0, n_eval::Int=200)
    mc_s = Vector{Float64}(mc_path.s)
    mc_E = Vector{Float64}(mc_path.energy)
    track_length = mc_path.s[end]

    # Create evaluation grid for MC KDE
    kde_s = collect(range(0.0, track_length, length=n_eval))

    # Compute MC KDE
    kde_f, _ = energy_weighted_kde(mc_s, mc_E, kde_s; bandwidth=bandwidth)

    (track_length = track_length, s = mc_s, E = mc_E, kde_s = kde_s, kde_f = kde_f)
end

"""
    find_peaks(kde_f, kde_s; prom_scale=0.1)

Find peaks in KDE profile using Peaks.jl with prominence filtering.
Pads the array with zeros at boundaries so Peaks.jl naturally detects
boundary peaks (important for Bragg peaks at track ends).

# Arguments
- `kde_f::Vector{Float64}`: KDE values (energy density)
- `kde_s::Vector{Float64}`: Arc-length positions
- `prom_scale::Float64=0.1`: Minimum prominence as fraction of y range

# Returns
NamedTuple with:
- `indices`: Peak indices sorted by prominence (descending)
- `proms`: Prominence values
- `positions`: Peak positions in s coordinates
- `widths`: Peak widths at half prominence
- `lefts`: Left edge of each peak
- `rights`: Right edge of each peak
"""
function find_peaks(kde_f::Vector{Float64}, kde_s::AbstractVector;
                    prom_scale::Float64=0.1)
    n = length(kde_f)

    # Pad with zeros so boundary peaks become interior peaks for Peaks.jl
    kde_f_padded = vcat(0.0, kde_f, 0.0)

    # Find all peaks in padded array
    indices_padded = argmaxima(kde_f_padded)

    if isempty(indices_padded)
        return (indices = Int[], proms = Float64[], positions = Float64[],
                widths = Float64[], lefts = Float64[], rights = Float64[])
    end

    # Get prominences with minimum filter (work in padded space)
    y_range = maximum(kde_f) - minimum(kde_f)
    min_prom = prom_scale * y_range

    indices_padded, proms = peakproms(indices_padded, kde_f_padded; min=min_prom)

    if isempty(indices_padded)
        return (indices = Int[], proms = Float64[], positions = Float64[],
                widths = Float64[], lefts = Float64[], rights = Float64[])
    end

    # Get widths (in padded space)
    indices_padded, widths, left_padded, right_padded = peakwidths(indices_padded, kde_f_padded, proms)

    # Convert from padded indices back to original indices (subtract 1)
    indices = indices_padded .- 1
    left = left_padded .- 1.0
    right = right_padded .- 1.0

    # Clamp left/right to valid range in original coordinates
    left = clamp.(left, 1.0, Float64(n))
    right = clamp.(right, 1.0, Float64(n))

    # Convert indices to s positions
    peak_s = collect(kde_s)[indices]

    # Sort by prominence (descending) to get Bragg peaks first
    order = sortperm(proms, rev=true)

    (indices = indices[order],
     proms = proms[order],
     positions = peak_s[order],
     widths = widths[order],
     lefts = left[order],
     rights = right[order])
end

"""
    plot_kde_peaks(kde_s, kde_f, peaks; kwargs...)

Plot KDE profile with peaks marked.

# Arguments
- `kde_s`: Arc-length positions
- `kde_f`: KDE values
- `peaks`: Result from find_peaks()
- `xlabel`: X-axis label (default: "Arc length s (mm)")
- `ylabel`: Y-axis label (default: "Energy density f(s)")
- `title`: Plot title (default: "KDE Energy Profile")
- `linewidth`: KDE line width (default: 2)
- `color`: KDE line color (default: :blue)

# Returns
- Plot object
"""
function plot_kde_peaks(kde_s::AbstractVector, kde_f::Vector{Float64}, peaks::NamedTuple;
                        xlabel::String="Arc length s (mm)",
                        ylabel::String="Energy density f(s)",
                        title::String="KDE Energy Profile",
                        linewidth::Real=2,
                        color=:blue)
    # Plot KDE curve
    p = plot(kde_s, kde_f,
             xlabel=xlabel, ylabel=ylabel, title=title,
             linewidth=linewidth, color=color, label="KDE")

    # Mark peak maxima with dots
    if !isempty(peaks.positions)
        peak_f = kde_f[peaks.indices]
        scatter!(p, peaks.positions, peak_f,
                 markersize=4, color=:red, label="peaks")
    end

    # Draw vertical dashed lines at left and right edges
    for i in eachindex(peaks.lefts)
        # Left edge
        left_s = collect(kde_s)[max(1, floor(Int, peaks.lefts[i]))]
        vline!(p, [left_s], linestyle=:dash, linewidth=1, color=:gray, alpha=0.7, label=(i==1 ? "width edges" : ""))

        # Right edge
        right_s = collect(kde_s)[min(length(kde_s), ceil(Int, peaks.rights[i]))]
        vline!(p, [right_s], linestyle=:dash, linewidth=1, color=:gray, alpha=0.7, label="")
    end

    return p
end
