using Statistics: mean
using Graphs: nv, edges, src, dst

function plot_hits_trk(trk::Tracks; nbins::Int=100)
	plot_hits(trk.voxels; nbins)
end

function plot_track_with_extremes(track::Tracks, walk_result=nothing;
                                 markersize_voxels::Float64=4.0,
                                 markersize_extremes::Float64=10.0)
    """
    Plot a track showing voxels, path, and extremes in three projections.

    Parameters:
    - track: A Tracks object to visualize
    - walk_result: Result from walk_track_from_extremes (optional)
    - markersize_voxels: Size for voxel markers
    - markersize_extremes: Size for extreme point markers

    Returns:
    - A plot with XY, XZ, YZ projections and 3D view
    """

    # Extract voxel positions and energies
    x = Float64.(track.voxels.x)
    y = Float64.(track.voxels.y)
    z = Float64.(track.voxels.z)
    e = Float64.(track.voxels.energy)

    # Compute plot limits with padding
    xmid, xrange = mean((minimum(x), maximum(x))), maximum(x) - minimum(x)
    ymid, yrange = mean((minimum(y), maximum(y))), maximum(y) - minimum(y)
    zmid, zrange = mean((minimum(z), maximum(z))), maximum(z) - minimum(z)

    # Add 20% padding
    xlim = (xmid - 0.6 * xrange, xmid + 0.6 * xrange)
    ylim = (ymid - 0.6 * yrange, ymid + 0.6 * yrange)
    zlim = (zmid - 0.6 * zrange, zmid + 0.6 * zrange)

    # Color map for energy
    cmap = cgrad(:viridis)

    # === XY Projection ===
    p1 = scatter(x, y, marker_z=e, ms=markersize_voxels,
                xlabel="X (mm)", ylabel="Y (mm)", title="XY Projection",
                xlims=xlim, ylims=ylim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="Voxels",
                markerstrokewidth=0)

    # Add extremes if walk_result provided
    if !isnothing(walk_result) && !isnothing(walk_result.extremes[1])
        start_voxel, end_voxel = walk_result.extremes

        # Mark extremes in red with prominent circles
        scatter!(p1, [start_voxel.x], [start_voxel.y],
                ms=markersize_extremes, color=:red, markershape=:circle,
                label="Start", markerstrokewidth=3, markerstrokecolor=:white)
        scatter!(p1, [end_voxel.x], [end_voxel.y],
                ms=markersize_extremes, color=:darkred, markershape=:circle,
                label="End", markerstrokewidth=3, markerstrokecolor=:white)

        # Draw path if available
        if length(walk_result.path_indices) > 1
            path_x = x[walk_result.path_indices]
            path_y = y[walk_result.path_indices]
            plot!(p1, path_x, path_y, color=:orange, linewidth=2,
                 alpha=0.7, label="Path")
        end
    end

    # === XZ Projection ===
    p2 = scatter(x, z, marker_z=e, ms=markersize_voxels,
                xlabel="X (mm)", ylabel="Z (mm)", title="XZ Projection",
                xlims=xlim, ylims=zlim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="Voxels",
                markerstrokewidth=0)

    if !isnothing(walk_result) && !isnothing(walk_result.extremes[1])
        start_voxel, end_voxel = walk_result.extremes
        scatter!(p2, [start_voxel.x], [start_voxel.z],
                ms=markersize_extremes, color=:red, markershape=:circle,
                label="", markerstrokewidth=3, markerstrokecolor=:white)
        scatter!(p2, [end_voxel.x], [end_voxel.z],
                ms=markersize_extremes, color=:darkred, markershape=:circle,
                label="", markerstrokewidth=3, markerstrokecolor=:white)

        if length(walk_result.path_indices) > 1
            path_x = x[walk_result.path_indices]
            path_z = z[walk_result.path_indices]
            plot!(p2, path_x, path_z, color=:orange, linewidth=2,
                 alpha=0.7, label="")
        end
    end

    # === YZ Projection ===
    p3 = scatter(y, z, marker_z=e, ms=markersize_voxels,
                xlabel="Y (mm)", ylabel="Z (mm)", title="YZ Projection",
                xlims=ylim, ylims=zlim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="Voxels",
                markerstrokewidth=0)

    if !isnothing(walk_result) && !isnothing(walk_result.extremes[1])
        start_voxel, end_voxel = walk_result.extremes
        scatter!(p3, [start_voxel.y], [start_voxel.z],
                ms=markersize_extremes, color=:red, markershape=:circle,
                label="", markerstrokewidth=3, markerstrokecolor=:white)
        scatter!(p3, [end_voxel.y], [end_voxel.z],
                ms=markersize_extremes, color=:darkred, markershape=:circle,
                label="", markerstrokewidth=3, markerstrokecolor=:white)

        if length(walk_result.path_indices) > 1
            path_y = y[walk_result.path_indices]
            path_z = z[walk_result.path_indices]
            plot!(p3, path_y, path_z, color=:orange, linewidth=2,
                 alpha=0.7, label="")
        end
    end

    # === 3D View ===
    p4 = scatter(x, y, z, marker_z=e, ms=markersize_voxels,
                xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)",
                title="3D View",
                xlims=xlim, ylims=ylim, zlims=zlim,
                colorbar_title="Energy (MeV)", color=cmap, label="Voxels",
                markerstrokewidth=0)

    # Add 3D extremes
    if !isnothing(walk_result) && !isnothing(walk_result.extremes[1])
        start_voxel, end_voxel = walk_result.extremes
        scatter!(p4, [start_voxel.x], [start_voxel.y], [start_voxel.z],
                ms=markersize_extremes, color=:red, markershape=:circle,
                label="Start", markerstrokewidth=3, markerstrokecolor=:white)
        scatter!(p4, [end_voxel.x], [end_voxel.y], [end_voxel.z],
                ms=markersize_extremes, color=:darkred, markershape=:circle,
                label="End", markerstrokewidth=3, markerstrokecolor=:white)

        # Draw 3D path
        if length(walk_result.path_indices) > 1
            path_x = x[walk_result.path_indices]
            path_y = y[walk_result.path_indices]
            path_z = z[walk_result.path_indices]
            plot!(p4, path_x, path_y, path_z, color=:orange, linewidth=2,
                 alpha=0.7, label="Path")
        end
    end

    # Add info text if walk_result provided
    title_text = "Track Visualization"
    if !isnothing(walk_result) && !isnothing(walk_result.extremes[1])
        title_text *= "\nPath length: $(round(walk_result.total_length, digits=2)) mm"
        title_text *= " | Voxels: $(length(walk_result.path_indices))"
        title_text *= " | Total energy: $(round(sum(e), digits=3))"
    end

    # Combine all plots
    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 1000),
               plot_title=title_text, titlefontsize=10)
end


"""
    plot_reco_track(reco_df; markersize=4.0, markersize_extremes=10.0, linewidth=2.0)

Plot a reconstructed track (from reconstruct_central_path) showing smoothed
path in three projections and 3D view.

# Arguments
- `reco_df::DataFrame`: Reconstructed track with columns x, y, z, energy_sum
- `markersize::Float64`: Size for path point markers (default: 4.0)
- `markersize_extremes::Float64`: Size for extreme point markers (default: 10.0)
- `linewidth::Float64`: Width of path line (default: 2.0)

# Returns
- A plot with XY, XZ, YZ projections and 3D view
"""
function plot_reco_track(reco_df::DataFrame;
                         markersize::Float64=4.0,
                         markersize_extremes::Float64=10.0,
                         linewidth::Float64=2.0)

    nrow(reco_df) == 0 && error("Empty reconstructed track DataFrame")

    # Extract positions and energy
    x = Float64.(reco_df.x)
    y = Float64.(reco_df.y)
    z = Float64.(reco_df.z)
    e = Float64.(reco_df.energy_sum)

    # Compute plot limits with padding
    xmid, xrange = mean((minimum(x), maximum(x))), maximum(x) - minimum(x)
    ymid, yrange = mean((minimum(y), maximum(y))), maximum(y) - minimum(y)
    zmid, zrange = mean((minimum(z), maximum(z))), maximum(z) - minimum(z)

    xlim = (xmid - 0.6 * xrange, xmid + 0.6 * xrange)
    ylim = (ymid - 0.6 * yrange, ymid + 0.6 * yrange)
    zlim = (zmid - 0.6 * zrange, zmid + 0.6 * zrange)

    cmap = cgrad(:viridis)

    # Calculate total path length
    total_length = 0.0
    for i in 1:length(x)-1
        total_length += sqrt((x[i+1]-x[i])^2 + (y[i+1]-y[i])^2 + (z[i+1]-z[i])^2)
    end

    # === XY Projection ===
    p1 = scatter(x, y, marker_z=e, ms=markersize,
                xlabel="X (mm)", ylabel="Y (mm)", title="XY Projection",
                xlims=xlim, ylims=ylim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    # Draw path line
    plot!(p1, x, y, color=:orange, linewidth=linewidth, alpha=0.7, label="")
    # Mark extremes
    scatter!(p1, [x[1]], [y[1]], ms=markersize_extremes, color=:red,
            markershape=:circle, label="Start", markerstrokewidth=3, markerstrokecolor=:white)
    scatter!(p1, [x[end]], [y[end]], ms=markersize_extremes, color=:darkred,
            markershape=:circle, label="End", markerstrokewidth=3, markerstrokecolor=:white)

    # === XZ Projection ===
    p2 = scatter(x, z, marker_z=e, ms=markersize,
                xlabel="X (mm)", ylabel="Z (mm)", title="XZ Projection",
                xlims=xlim, ylims=zlim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    plot!(p2, x, z, color=:orange, linewidth=linewidth, alpha=0.7, label="")
    scatter!(p2, [x[1]], [z[1]], ms=markersize_extremes, color=:red,
            markershape=:circle, label="", markerstrokewidth=3, markerstrokecolor=:white)
    scatter!(p2, [x[end]], [z[end]], ms=markersize_extremes, color=:darkred,
            markershape=:circle, label="", markerstrokewidth=3, markerstrokecolor=:white)

    # === YZ Projection ===
    p3 = scatter(y, z, marker_z=e, ms=markersize,
                xlabel="Y (mm)", ylabel="Z (mm)", title="YZ Projection",
                xlims=ylim, ylims=zlim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    plot!(p3, y, z, color=:orange, linewidth=linewidth, alpha=0.7, label="")
    scatter!(p3, [y[1]], [z[1]], ms=markersize_extremes, color=:red,
            markershape=:circle, label="", markerstrokewidth=3, markerstrokecolor=:white)
    scatter!(p3, [y[end]], [z[end]], ms=markersize_extremes, color=:darkred,
            markershape=:circle, label="", markerstrokewidth=3, markerstrokecolor=:white)

    # === 3D View ===
    p4 = scatter(x, y, z, marker_z=e, ms=markersize,
                xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)",
                title="3D View",
                xlims=xlim, ylims=ylim, zlims=zlim,
                colorbar_title="Energy (MeV)", color=cmap, label="",
                markerstrokewidth=0)
    plot!(p4, x, y, z, color=:orange, linewidth=linewidth, alpha=0.7, label="")
    scatter!(p4, [x[1]], [y[1]], [z[1]], ms=markersize_extremes, color=:red,
            markershape=:circle, label="Start", markerstrokewidth=3, markerstrokecolor=:white)
    scatter!(p4, [x[end]], [y[end]], [z[end]], ms=markersize_extremes, color=:darkred,
            markershape=:circle, label="End", markerstrokewidth=3, markerstrokecolor=:white)

    # Title with info
    title_text = "Reconstructed Track"
    title_text *= "\nPath length: $(round(total_length, digits=2)) mm"
    title_text *= " | Points: $(length(x))"
    title_text *= " | Total energy: $(round(sum(e), digits=3)) MeV"

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 1000),
               plot_title=title_text, titlefontsize=10)
end


"""
    plot_reco_track_with_voxels(track, reco_df;
                                markersize_voxels=3.0, markersize_reco=5.0,
                                linewidth=2.0, alpha_voxels=0.5)

Plot a reconstructed track overlaid on the original track voxels.

# Arguments
- `track::Tracks`: Original track with all voxels
- `reco_df::DataFrame`: Reconstructed track from `reconstruct_central_path`
- `markersize_voxels::Float64`: Size for original voxel markers (default: 3.0)
- `markersize_reco::Float64`: Size for reconstructed path markers (default: 5.0)
- `linewidth::Float64`: Width of reconstructed path line (default: 2.0)
- `alpha_voxels::Float64`: Transparency for original voxels (default: 0.5)

# Returns
- A plot with XY, XZ, YZ projections and 3D view
"""
function plot_reco_track_with_voxels(track::Tracks, reco_df::DataFrame;
                                     markersize_voxels::Float64=3.0,
                                     markersize_reco::Float64=5.0,
                                     linewidth::Float64=2.0,
                                     alpha_voxels::Float64=0.5)

    nrow(reco_df) == 0 && error("Empty reconstructed track DataFrame")

    # Extract original voxel positions and energy
    vx = Float64.(track.voxels.x)
    vy = Float64.(track.voxels.y)
    vz = Float64.(track.voxels.z)
    ve = Float64.(track.voxels.energy)

    # Extract reconstructed path positions and energy
    rx = Float64.(reco_df.x)
    ry = Float64.(reco_df.y)
    rz = Float64.(reco_df.z)
    re = Float64.(reco_df.energy_sum)

    # Compute plot limits with padding (use voxels extent)
    xmid, xrange = mean((minimum(vx), maximum(vx))), maximum(vx) - minimum(vx)
    ymid, yrange = mean((minimum(vy), maximum(vy))), maximum(vy) - minimum(vy)
    zmid, zrange = mean((minimum(vz), maximum(vz))), maximum(vz) - minimum(vz)

    xlim = (xmid - 0.6 * xrange, xmid + 0.6 * xrange)
    ylim = (ymid - 0.6 * yrange, ymid + 0.6 * yrange)
    zlim = (zmid - 0.6 * zrange, zmid + 0.6 * zrange)

    cmap = cgrad(:viridis)

    # Calculate reconstructed path length
    total_length = 0.0
    for i in 1:length(rx)-1
        total_length += sqrt((rx[i+1]-rx[i])^2 + (ry[i+1]-ry[i])^2 + (rz[i+1]-rz[i])^2)
    end

    # === XY Projection ===
    p1 = scatter(vx, vy, marker_z=ve, ms=markersize_voxels, alpha=alpha_voxels,
                xlabel="X (mm)", ylabel="Y (mm)", title="XY Projection",
                xlims=xlim, ylims=ylim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    # Overlay reconstructed path
    scatter!(p1, rx, ry, ms=markersize_reco, color=:orange, markershape=:circle,
            label="", markerstrokewidth=1, markerstrokecolor=:black)
    plot!(p1, rx, ry, color=:orange, linewidth=linewidth, alpha=0.8, label="")
    # Mark extremes
    scatter!(p1, [rx[1]], [ry[1]], ms=markersize_reco+3, color=:red,
            markershape=:circle, label="Start", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p1, [rx[end]], [ry[end]], ms=markersize_reco+3, color=:darkred,
            markershape=:circle, label="End", markerstrokewidth=2, markerstrokecolor=:white)

    # === XZ Projection ===
    p2 = scatter(vx, vz, marker_z=ve, ms=markersize_voxels, alpha=alpha_voxels,
                xlabel="X (mm)", ylabel="Z (mm)", title="XZ Projection",
                xlims=xlim, ylims=zlim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    scatter!(p2, rx, rz, ms=markersize_reco, color=:orange, markershape=:circle,
            label="", markerstrokewidth=1, markerstrokecolor=:black)
    plot!(p2, rx, rz, color=:orange, linewidth=linewidth, alpha=0.8, label="")
    scatter!(p2, [rx[1]], [rz[1]], ms=markersize_reco+3, color=:red,
            markershape=:circle, label="", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p2, [rx[end]], [rz[end]], ms=markersize_reco+3, color=:darkred,
            markershape=:circle, label="", markerstrokewidth=2, markerstrokecolor=:white)

    # === YZ Projection ===
    p3 = scatter(vy, vz, marker_z=ve, ms=markersize_voxels, alpha=alpha_voxels,
                xlabel="Y (mm)", ylabel="Z (mm)", title="YZ Projection",
                xlims=ylim, ylims=zlim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    scatter!(p3, ry, rz, ms=markersize_reco, color=:orange, markershape=:circle,
            label="", markerstrokewidth=1, markerstrokecolor=:black)
    plot!(p3, ry, rz, color=:orange, linewidth=linewidth, alpha=0.8, label="")
    scatter!(p3, [ry[1]], [rz[1]], ms=markersize_reco+3, color=:red,
            markershape=:circle, label="", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p3, [ry[end]], [rz[end]], ms=markersize_reco+3, color=:darkred,
            markershape=:circle, label="", markerstrokewidth=2, markerstrokecolor=:white)

    # === 3D View ===
    p4 = scatter(vx, vy, vz, marker_z=ve, ms=markersize_voxels, alpha=alpha_voxels,
                xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)",
                title="3D View",
                xlims=xlim, ylims=ylim, zlims=zlim,
                colorbar_title="Energy (MeV)", color=cmap, label="",
                markerstrokewidth=0)
    scatter!(p4, rx, ry, rz, ms=markersize_reco, color=:orange, markershape=:circle,
            label="Reco path", markerstrokewidth=1, markerstrokecolor=:black)
    plot!(p4, rx, ry, rz, color=:orange, linewidth=linewidth, alpha=0.8, label="")
    scatter!(p4, [rx[1]], [ry[1]], [rz[1]], ms=markersize_reco+3, color=:red,
            markershape=:circle, label="Start", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p4, [rx[end]], [ry[end]], [rz[end]], ms=markersize_reco+3, color=:darkred,
            markershape=:circle, label="End", markerstrokewidth=2, markerstrokecolor=:white)

    # Title with info
    title_text = "Reconstructed Track with Voxels"
    title_text *= "\nVoxels: $(length(vx)) | Reco points: $(length(rx))"
    title_text *= " | Path length: $(round(total_length, digits=2)) mm"
    title_text *= "\nTotal energy: $(round(sum(ve)*1e3, digits=1)) keV"

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 1000),
               plot_title=title_text, titlefontsize=10)
end


"""
    plot_reco_track_with_voxels_and_spheres(track, reco_df, blob_result, sphere_radius;
                                            markersize_voxels=3.0, markersize_reco=5.0,
                                            linewidth=2.0, alpha_voxels=0.5, alpha_spheres=0.3)

Plot a reconstructed track overlaid on original voxels with blob spheres.

# Arguments
- `track::Tracks`: Original track with all voxels
- `reco_df::DataFrame`: Reconstructed track from `reconstruct_central_path`
- `blob_result`: Output from `find_blob_energies` containing blob1, blob2, Eb1, Eb2
- `sphere_radius::Float64`: Radius of blob spheres in mm
- `markersize_voxels::Float64`: Size for original voxel markers (default: 3.0)
- `markersize_reco::Float64`: Size for reconstructed path markers (default: 5.0)
- `linewidth::Float64`: Width of reconstructed path line (default: 2.0)
- `alpha_voxels::Float64`: Transparency for original voxels (default: 0.5)
- `alpha_spheres::Float64`: Transparency for blob spheres (default: 0.3)

# Returns
- A plot with XY, XZ, YZ projections and 3D view showing voxels, reco path, and blob spheres
"""
function plot_reco_track_with_voxels_and_spheres(track::Tracks, reco_df::DataFrame,
                                                  blob_result, sphere_radius::Float64;
                                                  markersize_voxels::Float64=3.0,
                                                  markersize_reco::Float64=5.0,
                                                  linewidth::Float64=2.0,
                                                  alpha_voxels::Float64=0.5,
                                                  alpha_spheres::Float64=0.3)

    nrow(reco_df) == 0 && error("Empty reconstructed track DataFrame")

    # Helper: draw circle in 2D
    function draw_circle_2d(cx, cy, r)
        θ = range(0, 2π, length=100)
        cx .+ r * cos.(θ), cy .+ r * sin.(θ)
    end

    # Helper: draw wireframe sphere in 3D
    function draw_wireframe_sphere!(p, cx, cy, cz, r, color, n_lines=8)
        θ = range(0, 2π, length=50)
        for i in 1:n_lines
            ϕ = (i-1) * π / n_lines
            plot!(p, cx .+ r*sin(ϕ).*cos.(θ), cy .+ r*sin(ϕ).*sin.(θ),
                  cz .+ r*cos(ϕ).*ones(length(θ)), color=color, lw=2, alpha=0.6, label="")
            angle = (i-1) * π / (n_lines-1) - π/2
            r_circ = r * cos(angle)
            plot!(p, cx .+ r_circ.*cos.(θ), cy .+ r_circ.*sin.(θ),
                  (cz + r*sin(angle)).*ones(length(θ)), color=color, lw=2, alpha=0.6, label="")
        end
    end

    # Extract original voxel positions and energy
    vx = Float64.(track.voxels.x)
    vy = Float64.(track.voxels.y)
    vz = Float64.(track.voxels.z)
    ve = Float64.(track.voxels.energy)

    # Extract reconstructed path positions
    rx = Float64.(reco_df.x)
    ry = Float64.(reco_df.y)
    rz = Float64.(reco_df.z)

    # Get blob info
    b1 = blob_result.blob1
    b2 = blob_result.blob2
    Eb1 = blob_result.Eb1
    Eb2 = blob_result.Eb2

    # Colors: blob1 (higher energy) = red, blob2 (lower energy) = blue
    color_b1 = :red
    color_b2 = :blue

    # Compute plot limits with padding for spheres (use voxels extent)
    xmid, xrange = mean((minimum(vx), maximum(vx))), maximum(vx) - minimum(vx)
    ymid, yrange = mean((minimum(vy), maximum(vy))), maximum(vy) - minimum(vy)
    zmid, zrange = mean((minimum(vz), maximum(vz))), maximum(vz) - minimum(vz)

    padding = 0.7 + sphere_radius / max(1.0, min(xrange, yrange, zrange))
    xlim = (xmid - padding * xrange, xmid + padding * xrange)
    ylim = (ymid - padding * yrange, ymid + padding * yrange)
    zlim = (zmid - padding * zrange, zmid + padding * zrange)

    cmap = cgrad(:viridis)

    # Calculate reconstructed path length
    total_length = 0.0
    for i in 1:length(rx)-1
        total_length += sqrt((rx[i+1]-rx[i])^2 + (ry[i+1]-ry[i])^2 + (rz[i+1]-rz[i])^2)
    end

    # === XY Projection ===
    p1 = scatter(vx, vy, marker_z=ve, ms=markersize_voxels, alpha=alpha_voxels,
                xlabel="X (mm)", ylabel="Y (mm)", title="XY Projection",
                xlims=xlim, ylims=ylim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    # Overlay reconstructed path
    scatter!(p1, rx, ry, ms=markersize_reco, color=:orange, markershape=:circle,
            label="", markerstrokewidth=1, markerstrokecolor=:black)
    plot!(p1, rx, ry, color=:orange, linewidth=linewidth, alpha=0.8, label="")
    # Draw blob spheres
    cx, cy = draw_circle_2d(b1.x, b1.y, sphere_radius)
    plot!(p1, cx, cy, color=color_b1, fill=(0, alpha_spheres, color_b1), lw=2, label="")
    cx, cy = draw_circle_2d(b2.x, b2.y, sphere_radius)
    plot!(p1, cx, cy, color=color_b2, fill=(0, alpha_spheres, color_b2), lw=2, label="")
    # Mark blob centers
    scatter!(p1, [b1.x], [b1.y], ms=8, color=color_b1, markershape=:circle,
            label="Eb1=$(round(Eb1, digits=3)) MeV", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p1, [b2.x], [b2.y], ms=8, color=color_b2, markershape=:circle,
            label="Eb2=$(round(Eb2, digits=3)) MeV", markerstrokewidth=2, markerstrokecolor=:white)

    # === XZ Projection ===
    p2 = scatter(vx, vz, marker_z=ve, ms=markersize_voxels, alpha=alpha_voxels,
                xlabel="X (mm)", ylabel="Z (mm)", title="XZ Projection",
                xlims=xlim, ylims=zlim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    scatter!(p2, rx, rz, ms=markersize_reco, color=:orange, markershape=:circle,
            label="", markerstrokewidth=1, markerstrokecolor=:black)
    plot!(p2, rx, rz, color=:orange, linewidth=linewidth, alpha=0.8, label="")
    cx, cz = draw_circle_2d(b1.x, b1.z, sphere_radius)
    plot!(p2, cx, cz, color=color_b1, fill=(0, alpha_spheres, color_b1), lw=2, label="")
    cx, cz = draw_circle_2d(b2.x, b2.z, sphere_radius)
    plot!(p2, cx, cz, color=color_b2, fill=(0, alpha_spheres, color_b2), lw=2, label="")
    scatter!(p2, [b1.x], [b1.z], ms=8, color=color_b1, markershape=:circle,
            label="", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p2, [b2.x], [b2.z], ms=8, color=color_b2, markershape=:circle,
            label="", markerstrokewidth=2, markerstrokecolor=:white)

    # === YZ Projection ===
    p3 = scatter(vy, vz, marker_z=ve, ms=markersize_voxels, alpha=alpha_voxels,
                xlabel="Y (mm)", ylabel="Z (mm)", title="YZ Projection",
                xlims=ylim, ylims=zlim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    scatter!(p3, ry, rz, ms=markersize_reco, color=:orange, markershape=:circle,
            label="", markerstrokewidth=1, markerstrokecolor=:black)
    plot!(p3, ry, rz, color=:orange, linewidth=linewidth, alpha=0.8, label="")
    cy, cz = draw_circle_2d(b1.y, b1.z, sphere_radius)
    plot!(p3, cy, cz, color=color_b1, fill=(0, alpha_spheres, color_b1), lw=2, label="")
    cy, cz = draw_circle_2d(b2.y, b2.z, sphere_radius)
    plot!(p3, cy, cz, color=color_b2, fill=(0, alpha_spheres, color_b2), lw=2, label="")
    scatter!(p3, [b1.y], [b1.z], ms=8, color=color_b1, markershape=:circle,
            label="", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p3, [b2.y], [b2.z], ms=8, color=color_b2, markershape=:circle,
            label="", markerstrokewidth=2, markerstrokecolor=:white)

    # === 3D View ===
    p4 = scatter(vx, vy, vz, marker_z=ve, ms=markersize_voxels, alpha=alpha_voxels,
                xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)",
                title="3D View",
                xlims=xlim, ylims=ylim, zlims=zlim,
                colorbar=false, color=cmap, label="",
                markerstrokewidth=0)
    scatter!(p4, rx, ry, rz, ms=markersize_reco, color=:orange, markershape=:circle,
            label="Reco path", markerstrokewidth=1, markerstrokecolor=:black)
    plot!(p4, rx, ry, rz, color=:orange, linewidth=linewidth, alpha=0.8, label="")
    # Draw wireframe spheres
    draw_wireframe_sphere!(p4, b1.x, b1.y, b1.z, sphere_radius, color_b1)
    draw_wireframe_sphere!(p4, b2.x, b2.y, b2.z, sphere_radius, color_b2)
    scatter!(p4, [b1.x], [b1.y], [b1.z], ms=8, color=color_b1, markershape=:circle,
            label="Eb1=$(round(Eb1, digits=3)) MeV", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p4, [b2.x], [b2.y], [b2.z], ms=8, color=color_b2, markershape=:circle,
            label="Eb2=$(round(Eb2, digits=3)) MeV", markerstrokewidth=2, markerstrokecolor=:white)

    # Title with info
    asymmetry = blob_result.asymmetry
    title_text = "Rb=$(sphere_radius) mm | Eb1=$(round(Eb1, digits=3)) MeV | Eb2=$(round(Eb2, digits=3)) MeV"

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 1000),
               plot_title=title_text, titlefontsize=10)
end


"""
    plot_reco_track_with_blobs(reco_df, blob_result, sphere_radius;
                               markersize=4.0, linewidth=2.0, alpha_spheres=0.3)

Plot a reconstructed track with blob spheres showing Eb1 and Eb2 energies.

# Arguments
- `reco_df::DataFrame`: Reconstructed track from `reconstruct_central_path`
- `blob_result`: Output from `reconstruct_end_point_energy` containing blob1, blob2, Eb1, Eb2
- `sphere_radius::Float64`: Radius of blob spheres in mm
- `markersize::Float64`: Size for path point markers (default: 4.0)
- `linewidth::Float64`: Width of path line (default: 2.0)
- `alpha_spheres::Float64`: Transparency for sphere circles (default: 0.3)

# Returns
- A plot with XY, XZ, YZ projections and 3D view showing blobs
"""
function plot_reco_track_with_blobs(reco_df::DataFrame, blob_result, sphere_radius::Float64;
                                    markersize::Float64=4.0,
                                    linewidth::Float64=2.0,
                                    alpha_spheres::Float64=0.3)

    nrow(reco_df) == 0 && error("Empty reconstructed track DataFrame")

    # Helper: draw circle in 2D
    function draw_circle_2d(cx, cy, r)
        θ = range(0, 2π, length=100)
        cx .+ r * cos.(θ), cy .+ r * sin.(θ)
    end

    # Helper: draw wireframe sphere in 3D
    function draw_wireframe_sphere!(p, cx, cy, cz, r, color, n_lines=8)
        θ = range(0, 2π, length=50)
        for i in 1:n_lines
            ϕ = (i-1) * π / n_lines
            plot!(p, cx .+ r*sin(ϕ).*cos.(θ), cy .+ r*sin(ϕ).*sin.(θ),
                  cz .+ r*cos(ϕ).*ones(length(θ)), color=color, lw=2, alpha=0.6, label="")
            angle = (i-1) * π / (n_lines-1) - π/2
            r_circ = r * cos(angle)
            plot!(p, cx .+ r_circ.*cos.(θ), cy .+ r_circ.*sin.(θ),
                  (cz + r*sin(angle)).*ones(length(θ)), color=color, lw=2, alpha=0.6, label="")
        end
    end

    # Extract positions and energy
    x = Float64.(reco_df.x)
    y = Float64.(reco_df.y)
    z = Float64.(reco_df.z)
    e = Float64.(reco_df.energy_sum)

    # Get blob info
    b1 = blob_result.blob1
    b2 = blob_result.blob2
    Eb1 = blob_result.Eb1
    Eb2 = blob_result.Eb2

    # Colors: blob1 (higher energy) = red, blob2 (lower energy) = blue
    color_b1 = :red
    color_b2 = :blue

    # Compute plot limits with padding for spheres
    xmid, xrange = mean((minimum(x), maximum(x))), maximum(x) - minimum(x)
    ymid, yrange = mean((minimum(y), maximum(y))), maximum(y) - minimum(y)
    zmid, zrange = mean((minimum(z), maximum(z))), maximum(z) - minimum(z)

    padding = 0.7 + sphere_radius / max(1.0, min(xrange, yrange, zrange))
    xlim = (xmid - padding * xrange, xmid + padding * xrange)
    ylim = (ymid - padding * yrange, ymid + padding * yrange)
    zlim = (zmid - padding * zrange, zmid + padding * zrange)

    cmap = cgrad(:viridis)

    # Calculate total path length
    total_length = 0.0
    for i in 1:length(x)-1
        total_length += sqrt((x[i+1]-x[i])^2 + (y[i+1]-y[i])^2 + (z[i+1]-z[i])^2)
    end

    # === XY Projection ===
    p1 = scatter(x, y, marker_z=e, ms=markersize,
                xlabel="X (mm)", ylabel="Y (mm)", title="XY Projection",
                xlims=xlim, ylims=ylim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    plot!(p1, x, y, color=:orange, linewidth=linewidth, alpha=0.7, label="")
    # Draw blob spheres
    cx, cy = draw_circle_2d(b1.x, b1.y, sphere_radius)
    plot!(p1, cx, cy, color=color_b1, fill=(0, alpha_spheres, color_b1), lw=2, label="")
    cx, cy = draw_circle_2d(b2.x, b2.y, sphere_radius)
    plot!(p1, cx, cy, color=color_b2, fill=(0, alpha_spheres, color_b2), lw=2, label="")
    # Mark centers
    scatter!(p1, [b1.x], [b1.y], ms=8, color=color_b1, markershape=:circle,
            label="Eb1=$(round(Eb1*1e3, digits=1)) keV", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p1, [b2.x], [b2.y], ms=8, color=color_b2, markershape=:circle,
            label="Eb2=$(round(Eb2*1e3, digits=1)) keV", markerstrokewidth=2, markerstrokecolor=:white)

    # === XZ Projection ===
    p2 = scatter(x, z, marker_z=e, ms=markersize,
                xlabel="X (mm)", ylabel="Z (mm)", title="XZ Projection",
                xlims=xlim, ylims=zlim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    plot!(p2, x, z, color=:orange, linewidth=linewidth, alpha=0.7, label="")
    cx, cz = draw_circle_2d(b1.x, b1.z, sphere_radius)
    plot!(p2, cx, cz, color=color_b1, fill=(0, alpha_spheres, color_b1), lw=2, label="")
    cx, cz = draw_circle_2d(b2.x, b2.z, sphere_radius)
    plot!(p2, cx, cz, color=color_b2, fill=(0, alpha_spheres, color_b2), lw=2, label="")
    scatter!(p2, [b1.x], [b1.z], ms=8, color=color_b1, markershape=:circle,
            label="", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p2, [b2.x], [b2.z], ms=8, color=color_b2, markershape=:circle,
            label="", markerstrokewidth=2, markerstrokecolor=:white)

    # === YZ Projection ===
    p3 = scatter(y, z, marker_z=e, ms=markersize,
                xlabel="Y (mm)", ylabel="Z (mm)", title="YZ Projection",
                xlims=ylim, ylims=zlim, color=cmap,
                colorbar_title="Energy (MeV)", legend=false, label="",
                markerstrokewidth=0)
    plot!(p3, y, z, color=:orange, linewidth=linewidth, alpha=0.7, label="")
    cy, cz = draw_circle_2d(b1.y, b1.z, sphere_radius)
    plot!(p3, cy, cz, color=color_b1, fill=(0, alpha_spheres, color_b1), lw=2, label="")
    cy, cz = draw_circle_2d(b2.y, b2.z, sphere_radius)
    plot!(p3, cy, cz, color=color_b2, fill=(0, alpha_spheres, color_b2), lw=2, label="")
    scatter!(p3, [b1.y], [b1.z], ms=8, color=color_b1, markershape=:circle,
            label="", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p3, [b2.y], [b2.z], ms=8, color=color_b2, markershape=:circle,
            label="", markerstrokewidth=2, markerstrokecolor=:white)

    # === 3D View ===
    p4 = scatter(x, y, z, marker_z=e, ms=markersize,
                xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)",
                title="3D View",
                xlims=xlim, ylims=ylim, zlims=zlim,
                colorbar_title="Energy (MeV)", color=cmap, label="",
                markerstrokewidth=0)
    plot!(p4, x, y, z, color=:orange, linewidth=linewidth, alpha=0.7, label="")
    # Draw wireframe spheres
    draw_wireframe_sphere!(p4, b1.x, b1.y, b1.z, sphere_radius, color_b1)
    draw_wireframe_sphere!(p4, b2.x, b2.y, b2.z, sphere_radius, color_b2)
    scatter!(p4, [b1.x], [b1.y], [b1.z], ms=8, color=color_b1, markershape=:circle,
            label="Eb1=$(round(Eb1*1e3, digits=1)) keV", markerstrokewidth=2, markerstrokecolor=:white)
    scatter!(p4, [b2.x], [b2.y], [b2.z], ms=8, color=color_b2, markershape=:circle,
            label="Eb2=$(round(Eb2*1e3, digits=1)) keV", markerstrokewidth=2, markerstrokecolor=:white)

    # Title with info
    asymmetry = (Eb1 - Eb2) / (Eb1 + Eb2)
    title_text = "Reconstructed Track with Blobs (r=$(sphere_radius) mm)"
    title_text *= "\nEb1=$(round(Eb1*1e3, digits=1)) keV (red)"
    title_text *= " | Eb2=$(round(Eb2*1e3, digits=1)) keV (blue)"
    title_text *= " | Asymmetry=$(round(asymmetry, digits=3))"
    title_text *= "\nPath length: $(round(total_length, digits=2)) mm"
    title_text *= " | Total energy: $(round(sum(e)*1e3, digits=1)) keV"

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 1000),
               plot_title=title_text, titlefontsize=10)
end


function plot_track_blobs(track::Tracks, sphere_radius::Float64;
                         markersize_voxels::Float64=3.0,
                         show_connections::Bool=true,
                         alpha_connections::Float64=0.2,
                         alpha_spheres::Float64=0.3,
                         sphere_resolution::Int=20)
    """
    Plot a track showing voxels with energy-colored spheres around the extremes.

    Parameters:
    - track: A Tracks object to visualize
    - sphere_radius: Radius of spheres around extremes in mm
    - markersize_voxels: Size for voxel markers
    - show_connections: Whether to show graph edges
    - alpha_connections: Transparency for connection lines
    - alpha_spheres: Transparency for sphere surfaces
    - sphere_resolution: Resolution for sphere plotting (higher = smoother)

    Returns:
    - A plot with XY, XZ, YZ projections and 3D view with energy-colored spheres
    """

    # Get track extremes and energy in spheres
    walk_result = walk_track_from_extremes(track)
    sphere_energies = energy_in_spheres_around_extremes(track, walk_result, sphere_radius)

    # Extract voxel positions and energies
    x = Float64.(track.voxels.x)
    y = Float64.(track.voxels.y)
    z = Float64.(track.voxels.z)
    e = Float64.(track.voxels.energy)

    # Compute plot limits with padding
    xmid, xrange = mean((minimum(x), maximum(x))), maximum(x) - minimum(x)
    ymid, yrange = mean((minimum(y), maximum(y))), maximum(y) - minimum(y)
    zmid, zrange = mean((minimum(z), maximum(z))), maximum(z) - minimum(z)

    # Add extra padding to accommodate spheres
    padding_factor = 0.7 + sphere_radius / min(xrange, yrange, zrange)
    xlim = (xmid - padding_factor * xrange, xmid + padding_factor * xrange)
    ylim = (ymid - padding_factor * yrange, ymid + padding_factor * yrange)
    zlim = (zmid - padding_factor * zrange, zmid + padding_factor * zrange)

    # Color map for voxel energy
    cmap = cgrad(:viridis)

    # Create function to draw circle/sphere in 2D projections
    function draw_circle_2d(center_x, center_y, radius, n_points=100)
        θ = range(0, 2π, length=n_points)
        circle_x = center_x .+ radius * cos.(θ)
        circle_y = center_y .+ radius * sin.(θ)
        return circle_x, circle_y
    end

    # Normalize sphere energies for color mapping using blob1/blob2 format
    max_sphere_energy = max(sphere_energies.blob1_energy, sphere_energies.blob2_energy)
    if max_sphere_energy > 0
        blob1_color_intensity = sphere_energies.blob1_energy / max_sphere_energy
        blob2_color_intensity = sphere_energies.blob2_energy / max_sphere_energy
    else
        blob1_color_intensity = 0.0
        blob2_color_intensity = 0.0
    end

    # Choose colors based on energy (using a turbo colormap for better visibility)
    heat_cmap = cgrad(:turbo)
    blob1_color = heat_cmap[blob1_color_intensity]  # Higher energy gets brighter color
    blob2_color = heat_cmap[blob2_color_intensity]  # Lower energy gets darker but still visible color

    # Map blob colors to start/end positions based on which blob corresponds to which extreme
    if sphere_energies.blob1_center == (walk_result.extremes[1].x, walk_result.extremes[1].y, walk_result.extremes[1].z)
        start_sphere_color = blob1_color
        end_sphere_color = blob2_color
    else
        start_sphere_color = blob2_color
        end_sphere_color = blob1_color
    end

    # === XY Projection ===
    p1 = scatter(x, y, marker_z=e, ms=markersize_voxels,
                xlabel="X (mm)", ylabel="Y (mm)", title="XY Projection",
                xlims=xlim, ylims=ylim, color=cmap,
                colorbar_title="Voxel Energy", legend=false, label="",
                markerstrokewidth=0)

    # Add connections if requested
    if show_connections && nv(track.graph) > 0
        for edge in edges(track.graph)
            i, j = src(edge), dst(edge)
            plot!(p1, [x[i], x[j]], [y[i], y[j]],
                 color=:gray, alpha=alpha_connections, label="")
        end
    end

    # Add sphere projections if we have valid extremes
    if !isnothing(walk_result.extremes[1])
        start_voxel, end_voxel = walk_result.extremes

        # Draw start sphere
        circle_x, circle_y = draw_circle_2d(start_voxel.x, start_voxel.y, sphere_radius)
        plot!(p1, circle_x, circle_y, color=start_sphere_color,
             fill=(0, alpha_spheres, start_sphere_color), linewidth=2, label="")

        # Draw end sphere
        circle_x, circle_y = draw_circle_2d(end_voxel.x, end_voxel.y, sphere_radius)
        plot!(p1, circle_x, circle_y, color=end_sphere_color,
             fill=(0, alpha_spheres, end_sphere_color), linewidth=2, label="")

        # Mark centers
        scatter!(p1, [start_voxel.x], [start_voxel.y],
                ms=5, color=:white, markershape=:circle,
                markerstrokewidth=2, markerstrokecolor=start_sphere_color, label="")
        scatter!(p1, [end_voxel.x], [end_voxel.y],
                ms=5, color=:white, markershape=:circle,
                markerstrokewidth=2, markerstrokecolor=end_sphere_color, label="")
    end

    # === XZ Projection ===
    p2 = scatter(x, z, marker_z=e, ms=markersize_voxels,
                xlabel="X (mm)", ylabel="Z (mm)", title="XZ Projection",
                xlims=xlim, ylims=zlim, color=cmap,
                colorbar_title="Voxel Energy", legend=false, label="",
                markerstrokewidth=0)

    if show_connections && nv(track.graph) > 0
        for edge in edges(track.graph)
            i, j = src(edge), dst(edge)
            plot!(p2, [x[i], x[j]], [z[i], z[j]],
                 color=:gray, alpha=alpha_connections, label="")
        end
    end

    if !isnothing(walk_result.extremes[1])
        start_voxel, end_voxel = walk_result.extremes

        # Draw start sphere
        circle_x, circle_z = draw_circle_2d(start_voxel.x, start_voxel.z, sphere_radius)
        plot!(p2, circle_x, circle_z, color=start_sphere_color,
             fill=(0, alpha_spheres, start_sphere_color), linewidth=2, label="")

        # Draw end sphere
        circle_x, circle_z = draw_circle_2d(end_voxel.x, end_voxel.z, sphere_radius)
        plot!(p2, circle_x, circle_z, color=end_sphere_color,
             fill=(0, alpha_spheres, end_sphere_color), linewidth=2, label="")

        # Mark centers
        scatter!(p2, [start_voxel.x], [start_voxel.z],
                ms=5, color=:white, markershape=:circle,
                markerstrokewidth=2, markerstrokecolor=start_sphere_color, label="")
        scatter!(p2, [end_voxel.x], [end_voxel.z],
                ms=5, color=:white, markershape=:circle,
                markerstrokewidth=2, markerstrokecolor=end_sphere_color, label="")
    end

    # === YZ Projection ===
    p3 = scatter(y, z, marker_z=e, ms=markersize_voxels,
                xlabel="Y (mm)", ylabel="Z (mm)", title="YZ Projection",
                xlims=ylim, ylims=zlim, color=cmap,
                colorbar_title="Voxel Energy", legend=false, label="",
                markerstrokewidth=0)

    if show_connections && nv(track.graph) > 0
        for edge in edges(track.graph)
            i, j = src(edge), dst(edge)
            plot!(p3, [y[i], y[j]], [z[i], z[j]],
                 color=:gray, alpha=alpha_connections, label="")
        end
    end

    if !isnothing(walk_result.extremes[1])
        start_voxel, end_voxel = walk_result.extremes

        # Draw start sphere
        circle_y, circle_z = draw_circle_2d(start_voxel.y, start_voxel.z, sphere_radius)
        plot!(p3, circle_y, circle_z, color=start_sphere_color,
             fill=(0, alpha_spheres, start_sphere_color), linewidth=2, label="")

        # Draw end sphere
        circle_y, circle_z = draw_circle_2d(end_voxel.y, end_voxel.z, sphere_radius)
        plot!(p3, circle_y, circle_z, color=end_sphere_color,
             fill=(0, alpha_spheres, end_sphere_color), linewidth=2, label="")

        # Mark centers
        scatter!(p3, [start_voxel.y], [start_voxel.z],
                ms=5, color=:white, markershape=:circle,
                markerstrokewidth=2, markerstrokecolor=start_sphere_color, label="")
        scatter!(p3, [end_voxel.y], [end_voxel.z],
                ms=5, color=:white, markershape=:circle,
                markerstrokewidth=2, markerstrokecolor=end_sphere_color, label="")
    end

    # === 3D View ===
    p4 = scatter(x, y, z, marker_z=e, ms=markersize_voxels,
                xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)",
                title="3D View",
                xlims=xlim, ylims=ylim, zlims=zlim,
                colorbar_title="Voxel Energy", color=cmap, label="",
                markerstrokewidth=0)

    # Add 3D connections
    if show_connections && nv(track.graph) > 0
        for edge in edges(track.graph)
            i, j = src(edge), dst(edge)
            plot!(p4, [x[i], x[j]], [y[i], y[j]], [z[i], z[j]],
                 color=:gray, alpha=alpha_connections, label="")
        end
    end

    # Add 3D wireframe spheres
    if !isnothing(walk_result.extremes[1])
        start_voxel, end_voxel = walk_result.extremes

        # Function to draw wireframe sphere with latitude and longitude circles
        function draw_wireframe_sphere!(p, center_x, center_y, center_z, radius, color, n_lines=10)
            θ = range(0, 2π, length=50)

            # Draw longitude circles (vertical circles through poles)
            for i in 1:n_lines
                ϕ = (i-1) * π / n_lines
                x_circle = center_x .+ radius * sin(ϕ) .* cos.(θ)
                y_circle = center_y .+ radius * sin(ϕ) .* sin.(θ)
                z_circle = center_z .+ radius * cos(ϕ) .* ones(length(θ))
                plot!(p, x_circle, y_circle, z_circle,
                     color=color, linewidth=1.5, alpha=0.6, label="")
            end

            # Draw latitude circles (horizontal circles at different heights)
            for i in 1:n_lines
                angle = (i-1) * π / (n_lines-1) - π/2
                z_level = center_z + radius * sin(angle)
                r_circle = radius * cos(angle)
                x_circle = center_x .+ r_circle .* cos.(θ)
                y_circle = center_y .+ r_circle .* sin.(θ)
                z_circle = z_level .* ones(length(θ))
                plot!(p, x_circle, y_circle, z_circle,
                     color=color, linewidth=1.5, alpha=0.6, label="")
            end
        end

        # Draw start sphere wireframe
        draw_wireframe_sphere!(p4, start_voxel.x, start_voxel.y, start_voxel.z,
                              sphere_radius, start_sphere_color, sphere_resolution÷2)

        # Draw end sphere wireframe
        draw_wireframe_sphere!(p4, end_voxel.x, end_voxel.y, end_voxel.z,
                              sphere_radius, end_sphere_color, sphere_resolution÷2)

        # Mark centers
        scatter!(p4, [start_voxel.x], [start_voxel.y], [start_voxel.z],
                ms=5, color=:white, markershape=:circle,
                markerstrokewidth=2, markerstrokecolor=start_sphere_color, label="")
        scatter!(p4, [end_voxel.x], [end_voxel.y], [end_voxel.z],
                ms=5, color=:white, markershape=:circle,
                markerstrokewidth=2, markerstrokecolor=end_sphere_color, label="")
    end

    # Add info text
    title_text = "Track with Energy Blobs (radius=$sphere_radius mm)"
    if !isnothing(walk_result.extremes[1])
        title_text *= "\nStart sphere: $(round(sphere_energies.start_sphere_energy, digits=3)) MeV"
        title_text *= " ($(sphere_energies.start_voxel_count) voxels)"
        title_text *= " | End sphere: $(round(sphere_energies.end_sphere_energy, digits=3)) MeV"
        title_text *= " ($(sphere_energies.end_voxel_count) voxels)"
    end

    # Create a custom colorbar for sphere energies
    energy_colorbar = scatter([0], [0], zcolor=[0, max_sphere_energy],
                            clims=(0, max_sphere_energy), color=heat_cmap,
                            colorbar_title="Sphere Energy (MeV)",
                            markersize=0, showaxis=false, grid=false,
                            xlims=(1,0), ylims=(1,0), label="")

    # Create two separate plots: one for 2D projections, one for 3D view

    # 2D projections plot with colorbar
    plot_2d = plot(p1, p2, p3, energy_colorbar,
                   layout=(2, 2), size=(1200, 1000),
                   plot_title="2D Projections - $title_text",
                   titlefontsize=10,
                   left_margin=8Plots.mm, bottom_margin=5Plots.mm)

    # 3D plot (large and standalone)
    plot_3d = plot(p4, size=(1000, 1000),
                   plot_title="3D View - $title_text",
                   titlefontsize=10,
                   left_margin=8Plots.mm, bottom_margin=5Plots.mm)

    return plot_2d, plot_3d
end


"""
    plot_event(evtdf; nbins=100, energy_column=:energy, reco_path=nothing, mc_path=nothing)

Plot a single event DataFrame with XY, XZ, YZ heatmaps and 3D scatter.
Optionally overlay reconstructed path and/or MC path.

# Arguments
- `evtdf`: DataFrame containing event data with columns x, y, z, energy
- `nbins`: Number of bins for histograms (default: 100)
- `energy_column`: Column name for energy values (default: :energy)
- `reco_path`: Optional DataFrame with reco path (x, y, z, s columns)
- `mc_path`: Optional DataFrame with MC path (x, y, z, energy, s columns)

# Returns
- A plot with 4 panels: XY, XZ, YZ heatmaps and 3D scatter
"""
function plot_event(evtdf::DataFrame; nbins::Int=100, energy_column::Symbol=:energy,
                    reco_path::Union{DataFrame, Nothing}=nothing,
                    mc_path::Union{DataFrame, Nothing}=nothing)
    plot_hits(evtdf; nbins=nbins, energy_column=energy_column,
              reco_path=reco_path, mc_path=mc_path)
end


function plot_hits_evt(hitsdf::DataFrame, index::Int; nbins=100, energy_column::Symbol=:energy)
    """
    Plot a single event from a hits DataFrame.

    Parameters:
    - hitsdf: DataFrame containing hits from multiple events
    - index: Event ID to plot
    - nbins: Number of bins for histograms (default: 100)
    - energy_column: Column name for energy values (default: :energy)

    Returns:
    - A plot with 4 panels: XY, XZ, YZ heatmaps and 3D scatter
    """
    eventdf = get_event(hitsdf, index)
    plot_event(eventdf; nbins=nbins, energy_column=energy_column)
end


function plot_hits(df::DataFrame; nbins::Int=100, energy_column::Symbol=:energy,
                   reco_path::Union{DataFrame, Nothing}=nothing,
                   mc_path::Union{DataFrame, Nothing}=nothing)

    x = Float64.(df.x)
    y = Float64.(df.y)
    z = Float64.(df.z)
    e = Float64.(df[!, energy_column])

    # Compute padded limits (1.3x range)
    xmid, xrange = mean((minimum(x), maximum(x))), maximum(x) - minimum(x)
    ymid, yrange = mean((minimum(y), maximum(y))), maximum(y) - minimum(y)
    zmid, zrange = mean((minimum(z), maximum(z))), maximum(z) - minimum(z)

    xlim = (xmid - 0.65 * xrange, xmid + 0.65 * xrange)
    ylim = (ymid - 0.65 * yrange, ymid + 0.65 * yrange)
    zlim = (zmid - 0.65 * zrange, zmid + 0.65 * zrange)

    cmap = cgrad(:viridis, alpha=1.0)

    # === XY ===
    h_xy = fit(Histogram, (x, y), nbins=nbins)
    wxy = h_xy.weights
    wxy_masked = map(v -> v == 0.0 ? NaN : v, wxy)
    xcenters_xy = diff(h_xy.edges[1]) ./ 2 .+ h_xy.edges[1][1:end-1]
    ycenters_xy = diff(h_xy.edges[2]) ./ 2 .+ h_xy.edges[2][1:end-1]
    p1 = heatmap(xcenters_xy, ycenters_xy, wxy_masked';
        xlabel="x", ylabel="y", title="XY Heatmap",
        xlims=xlim, ylims=ylim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === XZ ===
    h_xz = fit(Histogram, (x, z), nbins=nbins)
    wxz_masked = map(v -> v == 0.0 ? NaN : v, h_xz.weights)
    xcenters_xz = diff(h_xz.edges[1]) ./ 2 .+ h_xz.edges[1][1:end-1]
    zcenters_xz = diff(h_xz.edges[2]) ./ 2 .+ h_xz.edges[2][1:end-1]
    p2 = heatmap(xcenters_xz, zcenters_xz, wxz_masked';
        xlabel="x", ylabel="z", title="XZ Heatmap",
        xlims=xlim, ylims=zlim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === YZ ===
    h_yz = fit(Histogram, (y, z), nbins=nbins)
    wyz_masked = map(v -> v == 0.0 ? NaN : v, h_yz.weights)
    ycenters_yz = diff(h_yz.edges[1]) ./ 2 .+ h_yz.edges[1][1:end-1]
    zcenters_yz = diff(h_yz.edges[2]) ./ 2 .+ h_yz.edges[2][1:end-1]
    p3 = heatmap(ycenters_yz, zcenters_yz, wyz_masked';
        xlabel="y", ylabel="z", title="YZ Heatmap",
        xlims=ylim, ylims=zlim, cgrad=cmap, nan_color=:white, colorbar_title="Counts")

    # === 3D scatter ===
    p4 = scatter(x, y, z, marker_z=e, ms=2,
        xlabel="x", ylabel="y", zlabel="z", title="3D Scatter",
        xlims=xlim, ylims=ylim, zlims=zlim,
        colorbar_title="Energy", legend=false, cgrad=cmap)

    # Overlay MC path if provided
    if !isnothing(mc_path) && nrow(mc_path) > 0
        xm, ym, zm = Float64.(mc_path.x), Float64.(mc_path.y), Float64.(mc_path.z)

        # Check for double-beta (two electrons)
        has_two_electrons = hasproperty(mc_path, :primary_electron) &&
                            (1 in mc_path.primary_electron) && (2 in mc_path.primary_electron)

        if has_two_electrons
            idx1 = findall(mc_path.primary_electron .== 1)
            idx2 = findall(mc_path.primary_electron .== 2)

            # Electron 1: blue
            plot!(p1, xm[idx1], ym[idx1], color=:blue, linewidth=2, alpha=0.8, label="MC e1")
            scatter!(p1, [xm[idx1[1]]], [ym[idx1[1]]], ms=6, color=:blue, markershape=:diamond, label="")
            plot!(p2, xm[idx1], zm[idx1], color=:blue, linewidth=2, alpha=0.8, label="")
            scatter!(p2, [xm[idx1[1]]], [zm[idx1[1]]], ms=6, color=:blue, markershape=:diamond, label="")
            plot!(p3, ym[idx1], zm[idx1], color=:blue, linewidth=2, alpha=0.8, label="")
            scatter!(p3, [ym[idx1[1]]], [zm[idx1[1]]], ms=6, color=:blue, markershape=:diamond, label="")
            plot!(p4, xm[idx1], ym[idx1], zm[idx1], color=:blue, linewidth=2, alpha=0.8, label="MC e1")
            scatter!(p4, [xm[idx1[1]]], [ym[idx1[1]]], [zm[idx1[1]]], ms=6, color=:blue, markershape=:diamond, label="")

            # Electron 2: green
            plot!(p1, xm[idx2], ym[idx2], color=:green, linewidth=2, alpha=0.8, label="MC e2")
            scatter!(p1, [xm[idx2[end]]], [ym[idx2[end]]], ms=6, color=:green, markershape=:diamond, label="")
            plot!(p2, xm[idx2], zm[idx2], color=:green, linewidth=2, alpha=0.8, label="")
            scatter!(p2, [xm[idx2[end]]], [zm[idx2[end]]], ms=6, color=:green, markershape=:diamond, label="")
            plot!(p3, ym[idx2], zm[idx2], color=:green, linewidth=2, alpha=0.8, label="")
            scatter!(p3, [ym[idx2[end]]], [zm[idx2[end]]], ms=6, color=:green, markershape=:diamond, label="")
            plot!(p4, xm[idx2], ym[idx2], zm[idx2], color=:green, linewidth=2, alpha=0.8, label="MC e2")
            scatter!(p4, [xm[idx2[end]]], [ym[idx2[end]]], [zm[idx2[end]]], ms=6, color=:green, markershape=:diamond, label="")

            # Vertex marker (purple star)
            scatter!(p1, [xm[idx1[end]]], [ym[idx1[end]]], ms=8, color=:purple, markershape=:star5, label="vtx")
            scatter!(p2, [xm[idx1[end]]], [zm[idx1[end]]], ms=8, color=:purple, markershape=:star5, label="")
            scatter!(p3, [ym[idx1[end]]], [zm[idx1[end]]], ms=8, color=:purple, markershape=:star5, label="")
            scatter!(p4, [xm[idx1[end]]], [ym[idx1[end]]], [zm[idx1[end]]], ms=8, color=:purple, markershape=:star5, label="")
        else
            # Single electron: all blue
            plot!(p1, xm, ym, color=:blue, linewidth=2, alpha=0.8, label="MC")
            scatter!(p1, [xm[1]], [ym[1]], ms=6, color=:blue, markershape=:diamond, label="")
            scatter!(p1, [xm[end]], [ym[end]], ms=6, color=:darkblue, markershape=:diamond, label="")

            plot!(p2, xm, zm, color=:blue, linewidth=2, alpha=0.8, label="")
            scatter!(p2, [xm[1]], [zm[1]], ms=6, color=:blue, markershape=:diamond, label="")
            scatter!(p2, [xm[end]], [zm[end]], ms=6, color=:darkblue, markershape=:diamond, label="")

            plot!(p3, ym, zm, color=:blue, linewidth=2, alpha=0.8, label="")
            scatter!(p3, [ym[1]], [zm[1]], ms=6, color=:blue, markershape=:diamond, label="")
            scatter!(p3, [ym[end]], [zm[end]], ms=6, color=:darkblue, markershape=:diamond, label="")

            plot!(p4, xm, ym, zm, color=:blue, linewidth=2, alpha=0.8, label="MC")
            scatter!(p4, [xm[1]], [ym[1]], [zm[1]], ms=6, color=:blue, markershape=:diamond, label="")
            scatter!(p4, [xm[end]], [ym[end]], [zm[end]], ms=6, color=:darkblue, markershape=:diamond, label="")
        end
    end

    # Overlay RECO path if provided (red)
    if !isnothing(reco_path) && nrow(reco_path) > 0
        xr, yr, zr = Float64.(reco_path.x), Float64.(reco_path.y), Float64.(reco_path.z)
        plot!(p1, xr, yr, color=:red, linewidth=2, alpha=0.8, label="RECO")
        scatter!(p1, [xr[1]], [yr[1]], ms=6, color=:red, markershape=:circle, label="")
        scatter!(p1, [xr[end]], [yr[end]], ms=6, color=:darkred, markershape=:circle, label="")

        plot!(p2, xr, zr, color=:red, linewidth=2, alpha=0.8, label="")
        scatter!(p2, [xr[1]], [zr[1]], ms=6, color=:red, markershape=:circle, label="")
        scatter!(p2, [xr[end]], [zr[end]], ms=6, color=:darkred, markershape=:circle, label="")

        plot!(p3, yr, zr, color=:red, linewidth=2, alpha=0.8, label="")
        scatter!(p3, [yr[1]], [zr[1]], ms=6, color=:red, markershape=:circle, label="")
        scatter!(p3, [yr[end]], [zr[end]], ms=6, color=:darkred, markershape=:circle, label="")

        plot!(p4, xr, yr, zr, color=:red, linewidth=2, alpha=0.8, label="RECO")
        scatter!(p4, [xr[1]], [yr[1]], [zr[1]], ms=6, color=:red, markershape=:circle, label="")
        scatter!(p4, [xr[end]], [yr[end]], [zr[end]], ms=6, color=:darkred, markershape=:circle, label="")
    end

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
end


"""
    plot_track_blobs(track, walk_result, sphere_radius; kwargs...)

Plot track with energy-colored spheres around extremes. Returns (plot_2d, plot_3d).
"""
function plot_track_blobs(track::Tracks, walk_result;
                          sphere_radius::Float64=10.0,
                          markersize_voxels::Float64=3.0,
                          alpha_spheres::Float64=0.3,
                          sphere_resolution::Int=20)

    # Helper: draw circle in 2D
    function draw_circle_2d(cx, cy, r)
        θ = range(0, 2π, length=100)
        cx .+ r * cos.(θ), cy .+ r * sin.(θ)
    end

    # Helper: draw wireframe sphere in 3D
    function draw_wireframe_sphere!(p, cx, cy, cz, r, color, n_lines=10)
        θ = range(0, 2π, length=50)
        for i in 1:n_lines
            ϕ = (i-1) * π / n_lines
            plot!(p, cx .+ r*sin(ϕ).*cos.(θ), cy .+ r*sin(ϕ).*sin.(θ),
                  cz .+ r*cos(ϕ).*ones(length(θ)), color=color, lw=1.5, alpha=0.6, label="")
            angle = (i-1) * π / (n_lines-1) - π/2
            r_circ = r * cos(angle)
            plot!(p, cx .+ r_circ.*cos.(θ), cy .+ r_circ.*sin.(θ),
                  (cz + r*sin(angle)).*ones(length(θ)), color=color, lw=1.5, alpha=0.6, label="")
        end
    end

    # Helper: add sphere circles and center markers to 2D plot
    function add_spheres_2d!(p, c1, c2, start_color, end_color, alpha)
        for (c, col) in [(c1, start_color), (c2, end_color)]
            cx, cy = draw_circle_2d(c[1], c[2], sphere_radius)
            plot!(p, cx, cy, color=col, fill=(0, alpha, col), lw=2, label="")
            scatter!(p, [c[1]], [c[2]], ms=5, color=:white, markershape=:circle,
                    markerstrokewidth=2, markerstrokecolor=col, label="")
        end
    end

    # Helper: add sphere wireframes and center markers to 3D plot
    function add_spheres_3d!(p, start_v, end_v, start_color, end_color, res)
        for (v, col) in [(start_v, start_color), (end_v, end_color)]
            draw_wireframe_sphere!(p, v.x, v.y, v.z, sphere_radius, col, res÷2)
            scatter!(p, [v.x], [v.y], [v.z], ms=5, color=:white, markershape=:circle,
                    markerstrokewidth=2, markerstrokecolor=col, label="")
        end
    end

    # Compute sphere energies and colors
    sphere_energies = energy_in_spheres_around_extremes(track, walk_result, sphere_radius)
    x, y, z = Float64.(track.voxels.x), Float64.(track.voxels.y), Float64.(track.voxels.z)
    e = Float64.(track.voxels.energy)

    # Plot limits with padding
    compute_lim(v, pad) = (mid = mean(extrema(v)); rng = v[end] - v[1];
                          rng = maximum(v) - minimum(v);
                          (mid - pad*rng, mid + pad*rng))
    padding = 0.7 + sphere_radius / min(maximum(x)-minimum(x), maximum(y)-minimum(y), maximum(z)-minimum(z))
    xlim = (mean(extrema(x)) - padding*(maximum(x)-minimum(x)), mean(extrema(x)) + padding*(maximum(x)-minimum(x)))
    ylim = (mean(extrema(y)) - padding*(maximum(y)-minimum(y)), mean(extrema(y)) + padding*(maximum(y)-minimum(y)))
    zlim = (mean(extrema(z)) - padding*(maximum(z)-minimum(z)), mean(extrema(z)) + padding*(maximum(z)-minimum(z)))

    cmap, heat_cmap = cgrad(:viridis), cgrad(:turbo)
    max_e = max(sphere_energies.blob1_energy, sphere_energies.blob2_energy)
    blob1_int = max_e > 0 ? sphere_energies.blob1_energy / max_e : 0.0
    blob2_int = max_e > 0 ? sphere_energies.blob2_energy / max_e : 0.0
    blob1_color, blob2_color = heat_cmap[blob1_int], heat_cmap[blob2_int]

    # Map colors to start/end based on blob centers
    ext1 = walk_result.extremes[1]
    if sphere_energies.blob1_center == (ext1.x, ext1.y, ext1.z)
        start_color, end_color = blob1_color, blob2_color
    else
        start_color, end_color = blob2_color, blob1_color
    end

    # Common scatter options
    scatter_opts = (marker_z=e, ms=markersize_voxels, color=cmap,
                   colorbar_title="Voxel Energy", legend=false, label="", markerstrokewidth=0)

    # 2D projections: (data_x, data_y, xlabel, ylabel, title, xlims, ylims, coord_idx for extremes)
    projections = [
        (x, y, "X (mm)", "Y (mm)", "XY Projection", xlim, ylim, (v -> v.x, v -> v.y)),
        (x, z, "X (mm)", "Z (mm)", "XZ Projection", xlim, zlim, (v -> v.x, v -> v.z)),
        (y, z, "Y (mm)", "Z (mm)", "YZ Projection", ylim, zlim, (v -> v.y, v -> v.z))
    ]

    plots_2d = map(projections) do (dx, dy, xlab, ylab, title, xl, yl, coords)
        p = scatter(dx, dy; scatter_opts..., xlabel=xlab, ylabel=ylab, title=title, xlims=xl, ylims=yl)
        if !isnothing(walk_result.extremes[1])
            sv, ev = walk_result.extremes
            add_spheres_2d!(p, (coords[1](sv), coords[2](sv)), (coords[1](ev), coords[2](ev)),
                           start_color, end_color, alpha_spheres)
        end
        p
    end

    # 3D view
    p4 = scatter(x, y, z; scatter_opts..., xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)",
                title="3D View", xlims=xlim, ylims=ylim, zlims=zlim)
    if !isnothing(walk_result.extremes[1])
        add_spheres_3d!(p4, walk_result.extremes[1], walk_result.extremes[2],
                       start_color, end_color, sphere_resolution)
    end

    # Colorbar for sphere energies
    energy_colorbar = scatter([0], [0], zcolor=[0, max_e], clims=(0, max_e), color=heat_cmap,
                             colorbar_title="Sphere Energy (MeV)", markersize=0,
                             showaxis=false, grid=false, xlims=(1,0), ylims=(1,0), label="")

    plot_2d = plot(plots_2d..., energy_colorbar, layout=(2,2), size=(1200,1000),
                   titlefontsize=10, left_margin=8Plots.mm, bottom_margin=5Plots.mm)
    plot_3d = plot(p4, size=(1000,1000), titlefontsize=10,
                   left_margin=8Plots.mm, bottom_margin=5Plots.mm)

    return (plot_2d, plot_3d)
end


"""
    plot_track_walk(walk_result; kwargs...)

Plot the walk result path voxels with electrons as intensity and red circles around extremes.

# Arguments
- `walk_result`: Result from walk_track_from_extremes containing extremes, path_indices, path_voxels, etc.
- `markersize_voxels::Float64`: Size for voxel markers (default: 4.0)
- `extreme_radius::Float64`: Radius for circles around extremes in mm (default: 5.0)
- `alpha_circles::Float64`: Transparency for circles (default: 0.3)

# Returns
- `(plot_2d, plot_3d)`: Tuple of 2D projections plot and 3D view plot
"""
function plot_track_walk(walk_result;
                         markersize_voxels::Float64=4.0,
                         extreme_radius::Float64=5.0,
                         alpha_circles::Float64=0.3)

    # Helper: draw circle in 2D
    function draw_circle_2d(cx, cy, r)
        θ = range(0, 2π, length=100)
        cx .+ r * cos.(θ), cy .+ r * sin.(θ)
    end

    # Helper: draw wireframe sphere in 3D
    function draw_wireframe_sphere!(p, cx, cy, cz, r, color, n_lines=8)
        θ = range(0, 2π, length=50)
        for i in 1:n_lines
            ϕ = (i-1) * π / n_lines
            plot!(p, cx .+ r*sin(ϕ).*cos.(θ), cy .+ r*sin(ϕ).*sin.(θ),
                  cz .+ r*cos(ϕ).*ones(length(θ)), color=color, lw=2, alpha=0.7, label="")
            angle = (i-1) * π / (n_lines-1) - π/2
            r_circ = r * cos(angle)
            plot!(p, cx .+ r_circ.*cos.(θ), cy .+ r_circ.*sin.(θ),
                  (cz + r*sin(angle)).*ones(length(θ)), color=color, lw=2, alpha=0.7, label="")
        end
    end

    # Check if we have valid path_voxels
    if isnothing(walk_result.path_voxels) || nrow(walk_result.path_voxels) == 0
        error("walk_result.path_voxels is empty or nothing")
    end

    # Extract path voxel positions and electrons
    path_df = walk_result.path_voxels
    x = Float64.(path_df.x)
    y = Float64.(path_df.y)
    z = Float64.(path_df.z)

    # Use electrons column if available, otherwise use energy
    if hasproperty(path_df, :electrons)
        intensity = Float64.(path_df.electrons)
        intensity_label = "Electrons"
    else
        intensity = Float64.(path_df.energy)
        intensity_label = "Energy"
    end

    # Compute plot limits with padding
    padding = 0.7 + extreme_radius / max(1.0, min(maximum(x)-minimum(x), maximum(y)-minimum(y), maximum(z)-minimum(z)))
    xlim = (mean(extrema(x)) - padding*(maximum(x)-minimum(x)), mean(extrema(x)) + padding*(maximum(x)-minimum(x)))
    ylim = (mean(extrema(y)) - padding*(maximum(y)-minimum(y)), mean(extrema(y)) + padding*(maximum(y)-minimum(y)))
    zlim = (mean(extrema(z)) - padding*(maximum(z)-minimum(z)), mean(extrema(z)) + padding*(maximum(z)-minimum(z)))

    cmap = cgrad(:viridis)

    # Common scatter options
    scatter_opts = (marker_z=intensity, ms=markersize_voxels, color=cmap,
                   colorbar_title=intensity_label, legend=false, label="", markerstrokewidth=0)

    # Get extremes
    start_voxel = walk_result.extremes[1]
    end_voxel = walk_result.extremes[2]

    # === XY Projection ===
    p1 = scatter(x, y; scatter_opts..., xlabel="X (mm)", ylabel="Y (mm)",
                title="XY Projection", xlims=xlim, ylims=ylim)
    if !isnothing(start_voxel)
        # Draw red circles around extremes
        cx, cy = draw_circle_2d(start_voxel.x, start_voxel.y, extreme_radius)
        plot!(p1, cx, cy, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p1, [start_voxel.x], [start_voxel.y], ms=6, color=:red, markershape=:circle, label="Start")

        cx, cy = draw_circle_2d(end_voxel.x, end_voxel.y, extreme_radius)
        plot!(p1, cx, cy, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p1, [end_voxel.x], [end_voxel.y], ms=6, color=:darkred, markershape=:circle, label="End")
    end

    # === XZ Projection ===
    p2 = scatter(x, z; scatter_opts..., xlabel="X (mm)", ylabel="Z (mm)",
                title="XZ Projection", xlims=xlim, ylims=zlim)
    if !isnothing(start_voxel)
        cx, cz = draw_circle_2d(start_voxel.x, start_voxel.z, extreme_radius)
        plot!(p2, cx, cz, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p2, [start_voxel.x], [start_voxel.z], ms=6, color=:red, markershape=:circle, label="")

        cx, cz = draw_circle_2d(end_voxel.x, end_voxel.z, extreme_radius)
        plot!(p2, cx, cz, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p2, [end_voxel.x], [end_voxel.z], ms=6, color=:darkred, markershape=:circle, label="")
    end

    # === YZ Projection ===
    p3 = scatter(y, z; scatter_opts..., xlabel="Y (mm)", ylabel="Z (mm)",
                title="YZ Projection", xlims=ylim, ylims=zlim)
    if !isnothing(start_voxel)
        cy, cz = draw_circle_2d(start_voxel.y, start_voxel.z, extreme_radius)
        plot!(p3, cy, cz, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p3, [start_voxel.y], [start_voxel.z], ms=6, color=:red, markershape=:circle, label="")

        cy, cz = draw_circle_2d(end_voxel.y, end_voxel.z, extreme_radius)
        plot!(p3, cy, cz, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p3, [end_voxel.y], [end_voxel.z], ms=6, color=:darkred, markershape=:circle, label="")
    end

    # === 3D View ===
    p4 = scatter(x, y, z; scatter_opts..., xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)",
                title="3D View", xlims=xlim, ylims=ylim, zlims=zlim)
    if !isnothing(start_voxel)
        # Draw red wireframe spheres
        draw_wireframe_sphere!(p4, start_voxel.x, start_voxel.y, start_voxel.z, extreme_radius, :red)
        scatter!(p4, [start_voxel.x], [start_voxel.y], [start_voxel.z],
                ms=6, color=:red, markershape=:circle, label="Start")

        draw_wireframe_sphere!(p4, end_voxel.x, end_voxel.y, end_voxel.z, extreme_radius, :darkred)
        scatter!(p4, [end_voxel.x], [end_voxel.y], [end_voxel.z],
                ms=6, color=:darkred, markershape=:circle, label="End")
    end

    # Title with info
    title_text = "Track Walk Path"
    if !isnothing(start_voxel)
        title_text *= "\nPath length: $(round(walk_result.total_length, digits=2)) mm"
        title_text *= " | Voxels: $(nrow(path_df))"
        title_text *= " | Confidence: $(round(walk_result.confidence, digits=2))"
    end

    # Combine plots
    plot_2d = plot(p1, p2, p3, layout=(2,2), size=(1200,1000),
                   plot_title=title_text, titlefontsize=10,
                   left_margin=8Plots.mm, bottom_margin=5Plots.mm)
    plot_3d = plot(p4, size=(1000,1000), plot_title=title_text,
                   titlefontsize=10, left_margin=8Plots.mm, bottom_margin=5Plots.mm)

    return (plot_2d, plot_3d)
end

"""
    plot_paths(reco_path, mc_path; kwargs...)

Plot reconstructed path and MC path together for comparison.

# Arguments
- `reco_path::DataFrame`: Reconstructed path with x, y, z, s columns
- `mc_path::DataFrame`: MC truth path with x, y, z, energy, s columns
- `markersize_reco::Float64`: Marker size for reco path (default: 4.0)
- `markersize_mc::Float64`: Marker size for MC path (default: 3.0)
- `linewidth::Float64`: Line width for paths (default: 2.0)
- `title::String`: Plot title (default: "RECO vs MC Path")

# Returns
- A 2x2 plot with XY, XZ, YZ projections and 3D view
"""
function plot_paths(reco_path::DataFrame, mc_path::DataFrame;
                    markersize_reco::Float64=4.0,
                    markersize_mc::Float64=3.0,
                    linewidth::Float64=2.0,
                    title::String="RECO vs MC Path")

    # Extract RECO path
    xr = Float64.(reco_path.x)
    yr = Float64.(reco_path.y)
    zr = Float64.(reco_path.z)

    # Extract MC path
    xm = Float64.(mc_path.x)
    ym = Float64.(mc_path.y)
    zm = Float64.(mc_path.z)
    em = hasproperty(mc_path, :energy) ? Float64.(mc_path.energy) : ones(nrow(mc_path))

    # Check for double-beta (two electrons)
    has_two_electrons = hasproperty(mc_path, :primary_electron) &&
                        (1 in mc_path.primary_electron) && (2 in mc_path.primary_electron)

    # Split MC path by electron if double-beta
    if has_two_electrons
        idx1 = findall(mc_path.primary_electron .== 1)
        idx2 = findall(mc_path.primary_electron .== 2)
    end

    # Compute combined limits
    all_x = vcat(xr, xm)
    all_y = vcat(yr, ym)
    all_z = vcat(zr, zm)

    xmid, xrange = mean(extrema(all_x)), maximum(all_x) - minimum(all_x)
    ymid, yrange = mean(extrema(all_y)), maximum(all_y) - minimum(all_y)
    zmid, zrange = mean(extrema(all_z)), maximum(all_z) - minimum(all_z)

    xlim = (xmid - 0.6 * max(xrange, 1.0), xmid + 0.6 * max(xrange, 1.0))
    ylim = (ymid - 0.6 * max(yrange, 1.0), ymid + 0.6 * max(yrange, 1.0))
    zlim = (zmid - 0.6 * max(zrange, 1.0), zmid + 0.6 * max(zrange, 1.0))

    # Calculate path lengths
    reco_length = nrow(reco_path) > 1 ? reco_path.s[end] : 0.0
    mc_length = nrow(mc_path) > 1 ? mc_path.s[end] : 0.0

    # Helper function to plot MC path (single or double-beta)
    function plot_mc_2d!(p, x_all, y_all, em_all, label_prefix; show_labels=true)
        if has_two_electrons
            # Electron 1: blue
            plot!(p, x_all[idx1], y_all[idx1], color=:blue, linewidth=linewidth, alpha=0.7,
                  label=show_labels ? "MC e1" : "")
            scatter!(p, x_all[idx1], y_all[idx1], ms=markersize_mc, color=:blue, alpha=0.6,
                     markerstrokewidth=0, label="")
            scatter!(p, [x_all[idx1[1]]], [y_all[idx1[1]]], ms=8, color=:blue,
                     markershape=:diamond, label=show_labels ? "e1 Bragg" : "")

            # Electron 2: green
            plot!(p, x_all[idx2], y_all[idx2], color=:green, linewidth=linewidth, alpha=0.7,
                  label=show_labels ? "MC e2" : "")
            scatter!(p, x_all[idx2], y_all[idx2], ms=markersize_mc, color=:green, alpha=0.6,
                     markerstrokewidth=0, label="")
            scatter!(p, [x_all[idx2[end]]], [y_all[idx2[end]]], ms=8, color=:green,
                     markershape=:diamond, label=show_labels ? "e2 Bragg" : "")

            # Vertex (junction between electrons) - use last of e1 or first of e2
            scatter!(p, [x_all[idx1[end]]], [y_all[idx1[end]]], ms=10, color=:purple,
                     markershape=:star5, label=show_labels ? "vertex" : "")
        else
            # Single electron: all blue
            plot!(p, x_all, y_all, color=:blue, linewidth=linewidth, alpha=0.7,
                  label=show_labels ? "MC" : "")
            scatter!(p, x_all, y_all, marker_z=em_all, ms=markersize_mc, color=:blues,
                     markerstrokewidth=0, label="", colorbar=false)
            scatter!(p, [x_all[1]], [y_all[1]], ms=8, color=:blue,
                     markershape=:diamond, label=show_labels ? "MC start" : "")
            scatter!(p, [x_all[end]], [y_all[end]], ms=8, color=:darkblue,
                     markershape=:diamond, label=show_labels ? "MC end" : "")
        end
    end

    function plot_mc_3d!(p, x_all, y_all, z_all, em_all; show_labels=true)
        if has_two_electrons
            # Electron 1: blue
            plot!(p, x_all[idx1], y_all[idx1], z_all[idx1], color=:blue,
                  linewidth=linewidth, alpha=0.7, label=show_labels ? "MC e1" : "")
            scatter!(p, x_all[idx1], y_all[idx1], z_all[idx1], ms=markersize_mc,
                     color=:blue, alpha=0.6, markerstrokewidth=0, label="")
            scatter!(p, [x_all[idx1[1]]], [y_all[idx1[1]]], [z_all[idx1[1]]], ms=8,
                     color=:blue, markershape=:diamond, label="")

            # Electron 2: green
            plot!(p, x_all[idx2], y_all[idx2], z_all[idx2], color=:green,
                  linewidth=linewidth, alpha=0.7, label=show_labels ? "MC e2" : "")
            scatter!(p, x_all[idx2], y_all[idx2], z_all[idx2], ms=markersize_mc,
                     color=:green, alpha=0.6, markerstrokewidth=0, label="")
            scatter!(p, [x_all[idx2[end]]], [y_all[idx2[end]]], [z_all[idx2[end]]], ms=8,
                     color=:green, markershape=:diamond, label="")

            # Vertex
            scatter!(p, [x_all[idx1[end]]], [y_all[idx1[end]]], [z_all[idx1[end]]], ms=10,
                     color=:purple, markershape=:star5, label="")
        else
            # Single electron
            plot!(p, x_all, y_all, z_all, color=:blue, linewidth=linewidth, alpha=0.7,
                  label=show_labels ? "MC" : "")
            scatter!(p, x_all, y_all, z_all, marker_z=em_all, ms=markersize_mc,
                     color=:blues, markerstrokewidth=0, label="", colorbar=false)
            scatter!(p, [x_all[1]], [y_all[1]], [z_all[1]], ms=8, color=:blue,
                     markershape=:diamond, label="")
            scatter!(p, [x_all[end]], [y_all[end]], [z_all[end]], ms=8, color=:darkblue,
                     markershape=:diamond, label="")
        end
    end

    # === XY Projection ===
    p1 = plot(xlabel="X (mm)", ylabel="Y (mm)", title="XY Projection",
              xlims=xlim, ylims=ylim)
    plot_mc_2d!(p1, xm, ym, em, "MC"; show_labels=true)

    plot!(p1, xr, yr, color=:red, linewidth=linewidth, alpha=0.7, label="RECO")
    scatter!(p1, xr, yr, ms=markersize_reco, color=:red, alpha=0.5,
             markerstrokewidth=0, label="")
    scatter!(p1, [xr[1]], [yr[1]], ms=8, color=:red, markershape=:circle, label="RECO start")
    scatter!(p1, [xr[end]], [yr[end]], ms=8, color=:darkred, markershape=:circle, label="RECO end")

    # === XZ Projection ===
    p2 = plot(xlabel="X (mm)", ylabel="Z (mm)", title="XZ Projection",
              xlims=xlim, ylims=zlim)
    plot_mc_2d!(p2, xm, zm, em, ""; show_labels=false)

    plot!(p2, xr, zr, color=:red, linewidth=linewidth, alpha=0.7, label="")
    scatter!(p2, xr, zr, ms=markersize_reco, color=:red, alpha=0.5,
             markerstrokewidth=0, label="")
    scatter!(p2, [xr[1]], [zr[1]], ms=8, color=:red, markershape=:circle, label="")
    scatter!(p2, [xr[end]], [zr[end]], ms=8, color=:darkred, markershape=:circle, label="")

    # === YZ Projection ===
    p3 = plot(xlabel="Y (mm)", ylabel="Z (mm)", title="YZ Projection",
              xlims=ylim, ylims=zlim)
    plot_mc_2d!(p3, ym, zm, em, ""; show_labels=false)

    plot!(p3, yr, zr, color=:red, linewidth=linewidth, alpha=0.7, label="")
    scatter!(p3, yr, zr, ms=markersize_reco, color=:red, alpha=0.5,
             markerstrokewidth=0, label="")
    scatter!(p3, [yr[1]], [zr[1]], ms=8, color=:red, markershape=:circle, label="")
    scatter!(p3, [yr[end]], [zr[end]], ms=8, color=:darkred, markershape=:circle, label="")

    # === 3D View ===
    p4 = plot(xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)", title="3D View",
              xlims=xlim, ylims=ylim, zlims=zlim)
    plot_mc_3d!(p4, xm, ym, zm, em; show_labels=true)

    plot!(p4, xr, yr, zr, color=:red, linewidth=linewidth, alpha=0.7, label="RECO")
    scatter!(p4, xr, yr, zr, ms=markersize_reco, color=:red, alpha=0.5,
             markerstrokewidth=0, label="")
    scatter!(p4, [xr[1]], [yr[1]], [zr[1]], ms=8, color=:red, markershape=:circle, label="")
    scatter!(p4, [xr[end]], [yr[end]], [zr[end]], ms=8, color=:darkred, markershape=:circle, label="")

    # Title with info
    title_text = title
    title_text *= "\nRECO: $(round(reco_length, digits=1)) mm ($(nrow(reco_path)) pts)"
    title_text *= " | MC: $(round(mc_length, digits=1)) mm ($(nrow(mc_path)) pts)"
    if has_two_electrons
        title_text *= " [double-beta]"
    end

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 1000),
                plot_title=title_text, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=5Plots.mm)
end

"""
    plot_track_with_paths(track, reco_path, mc_path; kwargs...)

Plot track voxels with both RECO and MC paths overlaid.

Shows voxels colored by energy, with:
- RECO path in red (with start/end markers)
- MC path in blue/green (two colors for double-beta)
- Extreme distance annotations (d1, d2)

# Arguments
- `track`: Tracks object with voxels
- `reco_path`: DataFrame with reco path (x, y, z, s)
- `mc_path`: DataFrame with MC path (x, y, z, energy, s, primary_electron)
- `markersize_voxels`: Size for voxel markers (default: 3.0)
- `linewidth`: Line width for paths (default: 2.0)
- `show_distances`: Show d1, d2 annotations (default: true)
- `title`: Plot title (default: "Track with RECO & MC Paths")

# Returns
- A 2x2 plot with XY, XZ, YZ projections and 3D view
"""
function plot_track_with_paths(track::Tracks, reco_path::DataFrame, mc_path::DataFrame;
                               markersize_voxels::Float64=3.0,
                               linewidth::Float64=2.0,
                               show_distances::Bool=true,
                               title::String="Track with RECO & MC Paths")

    # Extract voxel positions and energies
    xv = Float64.(track.voxels.x)
    yv = Float64.(track.voxels.y)
    zv = Float64.(track.voxels.z)
    ev = Float64.(track.voxels.energy)

    # Extract RECO path
    xr = Float64.(reco_path.x)
    yr = Float64.(reco_path.y)
    zr = Float64.(reco_path.z)

    # Extract MC path
    xm = Float64.(mc_path.x)
    ym = Float64.(mc_path.y)
    zm = Float64.(mc_path.z)

    # Check for double-beta
    has_two_electrons = hasproperty(mc_path, :primary_electron) &&
                        (1 in mc_path.primary_electron) && (2 in mc_path.primary_electron)

    if has_two_electrons
        idx1 = findall(mc_path.primary_electron .== 1)
        idx2 = findall(mc_path.primary_electron .== 2)
    end

    # Compute distances
    dists = compute_extreme_distances(reco_path, mc_path)

    # Compute plot limits with padding (use all points)
    all_x = vcat(xv, xr, xm)
    all_y = vcat(yv, yr, ym)
    all_z = vcat(zv, zr, zm)

    xmid, xrange = mean(extrema(all_x)), maximum(all_x) - minimum(all_x)
    ymid, yrange = mean(extrema(all_y)), maximum(all_y) - minimum(all_y)
    zmid, zrange = mean(extrema(all_z)), maximum(all_z) - minimum(all_z)

    xlim = (xmid - 0.6 * max(xrange, 1.0), xmid + 0.6 * max(xrange, 1.0))
    ylim = (ymid - 0.6 * max(yrange, 1.0), ymid + 0.6 * max(yrange, 1.0))
    zlim = (zmid - 0.6 * max(zrange, 1.0), zmid + 0.6 * max(zrange, 1.0))

    cmap = cgrad(:viridis)

    # Helper to plot MC path (handles single/double-beta)
    function add_mc_path_2d!(p, x_all, y_all; show_labels=true)
        if has_two_electrons
            # Electron 1: blue
            plot!(p, x_all[idx1], y_all[idx1], color=:blue, linewidth=linewidth,
                  alpha=0.8, label=show_labels ? "MC e1" : "")
            scatter!(p, [x_all[idx1[1]]], [y_all[idx1[1]]], ms=7, color=:blue,
                     markershape=:diamond, label="")
            # Electron 2: green
            plot!(p, x_all[idx2], y_all[idx2], color=:green, linewidth=linewidth,
                  alpha=0.8, label=show_labels ? "MC e2" : "")
            scatter!(p, [x_all[idx2[end]]], [y_all[idx2[end]]], ms=7, color=:green,
                     markershape=:diamond, label="")
            # Vertex
            scatter!(p, [x_all[idx1[end]]], [y_all[idx1[end]]], ms=9, color=:purple,
                     markershape=:star5, label=show_labels ? "vertex" : "")
        else
            plot!(p, x_all, y_all, color=:blue, linewidth=linewidth,
                  alpha=0.8, label=show_labels ? "MC" : "")
            scatter!(p, [x_all[1]], [y_all[1]], ms=7, color=:blue,
                     markershape=:diamond, label="")
            scatter!(p, [x_all[end]], [y_all[end]], ms=7, color=:darkblue,
                     markershape=:diamond, label="")
        end
    end

    function add_mc_path_3d!(p, x_all, y_all, z_all; show_labels=true)
        if has_two_electrons
            plot!(p, x_all[idx1], y_all[idx1], z_all[idx1], color=:blue,
                  linewidth=linewidth, alpha=0.8, label=show_labels ? "MC e1" : "")
            scatter!(p, [x_all[idx1[1]]], [y_all[idx1[1]]], [z_all[idx1[1]]], ms=7,
                     color=:blue, markershape=:diamond, label="")
            plot!(p, x_all[idx2], y_all[idx2], z_all[idx2], color=:green,
                  linewidth=linewidth, alpha=0.8, label=show_labels ? "MC e2" : "")
            scatter!(p, [x_all[idx2[end]]], [y_all[idx2[end]]], [z_all[idx2[end]]], ms=7,
                     color=:green, markershape=:diamond, label="")
            scatter!(p, [x_all[idx1[end]]], [y_all[idx1[end]]], [z_all[idx1[end]]], ms=9,
                     color=:purple, markershape=:star5, label="")
        else
            plot!(p, x_all, y_all, z_all, color=:blue, linewidth=linewidth,
                  alpha=0.8, label=show_labels ? "MC" : "")
            scatter!(p, [x_all[1]], [y_all[1]], [z_all[1]], ms=7, color=:blue,
                     markershape=:diamond, label="")
            scatter!(p, [x_all[end]], [y_all[end]], [z_all[end]], ms=7, color=:darkblue,
                     markershape=:diamond, label="")
        end
    end

    # === XY Projection ===
    p1 = scatter(xv, yv, marker_z=ev, ms=markersize_voxels, alpha=0.6,
                 xlabel="X (mm)", ylabel="Y (mm)", title="XY Projection",
                 xlims=xlim, ylims=ylim, color=cmap, markerstrokewidth=0,
                 colorbar=false, label="voxels")

    add_mc_path_2d!(p1, xm, ym; show_labels=true)

    plot!(p1, xr, yr, color=:red, linewidth=linewidth, alpha=0.8, label="RECO")
    scatter!(p1, [xr[1]], [yr[1]], ms=7, color=:red, markershape=:circle, label="")
    scatter!(p1, [xr[end]], [yr[end]], ms=7, color=:darkred, markershape=:circle, label="")

    # === XZ Projection ===
    p2 = scatter(xv, zv, marker_z=ev, ms=markersize_voxels, alpha=0.6,
                 xlabel="X (mm)", ylabel="Z (mm)", title="XZ Projection",
                 xlims=xlim, ylims=zlim, color=cmap, markerstrokewidth=0,
                 colorbar=false, label="")

    add_mc_path_2d!(p2, xm, zm; show_labels=false)

    plot!(p2, xr, zr, color=:red, linewidth=linewidth, alpha=0.8, label="")
    scatter!(p2, [xr[1]], [zr[1]], ms=7, color=:red, markershape=:circle, label="")
    scatter!(p2, [xr[end]], [zr[end]], ms=7, color=:darkred, markershape=:circle, label="")

    # === YZ Projection ===
    p3 = scatter(yv, zv, marker_z=ev, ms=markersize_voxels, alpha=0.6,
                 xlabel="Y (mm)", ylabel="Z (mm)", title="YZ Projection",
                 xlims=ylim, ylims=zlim, color=cmap, markerstrokewidth=0,
                 colorbar=false, label="")

    add_mc_path_2d!(p3, ym, zm; show_labels=false)

    plot!(p3, yr, zr, color=:red, linewidth=linewidth, alpha=0.8, label="")
    scatter!(p3, [yr[1]], [zr[1]], ms=7, color=:red, markershape=:circle, label="")
    scatter!(p3, [yr[end]], [zr[end]], ms=7, color=:darkred, markershape=:circle, label="")

    # === 3D View ===
    p4 = scatter(xv, yv, zv, marker_z=ev, ms=markersize_voxels, alpha=0.6,
                 xlabel="X", ylabel="Y", zlabel="Z", title="3D View",
                 xlims=xlim, ylims=ylim, zlims=zlim, color=cmap, markerstrokewidth=0,
                 colorbar=false, label="voxels")

    add_mc_path_3d!(p4, xm, ym, zm; show_labels=true)

    plot!(p4, xr, yr, zr, color=:red, linewidth=linewidth, alpha=0.8, label="RECO")
    scatter!(p4, [xr[1]], [yr[1]], [zr[1]], ms=7, color=:red, markershape=:circle, label="")
    scatter!(p4, [xr[end]], [yr[end]], [zr[end]], ms=7, color=:darkred, markershape=:circle, label="")

    # Build title with info
    reco_length = nrow(reco_path) > 1 ? reco_path.s[end] : 0.0
    mc_length = nrow(mc_path) > 1 ? mc_path.s[end] : 0.0

    title_text = title
    title_text *= "\nRECO: $(round(reco_length, digits=1)) mm | MC: $(round(mc_length, digits=1)) mm"
    if show_distances && !isnan(dists.d1)
        title_text *= " | d1=$(round(dists.d1, digits=1)), d2=$(round(dists.d2, digits=1)) mm"
    end
    if has_two_electrons
        title_text *= " [bb]"
    end

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 1000),
                plot_title=title_text, titlefontsize=10,
                left_margin=8Plots.mm, bottom_margin=5Plots.mm)
end
