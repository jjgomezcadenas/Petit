using Statistics: mean
using Graphs: nv, edges, src, dst

function plot_hits_trk(trk::Tracks; nbins::Int=100)
	plot_hits(trk.voxels; nbins)
end

function plot_track_with_extremes(track::Tracks, walk_result=nothing;
                                 markersize_voxels::Float64=4.0,
                                 markersize_extremes::Float64=10.0,
                                 show_connections::Bool=true,
                                 alpha_connections::Float64=0.3)
    """
    Plot a track showing voxels, connections, and extremes in three projections.

    Parameters:
    - track: A Tracks object to visualize
    - walk_result: Result from walk_track_from_extremes (optional)
    - markersize_voxels: Size for voxel markers
    - markersize_extremes: Size for extreme point markers
    - show_connections: Whether to show graph edges
    - alpha_connections: Transparency for connection lines

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

    # Add connections if requested
    if show_connections && nv(track.graph) > 0
        for edge in edges(track.graph)
            i, j = src(edge), dst(edge)
            plot!(p1, [x[i], x[j]], [y[i], y[j]],
                 color=:gray, alpha=alpha_connections, label="")
        end
    end

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

    if show_connections && nv(track.graph) > 0
        for edge in edges(track.graph)
            i, j = src(edge), dst(edge)
            plot!(p2, [x[i], x[j]], [z[i], z[j]],
                 color=:gray, alpha=alpha_connections, label="")
        end
    end

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

    if show_connections && nv(track.graph) > 0
        for edge in edges(track.graph)
            i, j = src(edge), dst(edge)
            plot!(p3, [y[i], y[j]], [z[i], z[j]],
                 color=:gray, alpha=alpha_connections, label="")
        end
    end

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

    # Add 3D connections
    if show_connections && nv(track.graph) > 0
        for edge in edges(track.graph)
            i, j = src(edge), dst(edge)
            plot!(p4, [x[i], x[j]], [y[i], y[j]], [z[i], z[j]],
                 color=:gray, alpha=alpha_connections, label="")
        end
    end

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

    # Choose colors based on energy (using a heat colormap)
    heat_cmap = cgrad(:heat, rev=true)
    blob1_color = heat_cmap[blob1_color_intensity]  # Higher energy gets "hotter" color
    blob2_color = heat_cmap[blob2_color_intensity]  # Lower energy gets "cooler" color

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

    # Add 3D spheres
    if !isnothing(walk_result.extremes[1])
        start_voxel, end_voxel = walk_result.extremes

        # Generate sphere surface points
        u = range(0, 2π, length=sphere_resolution)
        v = range(0, π, length=sphere_resolution÷2)

        # Start sphere
        sphere_x = [start_voxel.x + sphere_radius * sin(vj) * cos(ui) for ui in u, vj in v]
        sphere_y = [start_voxel.y + sphere_radius * sin(vj) * sin(ui) for ui in u, vj in v]
        sphere_z = [start_voxel.z + sphere_radius * cos(vj) for ui in u, vj in v]

        surface!(p4, sphere_x, sphere_y, sphere_z,
                color=fill(start_sphere_color, sphere_resolution, sphere_resolution÷2),
                alpha=alpha_spheres, colorbar=false, label="")

        # End sphere
        sphere_x = [end_voxel.x + sphere_radius * sin(vj) * cos(ui) for ui in u, vj in v]
        sphere_y = [end_voxel.y + sphere_radius * sin(vj) * sin(ui) for ui in u, vj in v]
        sphere_z = [end_voxel.z + sphere_radius * cos(vj) for ui in u, vj in v]

        surface!(p4, sphere_x, sphere_y, sphere_z,
                color=fill(end_sphere_color, sphere_resolution, sphere_resolution÷2),
                alpha=alpha_spheres, colorbar=false, label="")

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

    # Combine all plots with colorbar
    l = @layout [
        a b
        c d
        e{0.05h}
    ]

    return plot(p1, p2, p3, p4, energy_colorbar, layout=l, size=(1200, 1100),
               plot_title=title_text, titlefontsize=10)
end


function plot_hits_evt(hitsdf::DataFrame, index::Int; nbins=100)
    eventdf = get_event(hitsdf, index)
	plot_hits(eventdf; nbins)
end


function plot_hits(df::DataFrame; nbins::Int=100)
    
    x = Float64.(df.x)
    y = Float64.(df.y)
    z = Float64.(df.z)
    e = Float64.(df.energy)

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

    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
end