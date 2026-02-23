#!/usr/bin/env julia

"""
Single-track statistics visualization script.

Reads single_track_stats.csv from a results directory and produces visualization plots.

Usage:
    julia single_track_stats.jl <results_dir> [options]

Options:
    --no-save           Display only, don't save plots (default: save plots)
    --output-dir=X      Output directory for plots (default: same as results_dir)

Output files (saved by default):
    - single_track_energy_voxels.png (Frame 1: energy and voxels)
    - single_track_geometry.png (Frame 2: z-length and x/y extent)
    - single_track_occupancy.png (Frame 3: occupancy and time slots)
    - single_track_pixel_signals.png (Frame 4: pixel signal statistics)
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using CSV
using DataFrames
using Statistics
using Printf
using Plots

# =============================================================================
# Plotting Functions
# =============================================================================

"""
    add_stats_legend!(p, data; digits=1)

Add μ and σ as legend entries using invisible scatter plots.
"""
function add_stats_legend!(p, data; digits=1)
    μ = mean(data)
    σ = std(data)
    μ_str = round(μ, digits=digits)
    σ_str = round(σ, digits=digits)
    plot!(p, [], [], label="μ = $μ_str", seriestype=:scatter, markersize=0)
    plot!(p, [], [], label="σ = $σ_str", seriestype=:scatter, markersize=0)
    return p
end

"""
    plot_frame1_energy_voxels(df::DataFrame)

Frame 1 (1x2): Energy and Voxels histograms.
"""
function plot_frame1_energy_voxels(df::DataFrame)
    n_events = nrow(df)

    # Total energy histogram
    p1 = histogram(df.total_energy_keV, bins=50,
                   xlabel="Energy (keV)", ylabel="Events",
                   title="Total Energy", titlefontsize=10,
                   label="", color=:gray, bottom_margin=5Plots.mm)
    add_stats_legend!(p1, df.total_energy_keV; digits=1)

    # Number of voxels histogram
    p2 = histogram(df.n_voxels, bins=50,
                   xlabel="N Voxels", ylabel="Events",
                   title="Number of Voxels", titlefontsize=10,
                   label="", color=:gray, bottom_margin=5Plots.mm)
    add_stats_legend!(p2, df.n_voxels; digits=1)

    plt = plot(p1, p2, layout=(1, 2), size=(900, 400),
               plot_title="Single-Track Events (N=$n_events)", plot_titlefontsize=11)

    return plt
end

"""
    plot_frame2_geometry(df::DataFrame)

Frame 2 (1x2): Z-length histogram and X/Y extent 2D histogram.
"""
function plot_frame2_geometry(df::DataFrame)
    n_events = nrow(df)

    # Z-length histogram (30 bins, 0-120 mm)
    p1 = histogram(df.z_length_mm, bins=range(0, 120, length=31),
                   xlabel="Z Length (mm)", ylabel="Events",
                   title="Track Z-Extension", titlefontsize=10,
                   label="", color=:gray, bottom_margin=5Plots.mm,
                   xlims=(0, 120))
    add_stats_legend!(p1, df.z_length_mm; digits=1)

    # X extent vs Y extent 2D histogram (fixed 0-160 mm on both axes)
    edges = range(0, 160, length=51)
    μx, σx = round(mean(df.x_extent_mm), digits=1), round(std(df.x_extent_mm), digits=1)
    μy, σy = round(mean(df.y_extent_mm), digits=1), round(std(df.y_extent_mm), digits=1)
    p2 = histogram2d(df.x_extent_mm, df.y_extent_mm, bins=(edges, edges),
                     xlabel="X Extent (mm)", ylabel="Y Extent (mm)",
                     title="Track Bounding Box\n(μx=$μx, σx=$σx | μy=$μy, σy=$σy)", titlefontsize=9,
                     c=:viridis, colorbar=true, aspect_ratio=:equal,
                     xlims=(0, 160), ylims=(0, 160), bottom_margin=5Plots.mm)

    plt = plot(p1, p2, layout=(1, 2), size=(900, 400),
               plot_title="Track Geometry (N=$n_events)", plot_titlefontsize=11)

    return plt
end

"""
    plot_frame3_occupancy(df::DataFrame)

Frame 3 (1x2): Occupancy and time slots histograms.
"""
function plot_frame3_occupancy(df::DataFrame)
    n_events = nrow(df)

    # Average occupancy histogram (as percentage)
    occ_pct = 100.0 .* df.avg_occupancy
    p1 = histogram(occ_pct, bins=30,
                   xlabel="Avg Occupancy (%)", ylabel="Events",
                   title="Pixel Occupancy", titlefontsize=10,
                   label="", color=:gray, bottom_margin=5Plots.mm)
    add_stats_legend!(p1, occ_pct; digits=2)

    # Number of time slots histogram
    p2 = histogram(df.n_time_slots, bins=30,
                   xlabel="N Time Slots", ylabel="Events",
                   title="Time Slots (2mm bins)", titlefontsize=10,
                   label="", color=:gray, bottom_margin=5Plots.mm)
    add_stats_legend!(p2, df.n_time_slots; digits=1)

    plt = plot(p1, p2, layout=(1, 2), size=(900, 400),
               plot_title="Occupancy Analysis (N=$n_events)", plot_titlefontsize=11)

    return plt
end

"""
    plot_frame4_pixel_signals(df::DataFrame)

Frame 4 (2x2): Pixel signal statistics histograms.
"""
function plot_frame4_pixel_signals(df::DataFrame)
    n_events = nrow(df)

    # Number of pixels above threshold (range 0-60)
    p1 = histogram(df.avg_n_pixels_above_threshold, bins=range(0, 60, length=31),
                   xlabel="N Pixels (>10 ions)", ylabel="Events",
                   title="Avg Pixels Above Threshold", titlefontsize=9,
                   label="", color=:gray, xlims=(0, 60), bottom_margin=5Plots.mm)
    add_stats_legend!(p1, df.avg_n_pixels_above_threshold; digits=1)

    # Max signal (range 0-2500)
    p2 = histogram(df.avg_max_signal, bins=range(0, 2500, length=31),
                   xlabel="Max Signal (ions)", ylabel="Events",
                   title="Avg Max Signal", titlefontsize=9,
                   label="", color=:gray, xlims=(0, 2500), bottom_margin=5Plots.mm)
    add_stats_legend!(p2, df.avg_max_signal; digits=0)

    # Min signal (range 0-100, filter out zeros)
    min_signals = filter(x -> x > 0, df.avg_min_signal)
    if !isempty(min_signals)
        p3 = histogram(min_signals, bins=range(0, 100, length=31),
                       xlabel="Min Signal (ions)", ylabel="Events",
                       title="Avg Min Signal (>10 ions)", titlefontsize=9,
                       label="", color=:gray, xlims=(0, 100), bottom_margin=5Plots.mm)
        add_stats_legend!(p3, min_signals; digits=1)
    else
        p3 = plot(title="Avg Min Signal\n(no data)", titlefontsize=9, legend=false)
    end

    # Mean signal (range 0-500)
    p4 = histogram(df.avg_mean_signal, bins=range(0, 500, length=31),
                   xlabel="Mean Signal (ions)", ylabel="Events",
                   title="Avg Mean Signal", titlefontsize=9,
                   label="", color=:gray, xlims=(0, 500), bottom_margin=5Plots.mm)
    add_stats_legend!(p4, df.avg_mean_signal; digits=1)

    plt = plot(p1, p2, p3, p4, layout=(2, 2), size=(900, 700),
               plot_title="Pixel Signal Statistics (N=$n_events)", plot_titlefontsize=11)

    return plt
end

# =============================================================================
# Summary Statistics
# =============================================================================

"""
    print_summary_statistics(df::DataFrame)

Print summary statistics for all variables.
"""
function print_summary_statistics(df::DataFrame)
    n_events = nrow(df)

    println()
    println("=" ^ 70)
    println("           SINGLE-TRACK STATISTICS SUMMARY")
    println("=" ^ 70)
    println("Number of single-track events: $n_events")
    println()

    # Define variables to summarize (exclude event_id)
    vars = [
        (:total_energy_keV, "Total Energy (keV)"),
        (:n_voxels, "N Voxels"),
        (:z_length_mm, "Z Length (mm)"),
        (:x_extent_mm, "X Extent (mm)"),
        (:y_extent_mm, "Y Extent (mm)"),
        (:avg_occupancy, "Avg Occupancy"),
        (:n_time_slots, "N Time Slots"),
        (:avg_n_pixels_above_threshold, "Avg N Pixels >10 ions"),
        (:avg_max_signal, "Avg Max Signal (ions)"),
        (:avg_min_signal, "Avg Min Signal (ions)"),
        (:avg_mean_signal, "Avg Mean Signal (ions)")
    ]

    println("Variable                        Mean          Std          Min          Max")
    println("-" ^ 70)

    for (col, label) in vars
        data = df[!, col]
        # Filter out zeros/NaNs for min signal
        if col == :avg_min_signal
            data = filter(x -> x > 0, data)
            if isempty(data)
                println("$(rpad(label, 28))  (no data)")
                continue
            end
        end

        m = mean(data)
        s = std(data)
        mi = minimum(data)
        ma = maximum(data)

        # Format based on magnitude
        if col == :avg_occupancy
            println("$(rpad(label, 28))  $(lpad(round(100*m, digits=2), 10))%  $(lpad(round(100*s, digits=2), 10))%  $(lpad(round(100*mi, digits=2), 10))%  $(lpad(round(100*ma, digits=2), 10))%")
        elseif col in (:n_voxels, :n_time_slots)
            println("$(rpad(label, 28))  $(lpad(round(m, digits=1), 10))  $(lpad(round(s, digits=1), 10))  $(lpad(Int(mi), 10))  $(lpad(Int(ma), 10))")
        else
            println("$(rpad(label, 28))  $(lpad(round(m, digits=1), 10))  $(lpad(round(s, digits=1), 10))  $(lpad(round(mi, digits=1), 10))  $(lpad(round(ma, digits=1), 10))")
        end
    end

    println("=" ^ 70)
end

# =============================================================================
# Main
# =============================================================================

function main()
    # Check arguments
    if length(ARGS) < 1
        println("Error: Missing required argument")
        println()
        println("Usage: julia single_track_stats.jl <results_dir> [options]")
        println()
        println("Arguments:")
        println("  results_dir       Directory containing single_track_stats.csv")
        println()
        println("Options:")
        println("  --no-save         Display only, don't save plots (default: save)")
        println("  --output-dir=X    Output directory for plots (default: results_dir)")
        println()
        println("Example:")
        println("  julia single_track_stats.jl bb0nu_voxel_2mm_maxd_4mm_evt_1_1000")
        exit(1)
    end

    # Parse arguments
    results_dir = ARGS[1]
    save_plots = true
    output_dir = results_dir

    for arg in ARGS[2:end]
        if arg == "--no-save"
            save_plots = false
        elseif startswith(arg, "--output-dir=")
            output_dir = split(arg, '=')[2]
        else
            println("Warning: Unknown argument: $arg")
        end
    end

    # Check if results directory exists
    if !isdir(results_dir)
        error("Results directory not found: $results_dir")
    end

    # Load CSV file
    csv_path = joinpath(results_dir, "single_track_stats.csv")
    if !isfile(csv_path)
        error("CSV file not found: $csv_path")
    end

    println("Loading: $csv_path")
    df = CSV.read(csv_path, DataFrame)
    println("Loaded $(nrow(df)) single-track events")

    if nrow(df) == 0
        println("No data to plot.")
        return
    end

    # Print summary statistics
    print_summary_statistics(df)

    # Generate plots
    println("\nGenerating plots...")

    plt1 = plot_frame1_energy_voxels(df)
    plt2 = plot_frame2_geometry(df)
    plt3 = plot_frame3_occupancy(df)
    plt4 = plot_frame4_pixel_signals(df)

    # Save or display plots
    if save_plots
        if !isdir(output_dir)
            mkpath(output_dir)
        end

        path1 = joinpath(output_dir, "single_track_energy_voxels.png")
        path2 = joinpath(output_dir, "single_track_geometry.png")
        path3 = joinpath(output_dir, "single_track_occupancy.png")
        path4 = joinpath(output_dir, "single_track_pixel_signals.png")

        savefig(plt1, path1)
        println("Saved: $path1")

        savefig(plt2, path2)
        println("Saved: $path2")

        savefig(plt3, path3)
        println("Saved: $path3")

        savefig(plt4, path4)
        println("Saved: $path4")
    end

    # Display plots
    display(plt1)
    println("\nPress Enter for geometry plots...")
    readline()

    display(plt2)
    println("\nPress Enter for occupancy plots...")
    readline()

    display(plt3)
    println("\nPress Enter for pixel signal plots...")
    readline()

    display(plt4)
    println("\nPress Enter to exit...")
    readline()

    println("\nDone!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
