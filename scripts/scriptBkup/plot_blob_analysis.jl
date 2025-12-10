#!/usr/bin/env julia

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

# Now load packages
using ArgParse
using CSV
using DataFrames
using Plots
using StatsPlots
using Statistics
using Printf

"""
    parse_commandline()

Parse command line arguments for plotting script.
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input", "-i"
            help = "Input CSV file path"
            arg_type = String
            default = "blob_analysis.csv"
        "--output", "-o"
            help = "Output plot file (PNG/PDF/SVG)"
            arg_type = String
            default = "blob_analysis_plots.png"
        "--nbins"
            help = "Number of bins for histograms"
            arg_type = Int
            default = 30
        "--figsize"
            help = "Figure size as 'width,height' in pixels"
            arg_type = String
            default = "1200,1000"
        "--title"
            help = "Main title for the figure"
            arg_type = String
            default = "Blob Analysis Results"
        "--eblob-cut"
            help = "Energy cut for blob 2 in keV (shown as horizontal line)"
            arg_type = Float64
            default = 600.0
    end

    return parse_args(s)
end

"""
    create_2x2_plots(df::DataFrame; nbins=30, title="Blob Analysis Results", eblob_cut=600.0)

Create a 2x2 subplot layout with:
- Top left: Track energy histogram
- Top right: Track length histogram
- Bottom left: Scatter plot of eblob1 vs eblob2 with energy cut line
- Bottom right: Scatter plot of nblob1 vs nblob2
"""
function create_2x2_plots(df::DataFrame; nbins=30, title="Blob Analysis Results", eblob_cut=600.0)

    # Create subplots
    p1 = histogram(df.etrk,
        bins = nbins,
        xlabel = "Track Energy (keV)",
        ylabel = "Counts",
        title = "Track Energy Distribution",
        label = nothing,
        fillcolor = :steelblue,
        linecolor = :black,
        alpha = 0.7,
        xlims = (minimum(df.etrk) * 0.95, maximum(df.etrk) * 1.05)
    )

    # Add mean and std annotation
    mean_etrk = mean(df.etrk)
    std_etrk = std(df.etrk)
    annotate!(p1,
        :topright,
        text(@sprintf("μ = %.1f keV\nσ = %.1f keV", mean_etrk, std_etrk), 10, :left)
    )

    p2 = histogram(df.trkl,
        bins = nbins,
        xlabel = "Track Length (mm)",
        ylabel = "Counts",
        title = "Track Length Distribution",
        label = nothing,
        fillcolor = :darkorange,
        linecolor = :black,
        alpha = 0.7,
        xlims = (minimum(df.trkl) * 0.95, maximum(df.trkl) * 1.05)
    )

    # Add mean and std annotation
    mean_trkl = mean(df.trkl)
    std_trkl = std(df.trkl)
    annotate!(p2,
        :topright,
        text(@sprintf("μ = %.1f mm\nσ = %.1f mm", mean_trkl, std_trkl), 10, :left)
    )

    p3 = scatter(df.eblob1, df.eblob2,
        xlabel = "Blob 1 Energy (keV)",
        ylabel = "Blob 2 Energy (keV)",
        title = "Blob Energy Correlation",
        label = nothing,
        markersize = 4,
        markercolor = :purple,
        markeralpha = 0.6,
        markerstrokewidth = 0.5,
        markerstrokecolor = :black,
        xlims = (0, maximum([maximum(df.eblob1), maximum(df.eblob2)]) * 1.1),
        ylims = (0, maximum([maximum(df.eblob1), maximum(df.eblob2)]) * 1.1),
        aspect_ratio = :equal,
        grid = true,
        gridstyle = :dot,
        gridalpha = 0.3
    )

    # Add diagonal reference line
    max_energy = maximum([maximum(df.eblob1), maximum(df.eblob2)])
    plot!(p3, [0, max_energy], [0, max_energy],
        line = :dash,
        linecolor = :red,
        linewidth = 1,
        label = "y = x",
        alpha = 0.5
    )

    # Add horizontal line at eblob_cut level
    plot!(p3, [0, max_energy], [eblob_cut, eblob_cut],
        line = :dash,
        linecolor = :blue,
        linewidth = 2,
        label = "Cut: $(eblob_cut) keV",
        alpha = 0.7
    )

    # Calculate and add correlation coefficient
    if nrow(df) > 1
        corr_energy = cor(df.eblob1, df.eblob2)
        annotate!(p3,
            :bottomright,
            text(@sprintf("r = %.3f", corr_energy), 10, :right)
        )
    end

    p4 = scatter(df.nblob1, df.nblob2,
        xlabel = "Blob 1 Voxel Count",
        ylabel = "Blob 2 Voxel Count",
        title = "Blob Voxel Count Correlation",
        label = nothing,
        markersize = 4,
        markercolor = :green,
        markeralpha = 0.6,
        markerstrokewidth = 0.5,
        markerstrokecolor = :black,
        xlims = (0, maximum([maximum(df.nblob1), maximum(df.nblob2)]) * 1.1),
        ylims = (0, maximum([maximum(df.nblob1), maximum(df.nblob2)]) * 1.1),
        aspect_ratio = :equal,
        grid = true,
        gridstyle = :dot,
        gridalpha = 0.3
    )

    # Add diagonal reference line
    max_count = maximum([maximum(df.nblob1), maximum(df.nblob2)])
    plot!(p4, [0, max_count], [0, max_count],
        line = :dash,
        linecolor = :red,
        linewidth = 1,
        label = "y = x",
        alpha = 0.5
    )

    # Calculate and add correlation coefficient
    if nrow(df) > 1
        corr_count = cor(df.nblob1, df.nblob2)
        annotate!(p4,
            :bottomright,
            text(@sprintf("r = %.3f", corr_count), 10, :right)
        )
    end

    # Combine all plots into 2x2 layout
    plot(p1, p2, p3, p4,
        layout = (2, 2),
        size = (1200, 1000),
        plot_title = title,
        plot_titlefontsize = 16,
        margin = 5Plots.mm
    )
end

"""
    main()

Main function for plotting blob analysis results.
"""
function main()
    # Parse command line arguments
    args = parse_commandline()

    # Parse figure size
    figsize = parse.(Int, split(args["figsize"], ","))
    if length(figsize) != 2
        figsize = [1200, 1000]
        println("Warning: Invalid figsize format. Using default 1200,1000")
    end

    println("="^60)
    println("Blob Analysis Plotting Script")
    println("="^60)
    println("Input file: $(args["input"])")
    println("Output file: $(args["output"])")
    println("Number of bins: $(args["nbins"])")
    println("Figure size: $(figsize[1])x$(figsize[2])")
    println("Blob 2 energy cut: $(args["eblob-cut"]) keV")
    println("="^60)

    # Check if input file exists
    if !isfile(args["input"])
        error("Input file $(args["input"]) does not exist!")
    end

    # Read CSV file
    println("\nReading data from $(args["input"])...")
    df = CSV.read(args["input"], DataFrame)
    println("Loaded $(nrow(df)) rows")

    # Print summary statistics
    println("\n" * "="^60)
    println("Data Summary:")
    println("="^60)
    println("Track Energy:")
    println("  Mean: $(round(mean(df.etrk), digits=1)) keV")
    println("  Std:  $(round(std(df.etrk), digits=1)) keV")
    println("  Min:  $(round(minimum(df.etrk), digits=1)) keV")
    println("  Max:  $(round(maximum(df.etrk), digits=1)) keV")

    println("\nTrack Length:")
    println("  Mean: $(round(mean(df.trkl), digits=1)) mm")
    println("  Std:  $(round(std(df.trkl), digits=1)) mm")
    println("  Min:  $(round(minimum(df.trkl), digits=1)) mm")
    println("  Max:  $(round(maximum(df.trkl), digits=1)) mm")

    println("\nBlob 1 Energy:")
    println("  Mean: $(round(mean(df.eblob1), digits=1)) keV")
    println("  Std:  $(round(std(df.eblob1), digits=1)) keV")

    println("\nBlob 2 Energy:")
    println("  Mean: $(round(mean(df.eblob2), digits=1)) keV")
    println("  Std:  $(round(std(df.eblob2), digits=1)) keV")

    println("\nConfidence:")
    println("  Mean: $(round(mean(df.confidence), digits=3))")
    println("  Std:  $(round(std(df.confidence), digits=3))")

    # Check if 'select' column exists and print selection statistics
    if "select" in names(df)
        println("\nSelection Statistics (eblob_cut = $(args["eblob-cut"]) keV):")
        n_selected = sum(df.select)
        n_total = nrow(df)
        println("  Tracks passing cut: $(n_selected)/$(n_total) ($(round(100*n_selected/n_total, digits=1))%)")
    end

    println("="^60)

    # Create plots
    println("\nCreating plots...")
    p = create_2x2_plots(df,
        nbins = args["nbins"],
        title = args["title"],
        eblob_cut = args["eblob-cut"]
    )

    # Adjust size
    plot!(p, size = (figsize[1], figsize[2]))

    # Save plot
    savefig(p, args["output"])
    println("Plot saved to: $(args["output"])")

    # Also create individual plots if needed
    if get(ENV, "CREATE_INDIVIDUAL_PLOTS", "false") == "true"
        println("\nCreating individual plots...")

        # Extract base name and extension
        base_name = splitext(args["output"])[1]
        ext = splitext(args["output"])[2]

        # Save individual plots
        p1 = histogram(df.etrk, bins=args["nbins"], title="Track Energy", xlabel="Energy (keV)", label=nothing)
        savefig(p1, "$(base_name)_track_energy$(ext)")

        p2 = histogram(df.trkl, bins=args["nbins"], title="Track Length", xlabel="Length (mm)", label=nothing)
        savefig(p2, "$(base_name)_track_length$(ext)")

        p3 = scatter(df.eblob1, df.eblob2, title="Blob Energies", xlabel="Blob 1 (keV)", ylabel="Blob 2 (keV)", label=nothing)
        savefig(p3, "$(base_name)_blob_energy_scatter$(ext)")

        p4 = scatter(df.nblob1, df.nblob2, title="Blob Voxel Counts", xlabel="Blob 1 Count", ylabel="Blob 2 Count", label=nothing)
        savefig(p4, "$(base_name)_blob_count_scatter$(ext)")

        println("Individual plots saved with prefix: $(base_name)_*$(ext)")
    end

    println("\n✓ Done!")
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end