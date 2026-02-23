#!/usr/bin/env julia

"""
ITACA Energy ROI Histograms.

Reads filtered event CSV files from event_stats_filtered/ and creates histograms
of main_track_energy_keV for each species (0nubb, bi214, tl208).

For each isotope, combines both σ=1mm and σ=10mm data.

Output:
    - 3 unsmeared histograms (one per isotope, combining both sigma values)
    - 3 smeared histograms (Gaussian smearing applied)
    - Combined plots

Usage:
    julia itaca_eroi.jl [options]

Options:
    --input-dir=X   Input directory with filtered CSVs (default: event_stats_filtered)
    --output-dir=X  Output directory for plots (default: eroi_plots)
    --sigma=X       Gaussian smearing sigma in keV (default: 25)
    --help, -h      Show this help message
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using CSV
using DataFrames
using Plots
using Printf
using Dates
using Random

# =============================================================================
# Configuration
# =============================================================================

const ISOTOPES = ["0nubb", "bi214", "tl208"]
const SIGMA_CONFIGS = [
    (label="1mm", pattern="sigmat_1mm"),
    (label="10mm", pattern="sigmat_10mm")
]

# =============================================================================
# Data Loading
# =============================================================================

"""
    find_csv_files(input_dir::String, isotope::String, sigma_pattern::String)

Find all CSV files matching the isotope and sigma pattern.
Returns a vector of file paths.
"""
function find_csv_files(input_dir::String, isotope::String, sigma_pattern::String)
    files = String[]
    for f in readdir(input_dir)
        if endswith(f, ".csv") && startswith(f, isotope) && contains(f, sigma_pattern)
            push!(files, joinpath(input_dir, f))
        end
    end
    return files
end

"""
    load_and_combine_energies(input_dir::String, isotope::String, sigma_pattern::String)

Load all CSV files for a given isotope/sigma combination and extract main_track_energy_keV.
Returns a vector of energy values.
"""
function load_and_combine_energies(input_dir::String, isotope::String, sigma_pattern::String)
    files = find_csv_files(input_dir, isotope, sigma_pattern)

    if isempty(files)
        return Float64[]
    end

    energies = Float64[]
    for f in files
        df = CSV.read(f, DataFrame; comment="#")
        if "main_track_energy_keV" in names(df)
            append!(energies, df.main_track_energy_keV)
        else
            @warn "Column main_track_energy_keV not found in $f"
        end
    end

    return energies
end

"""
    load_all_energies_for_isotope(input_dir::String, isotope::String; energy_cut::Union{Float64, Nothing}=nothing)

Load all CSV files for a given isotope (both σ=1mm and σ=10mm) and combine.
If energy_cut is provided, only include events with energy > energy_cut.
Returns a vector of energy values.
"""
function load_all_energies_for_isotope(input_dir::String, isotope::String; energy_cut::Union{Float64, Nothing}=nothing)
    all_energies = Float64[]

    for sigma_cfg in SIGMA_CONFIGS
        energies = load_and_combine_energies(input_dir, isotope, sigma_cfg.pattern)
        if !isempty(energies)
            println("    $(sigma_cfg.label): $(length(energies)) events")
            append!(all_energies, energies)
        end
    end

    # Apply energy cut if specified
    if energy_cut !== nothing && !isempty(all_energies)
        n_before = length(all_energies)
        all_energies = filter(e -> e > energy_cut, all_energies)
        n_after = length(all_energies)
        println("    Applied cut E > $energy_cut keV: $n_before → $n_after events")
    end

    return all_energies
end

# =============================================================================
# Smearing
# =============================================================================

using Statistics: mean, std

"""
    apply_gaussian_smearing(energies::Vector{Float64}, sigma::Float64)

Apply Gaussian smearing to energy values.
Each energy value E becomes E + N(0, sigma).
"""
function apply_gaussian_smearing(energies::Vector{Float64}, sigma::Float64)
    return energies .+ sigma .* randn(length(energies))
end

# =============================================================================
# Plotting
# =============================================================================

"""
    format_isotope_display(isotope::String)

Format isotope name for display with proper Unicode.
"""
function format_isotope_display(isotope::String)
    if isotope == "0nubb"
        return "0νββ"
    elseif isotope == "bi214"
        return "²¹⁴Bi"
    elseif isotope == "tl208"
        return "²⁰⁸Tl"
    else
        return isotope
    end
end

"""
    get_isotope_color(isotope::String)

Get color for each isotope.
"""
function get_isotope_color(isotope::String)
    if isotope == "0nubb"
        return :blue
    elseif isotope == "bi214"
        return :red
    elseif isotope == "tl208"
        return :green
    else
        return :gray
    end
end

"""
    plot_unsmeared_histogram(energies::Vector{Float64}, isotope::String)

Create a histogram of unsmeared main track energies.
"""
function plot_unsmeared_histogram(energies::Vector{Float64}, isotope::String)
    if isempty(energies)
        return nothing
    end

    isotope_display = format_isotope_display(isotope)
    n_events = length(energies)
    mean_e = mean(energies)
    std_e = std(energies)

    plt = histogram(energies,
        bins=50,
        xlabel="Main Track Energy (keV)",
        ylabel="Counts",
        title="$isotope_display - Unsmeared (N=$n_events)",
        label="",
        fillalpha=0.7,
        linewidth=1,
        color=get_isotope_color(isotope),
        size=(600, 400),
        xlims=(2380, 2520),
        titlefontsize=12,
        guidefontsize=10,
        tickfontsize=9
    )

    # Add statistics annotation
    annotate!(plt, 2390, 0.9 * ylims(plt)[2],
        text(@sprintf("μ = %.1f keV\nσ = %.1f keV", mean_e, std_e), 9, :left))

    return plt
end

"""
    plot_smeared_histogram(energies::Vector{Float64}, isotope::String, smear_sigma::Float64)

Create a histogram of smeared main track energies.
"""
function plot_smeared_histogram(energies::Vector{Float64}, isotope::String, smear_sigma::Float64)
    if isempty(energies)
        return nothing
    end

    # Apply smearing
    smeared = apply_gaussian_smearing(energies, smear_sigma)

    isotope_display = format_isotope_display(isotope)
    n_events = length(smeared)
    mean_e = mean(smeared)
    std_e = std(smeared)

    plt = histogram(smeared,
        bins=50,
        xlabel="Energy (keV)",
        ylabel="Counts",
        title="$isotope_display - Smeared (σ=$(smear_sigma) keV, N=$n_events)",
        label="",
        fillalpha=0.7,
        linewidth=1,
        color=get_isotope_color(isotope),
        size=(600, 400),
        xlims=(2300, 2600),
        titlefontsize=12,
        guidefontsize=10,
        tickfontsize=9
    )

    # Add statistics annotation
    annotate!(plt, 2320, 0.9 * ylims(plt)[2],
        text(@sprintf("μ = %.1f keV\nσ = %.1f keV", mean_e, std_e), 9, :left))

    return plt
end

"""
    create_unsmeared_combined_plot(histograms::Dict, output_path::String)

Create a 1x3 grid of unsmeared histograms.
"""
function create_unsmeared_combined_plot(histograms::Dict, output_path::String)
    plots_grid = []

    for isotope in ISOTOPES
        if haskey(histograms, isotope) && histograms[isotope] !== nothing
            push!(plots_grid, histograms[isotope])
        else
            push!(plots_grid, plot(title="No data", grid=false, axis=false))
        end
    end

    combined = plot(plots_grid...,
        layout=(1, 3),
        size=(1500, 450),
        margin=5Plots.mm,
        plot_title="Unsmeared Energy Data (Single-Track Events)"
    )

    savefig(combined, output_path)
    println("Saved unsmeared combined plot: $output_path")

    return combined
end

"""
    create_smeared_combined_plot(histograms::Dict, smear_sigma::Float64, output_path::String)

Create a 1x3 grid of smeared histograms.
"""
function create_smeared_combined_plot(histograms::Dict, smear_sigma::Float64, output_path::String)
    plots_grid = []

    for isotope in ISOTOPES
        if haskey(histograms, isotope) && histograms[isotope] !== nothing
            push!(plots_grid, histograms[isotope])
        else
            push!(plots_grid, plot(title="No data", grid=false, axis=false))
        end
    end

    combined = plot(plots_grid...,
        layout=(1, 3),
        size=(1500, 450),
        margin=5Plots.mm,
        plot_title="Energy Distribution for Signal and Background (σ_smear=$(smear_sigma) keV)"
    )

    savefig(combined, output_path)
    println("Saved smeared combined plot: $output_path")

    return combined
end

"""
    create_overlaid_smeared_plot(all_energies::Dict, smear_sigma::Float64, output_path::String)

Create a single plot with all three smeared histograms overlaid.
Uses different colors and alpha values for visibility.
"""
function create_overlaid_smeared_plot(all_energies::Dict, smear_sigma::Float64, output_path::String)
    plt = plot(
        xlabel="Energy (keV)",
        ylabel="Counts (normalized)",
        title="Signal and Background Energy Distribution (σ=$(smear_sigma) keV)",
        size=(800, 500),
        xlims=(2300, 2600),
        legend=:topright,
        titlefontsize=12,
        guidefontsize=10,
        tickfontsize=9
    )

    # Plot each isotope with different alpha and colors
    isotope_settings = [
        ("0nubb", :blue, 0.4, "0νββ Signal"),
        ("bi214", :red, 0.3, "²¹⁴Bi Background"),
        ("tl208", :green, 0.3, "²⁰⁸Tl Background")
    ]

    for (isotope, color, alpha, label) in isotope_settings
        if haskey(all_energies, isotope) && !isempty(all_energies[isotope])
            smeared = apply_gaussian_smearing(all_energies[isotope], smear_sigma)
            histogram!(plt, smeared,
                bins=50,
                fillalpha=alpha,
                linewidth=1,
                color=color,
                label=label,
                normalize=:probability
            )
        end
    end

    savefig(plt, output_path)
    println("Saved overlaid smeared plot: $output_path")

    return plt
end

# =============================================================================
# Main
# =============================================================================

function main()
    # Default values
    input_dir = "event_stats_filtered"
    output_dir = "eroi_plots"
    smear_sigma = 25.0  # keV

    # Parse arguments
    for arg in ARGS
        if startswith(arg, "--input-dir=")
            input_dir = String(split(arg, '=')[2])
        elseif startswith(arg, "--output-dir=")
            output_dir = String(split(arg, '=')[2])
        elseif startswith(arg, "--sigma=")
            smear_sigma = parse(Float64, split(arg, '=')[2])
        elseif arg == "--help" || arg == "-h"
            println("""
Usage: julia itaca_eroi.jl [options]

Options:
    --input-dir=X   Input directory with filtered CSVs (default: event_stats_filtered)
    --output-dir=X  Output directory for plots (default: eroi_plots)
    --sigma=X       Gaussian smearing sigma in keV (default: 25)
    --help, -h      Show this help message

Output:
    Creates energy histograms for 0nubb, bi214, tl208:
    - 3 unsmeared plots (combining σ=1mm and σ=10mm data)
    - 3 smeared plots (Gaussian smearing applied)
    - Combined plot views
""")
            return
        else
            @warn "Unknown argument: $arg"
        end
    end

    # Resolve paths relative to script directory
    script_dir = dirname(@__FILE__)
    full_input_dir = isabspath(input_dir) ? input_dir : joinpath(script_dir, input_dir)
    full_output_dir = isabspath(output_dir) ? output_dir : joinpath(script_dir, output_dir)

    # Validate input directory
    if !isdir(full_input_dir)
        error("Input directory not found: $full_input_dir")
    end

    # Create output directory
    if !isdir(full_output_dir)
        mkpath(full_output_dir)
        println("Created output directory: $full_output_dir")
    end

    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")

    println("=" ^ 60)
    println("ITACA Energy ROI Histograms")
    println("=" ^ 60)
    println("Generated: $timestamp")
    println("Input directory: $full_input_dir")
    println("Output directory: $full_output_dir")
    println("Smearing sigma: $smear_sigma keV")
    println()

    # Load all data for each isotope
    # Apply E > 2440 keV cut for 0nubb only
    println("Loading data...")
    all_energies = Dict{String, Vector{Float64}}()

    for isotope in ISOTOPES
        println("  $isotope:")
        # Apply energy cut only for 0nubb signal
        energy_cut = isotope == "0nubb" ? 2440.0 : nothing
        energies = load_all_energies_for_isotope(full_input_dir, isotope; energy_cut=energy_cut)
        all_energies[isotope] = energies
        println("    Total: $(length(energies)) events")
    end

    # Create unsmeared histograms
    println("\nCreating unsmeared histograms...")
    unsmeared_hists = Dict{String, Any}()

    for isotope in ISOTOPES
        energies = all_energies[isotope]
        if !isempty(energies)
            plt = plot_unsmeared_histogram(energies, isotope)
            unsmeared_hists[isotope] = plt

            # Save individual histogram
            individual_path = joinpath(full_output_dir, "unsmeared_$(isotope).png")
            savefig(plt, individual_path)
            println("  Saved: $individual_path")
        else
            unsmeared_hists[isotope] = nothing
        end
    end

    # Create unsmeared combined plot
    unsmeared_combined_path = joinpath(full_output_dir, "unsmeared_combined.png")
    create_unsmeared_combined_plot(unsmeared_hists, unsmeared_combined_path)

    # Create smeared histograms
    println("\nCreating smeared histograms (σ=$smear_sigma keV)...")
    smeared_hists = Dict{String, Any}()

    for isotope in ISOTOPES
        energies = all_energies[isotope]
        if !isempty(energies)
            plt = plot_smeared_histogram(energies, isotope, smear_sigma)
            smeared_hists[isotope] = plt

            # Save individual histogram
            individual_path = joinpath(full_output_dir, "smeared_$(isotope)_sigma$(smear_sigma)keV.png")
            savefig(plt, individual_path)
            println("  Saved: $individual_path")
        else
            smeared_hists[isotope] = nothing
        end
    end

    # Create smeared combined plot (side by side)
    smeared_combined_path = joinpath(full_output_dir, "smeared_combined_sigma$(smear_sigma)keV.png")
    create_smeared_combined_plot(smeared_hists, smear_sigma, smeared_combined_path)

    # Create overlaid smeared plot (all three on same axes)
    overlaid_path = joinpath(full_output_dir, "smeared_overlaid_sigma$(smear_sigma)keV.png")
    create_overlaid_smeared_plot(all_energies, smear_sigma, overlaid_path)

    # ==========================================================================
    # ROI Cut Analysis
    # ==========================================================================
    println("\n" * "=" ^ 60)
    println("ROI Cut Analysis")
    println("=" ^ 60)

    # Define ROI: Qbb ± 1 FWHM, where FWHM = 2.35 * sigma
    Qbb = 2457.0  # keV
    fwhm = 2.35 * smear_sigma
    roi_min = Qbb - fwhm
    roi_max = Qbb + fwhm

    println("Qββ = $Qbb keV")
    println("FWHM = 2.35 × σ = $(round(fwhm, digits=2)) keV")
    println("ROI = [$(@sprintf("%.2f", roi_min)), $(@sprintf("%.2f", roi_max))] keV")
    println()

    # Apply smearing and count events in ROI for each isotope
    roi_results = Dict{String, NamedTuple{(:total, :in_roi, :efficiency), Tuple{Int, Int, Float64}}}()

    for isotope in ISOTOPES
        energies = all_energies[isotope]
        if isempty(energies)
            continue
        end

        # Apply smearing
        smeared = apply_gaussian_smearing(energies, smear_sigma)
        total = length(smeared)

        # Count events in ROI
        in_roi = count(e -> roi_min <= e <= roi_max, smeared)
        efficiency = in_roi / total

        roi_results[isotope] = (total=total, in_roi=in_roi, efficiency=efficiency)

        isotope_display = format_isotope_display(isotope)
        println("  $isotope_display:")
        println("    Total events: $total")
        println("    Events in ROI: $in_roi")
        println("    Efficiency: $(@sprintf("%.4f", efficiency)) ($(@sprintf("%.2f", efficiency * 100))%)")
    end

    # Write ROI efficiency table (LaTeX)
    roi_tex_path = joinpath(full_output_dir, "roi_efficiency_sigma$(smear_sigma)keV.tex")
    write_roi_table_tex(roi_results, smear_sigma, Qbb, fwhm, roi_tex_path)

    # Write ROI efficiency table (CSV)
    roi_csv_path = joinpath(full_output_dir, "roi_efficiency_sigma$(smear_sigma)keV.csv")
    write_roi_table_csv(roi_results, smear_sigma, Qbb, fwhm, roi_csv_path)

    println()
    println("=" ^ 60)
    println("Done!")
    println("=" ^ 60)
end

# =============================================================================
# ROI Table Output
# =============================================================================

"""
    write_roi_table_tex(results::Dict, smear_sigma::Float64, Qbb::Float64, fwhm::Float64, output_path::String)

Write ROI efficiency table in LaTeX format.
"""
function write_roi_table_tex(results::Dict, smear_sigma::Float64, Qbb::Float64, fwhm::Float64, output_path::String)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    roi_min = Qbb - fwhm
    roi_max = Qbb + fwhm

    open(output_path, "w") do f
        println(f, "% Generated: $timestamp")
        println(f, "\\documentclass{article}")
        println(f, "\\usepackage{booktabs}")
        println(f, "\\usepackage{siunitx}")
        println(f, "\\begin{document}")
        println(f)
        println(f, "\\begin{table}[htbp]")
        println(f, "\\centering")
        println(f, "\\caption{ROI Efficiency (\$\\sigma\$ = $smear_sigma keV, FWHM = $(@sprintf("%.2f", fwhm)) keV, ROI = [$(@sprintf("%.1f", roi_min)), $(@sprintf("%.1f", roi_max))] keV)}")
        println(f, "\\label{tab:roi_efficiency}")
        println(f, "\\begin{tabular}{lrrr}")
        println(f, "\\toprule")
        println(f, "Species & Total Events & Events in ROI & Efficiency \\\\")
        println(f, "\\midrule")

        for isotope in ISOTOPES
            if !haskey(results, isotope)
                continue
            end
            r = results[isotope]

            # Format isotope name for LaTeX
            label = if isotope == "0nubb"
                "\$0\\nu\\beta\\beta\$"
            elseif isotope == "bi214"
                "\$^{214}\$Bi"
            elseif isotope == "tl208"
                "\$^{208}\$Tl"
            else
                isotope
            end

            @printf(f, "%s & %d & %d & %.4f \\\\\n", label, r.total, r.in_roi, r.efficiency)
        end

        println(f, "\\bottomrule")
        println(f, "\\end{tabular}")
        println(f, "\\end{table}")
        println(f)
        println(f, "\\end{document}")
    end

    println("  Saved LaTeX table: $output_path")
end

"""
    write_roi_table_csv(results::Dict, smear_sigma::Float64, Qbb::Float64, fwhm::Float64, output_path::String)

Write ROI efficiency table in CSV format.
"""
function write_roi_table_csv(results::Dict, smear_sigma::Float64, Qbb::Float64, fwhm::Float64, output_path::String)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    roi_min = Qbb - fwhm
    roi_max = Qbb + fwhm

    open(output_path, "w") do f
        println(f, "# Generated: $timestamp")
        println(f, "# Smearing sigma: $smear_sigma keV")
        println(f, "# FWHM: $(@sprintf("%.2f", fwhm)) keV")
        println(f, "# ROI: [$(@sprintf("%.2f", roi_min)), $(@sprintf("%.2f", roi_max))] keV")
        println(f, "species,total_events,events_in_roi,efficiency")

        for isotope in ISOTOPES
            if !haskey(results, isotope)
                continue
            end
            r = results[isotope]
            @printf(f, "%s,%d,%d,%.6f\n", isotope, r.total, r.in_roi, r.efficiency)
        end
    end

    println("  Saved CSV table: $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
