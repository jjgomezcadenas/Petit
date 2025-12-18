#!/usr/bin/env julia

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)

using DataFrames
using CSV
using Plots
using ArgParse
using Statistics

#=============================================================================
# Plotting Functions
=============================================================================#

"""
    plot_variable_vs_drift(df, var_col, ylabel; title=nothing, aggregate=false)

Plot a variable vs drift length for all events.
If aggregate=true, show mean ± std instead of individual events.
"""
function plot_variable_vs_drift(df, var_col::Symbol, ylabel::String;
                                 title::String="", aggregate::Bool=false)
    # Filter out NaN values for this variable
    df_valid = filter(row -> !isnan(row[var_col]), df)

    if nrow(df_valid) == 0
        println("  Warning: No valid data for $var_col")
        return plot(title="No valid data for $var_col")
    end

    events = unique(df_valid.event)

    if aggregate && length(events) > 1
        # Group by drift length and compute mean ± std (filtering NaN)
        grouped = combine(groupby(df_valid, :ldrft_cm),
                         var_col => (x -> mean(filter(!isnan, x))) => :mean,
                         var_col => (x -> std(filter(!isnan, x))) => :std)

        plt = plot(grouped.ldrft_cm, grouped.mean,
                   ribbon=grouped.std,
                   fillalpha=0.3,
                   xlabel="Drift length (cm)",
                   ylabel=ylabel,
                   title=title,
                   marker=:circle,
                   markersize=5,
                   linewidth=2,
                   label="Mean ± σ",
                   legend=:topright)
    else
        plt = plot(xlabel="Drift length (cm)",
                   ylabel=ylabel,
                   title=title,
                   legend=:topright)

        for evt in events
            evt_data = filter(row -> row.event == evt && !isnan(row[var_col]), df_valid)
            if nrow(evt_data) > 0
                plot!(plt, evt_data.ldrft_cm, evt_data[!, var_col],
                      marker=:circle,
                      markersize=5,
                      linewidth=2,
                      label="Event $evt")
            end
        end
    end

    return plt
end

"""
    create_all_plots(df, outdir; aggregate=false)

Create all variable plots and save to outdir.
"""
function create_all_plots(df, outdir; aggregate::Bool=false)
    mkpath(outdir)

    # Define base variables to plot (always present)
    variables = [
        (:d1_mm, "d₁ (mm)", "Distance d₁ vs Drift Length"),
        (:d2_mm, "d₂ (mm)", "Distance d₂ vs Drift Length"),
        (:Eb1_keV, "Eb₁ (keV)", "Blob Energy Eb₁ vs Drift Length"),
        (:Eb2_keV, "Eb₂ (keV)", "Blob Energy Eb₂ vs Drift Length"),
        (:asymmetry, "Asymmetry", "Blob Asymmetry vs Drift Length"),
    ]

    # Add peak variables if they exist in the CSV
    has_peaks = hasproperty(df, :peak1_left)
    if has_peaks
        push!(variables, (:peak1_left, "Peak1 Left", "Peak1 Left Edge vs Drift Length"))
        push!(variables, (:peak1_right, "Peak1 Right", "Peak1 Right Edge vs Drift Length"))
        push!(variables, (:peak1_prom, "Peak1 Prominence", "Peak1 Prominence vs Drift Length"))
        push!(variables, (:peak2_left, "Peak2 Left", "Peak2 Left Edge vs Drift Length"))
        push!(variables, (:peak2_right, "Peak2 Right", "Peak2 Right Edge vs Drift Length"))
        push!(variables, (:peak2_prom, "Peak2 Prominence", "Peak2 Prominence vs Drift Length"))
    end

    # Add derived variables
    df.d_total = df.d1_mm .+ df.d2_mm
    df.E_total = df.Eb1_keV .+ df.Eb2_keV

    push!(variables, (:d_total, "d₁ + d₂ (mm)", "Total Distance vs Drift Length"))
    push!(variables, (:E_total, "Eb₁ + Eb₂ (keV)", "Total Blob Energy vs Drift Length"))

    if has_peaks
        df.peak1_width = df.peak1_right .- df.peak1_left
        df.peak2_width = df.peak2_right .- df.peak2_left
        push!(variables, (:peak1_width, "Peak1 Width", "Peak1 Width vs Drift Length"))
        push!(variables, (:peak2_width, "Peak2 Width", "Peak2 Width vs Drift Length"))
    end

    suffix = aggregate ? "_aggregate" : ""

    for (col, ylabel, title) in variables
        if hasproperty(df, col)
            println("Plotting: $col")
            plt = plot_variable_vs_drift(df, col, ylabel; title=title, aggregate=aggregate)
            filename = joinpath(outdir, "$(col)_vs_drift$(suffix).png")
            savefig(plt, filename)
            println("  Saved: $filename")
        end
    end

    # Create a combined plot with all distances
    println("Creating combined distance plot...")
    plt = plot(xlabel="Drift length (cm)",
               ylabel="Distance (mm)",
               title="Distances vs Drift Length",
               legend=:topright)

    events = unique(df.event)
    colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :brown]

    for (i, evt) in enumerate(events)
        evt_data = filter(row -> row.event == evt, df)
        c = colors[mod1(i, length(colors))]
        plot!(plt, evt_data.ldrft_cm, evt_data.d1_mm,
              marker=:circle, color=c, linestyle=:solid,
              label="Evt $evt d₁")
        plot!(plt, evt_data.ldrft_cm, evt_data.d2_mm,
              marker=:square, color=c, linestyle=:dash,
              label="Evt $evt d₂")
    end
    savefig(plt, joinpath(outdir, "distances_combined$(suffix).png"))

    # Create a combined plot with blob energies
    println("Creating combined blob energy plot...")
    plt = plot(xlabel="Drift length (cm)",
               ylabel="Energy (keV)",
               title="Blob Energies vs Drift Length",
               legend=:topright)

    for (i, evt) in enumerate(events)
        evt_data = filter(row -> row.event == evt, df)
        c = colors[mod1(i, length(colors))]
        plot!(plt, evt_data.ldrft_cm, evt_data.Eb1_keV,
              marker=:circle, color=c, linestyle=:solid,
              label="Evt $evt Eb₁")
        plot!(plt, evt_data.ldrft_cm, evt_data.Eb2_keV,
              marker=:square, color=c, linestyle=:dash,
              label="Evt $evt Eb₂")
    end
    savefig(plt, joinpath(outdir, "blob_energies_combined$(suffix).png"))

    # Create combined peak plots if peak columns exist
    if has_peaks
        # Combined plot for peak prominences
        println("Creating combined peak prominence plot...")
        plt = plot(xlabel="Drift length (cm)",
                   ylabel="Prominence (normalized)",
                   title="Peak Prominences vs Drift Length",
                   legend=:topright)

        for (i, evt) in enumerate(events)
            evt_data = filter(row -> row.event == evt && row.peak1_prom > 0, df)
            c = colors[mod1(i, length(colors))]
            if nrow(evt_data) > 0
                plot!(plt, evt_data.ldrft_cm, evt_data.peak1_prom,
                      marker=:circle, color=c, linestyle=:solid,
                      label="Evt $evt Peak1")
            end
            evt_data2 = filter(row -> row.event == evt && row.peak2_prom > 0, df)
            if nrow(evt_data2) > 0
                plot!(plt, evt_data2.ldrft_cm, evt_data2.peak2_prom,
                      marker=:square, color=c, linestyle=:dash,
                      label="Evt $evt Peak2")
            end
        end
        savefig(plt, joinpath(outdir, "peak_proms_combined$(suffix).png"))

        # Combined plot for peak left edges
        println("Creating combined peak left edges plot...")
        plt = plot(xlabel="Drift length (cm)",
                   ylabel="Left edge (index)",
                   title="Peak Left Edges vs Drift Length",
                   legend=:topright)

        for (i, evt) in enumerate(events)
            evt_data = filter(row -> row.event == evt && row.peak1_left > 0, df)
            c = colors[mod1(i, length(colors))]
            if nrow(evt_data) > 0
                plot!(plt, evt_data.ldrft_cm, evt_data.peak1_left,
                      marker=:circle, color=c, linestyle=:solid,
                      label="Evt $evt Peak1")
            end
            evt_data2 = filter(row -> row.event == evt && row.peak2_left > 0, df)
            if nrow(evt_data2) > 0
                plot!(plt, evt_data2.ldrft_cm, evt_data2.peak2_left,
                      marker=:square, color=c, linestyle=:dash,
                      label="Evt $evt Peak2")
            end
        end
        savefig(plt, joinpath(outdir, "peak_lefts_combined$(suffix).png"))
    end

    println("\n✓ All plots saved to: $outdir")
end

#=============================================================================
# Main
=============================================================================#

function main(; input_csv::String, outdir::String="plots", aggregate::Bool=false)
    println("╔═══════════════════════════════════════════════════════════════╗")
    println("║          Plot Resolution Scan Results                        ║")
    println("╚═══════════════════════════════════════════════════════════════╝")

    println("\nReading: $input_csv")
    df = CSV.read(input_csv, DataFrame)

    println("Found $(nrow(df)) rows, $(length(unique(df.event))) events")
    println("Drift lengths: $(sort(unique(df.ldrft_cm))) cm")

    create_all_plots(df, outdir; aggregate=aggregate)

    # Print summary statistics
    println("\n═══════════════════════════════════════")
    println("        Summary Statistics")
    println("═══════════════════════════════════════")

    for col in [:d1_mm, :d2_mm, :Eb1_keV, :Eb2_keV, :asymmetry,
                :peak1_left, :peak1_right, :peak1_prom,
                :peak2_left, :peak2_right, :peak2_prom]
        if hasproperty(df, col)
            vals = filter(!isnan, df[!, col])
            if length(vals) > 0
                μ = mean(vals)
                σ = std(vals)
                println("$col: $(round(μ, digits=2)) ± $(round(σ, digits=2))")
            end
        end
    end

    println("═══════════════════════════════════════")
end

#=============================================================================
# Command Line Interface
=============================================================================#

if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParseSettings(description="Plot Resolution Scan Results")

    @add_arg_table! s begin
        "input"
            help = "Input CSV file (summary_all_events.csv)"
            arg_type = String
            required = true
        "--outdir", "-o"
            help = "Output directory for plots"
            arg_type = String
            default = "plots"
        "--aggregate", "-a"
            help = "Show mean ± std instead of individual events"
            action = :store_true
    end

    args = parse_args(s)

    main(; input_csv=args["input"],
           outdir=args["outdir"],
           aggregate=args["aggregate"])
end
