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
using Statistics
using Printf
using Plots
using StatsPlots

"""
    parse_commandline()

Parse command line arguments for selection statistics script.
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input", "-i"
            help = "Input CSV file(s) - can specify multiple files"
            arg_type = String
            nargs = '+'
            required = true
        "--output-stats", "-o"
            help = "Output file for statistics summary (optional)"
            arg_type = String
            default = ""
        "--plot-output", "-p"
            help = "Output file for comparison plots (optional)"
            arg_type = String
            default = ""
        "--verbose", "-v"
            help = "Verbose output with detailed statistics"
            action = :store_true
        "--compare-mode"
            help = "Compare distributions between selected and non-selected"
            action = :store_true
    end

    return parse_args(s)
end

"""
    compute_statistics(df::DataFrame, column::Symbol)

Compute mean, std, min, max, and median for a given column.
"""
function compute_statistics(df::DataFrame, column::Symbol)
    if nrow(df) == 0
        return (mean=NaN, std=NaN, min=NaN, max=NaN, median=NaN, count=0)
    end

    values = df[!, column]
    return (
        mean = mean(values),
        std = std(values),
        min = minimum(values),
        max = maximum(values),
        median = median(values),
        count = length(values)
    )
end

"""
    analyze_selection(df::DataFrame; verbose=false)

Analyze statistics for selected and non-selected events.
"""
function analyze_selection(df::DataFrame; verbose=false)
    # Check if 'select' column exists
    if !("select" in names(df))
        error("Column 'select' not found in dataframe. Available columns: $(names(df))")
    end

    # Split data by selection
    df_selected = filter(row -> row.select == true, df)
    df_not_selected = filter(row -> row.select == false, df)

    n_total = nrow(df)
    n_selected = nrow(df_selected)
    n_not_selected = nrow(df_not_selected)

    # Calculate fractions
    frac_selected = n_selected / n_total
    frac_not_selected = n_not_selected / n_total

    # Compute statistics for key columns
    columns_to_analyze = [:etrk, :trkl, :eblob1, :eblob2, :nblob1, :nblob2, :confidence]

    stats_selected = Dict{Symbol, Any}()
    stats_not_selected = Dict{Symbol, Any}()
    stats_total = Dict{Symbol, Any}()

    for col in columns_to_analyze
        if col in propertynames(df)
            stats_selected[col] = compute_statistics(df_selected, col)
            stats_not_selected[col] = compute_statistics(df_not_selected, col)
            stats_total[col] = compute_statistics(df, col)
        end
    end

    # Create results dictionary
    results = Dict(
        :n_total => n_total,
        :n_selected => n_selected,
        :n_not_selected => n_not_selected,
        :frac_selected => frac_selected,
        :frac_not_selected => frac_not_selected,
        :stats_selected => stats_selected,
        :stats_not_selected => stats_not_selected,
        :stats_total => stats_total,
        :df_selected => df_selected,
        :df_not_selected => df_not_selected
    )

    return results
end

"""
    print_statistics(results::Dict, filename::String=""; verbose=false)

Print formatted statistics from analysis results.
"""
function print_statistics(results::Dict, filename::String=""; verbose=false)
    println("\n" * "="^70)
    if !isempty(filename)
        println("File: $filename")
        println("="^70)
    end

    println("Selection Summary:")
    println("-"^35)
    println("Total events:        $(results[:n_total])")
    println("Selected (true):     $(results[:n_selected]) ($(round(100*results[:frac_selected], digits=1))%)")
    println("Not selected (false): $(results[:n_not_selected]) ($(round(100*results[:frac_not_selected], digits=1))%)")

    # Print comparative statistics
    println("\n" * "="^70)
    println("Comparative Statistics (mean ± std):")
    println("="^70)

    # Header
    @printf("%-15s | %-20s | %-20s | %-20s\n", "Variable", "All Events", "Selected", "Not Selected")
    println("-"^80)

    for col in [:etrk, :trkl, :eblob1, :eblob2, :confidence]
        if haskey(results[:stats_total], col)
            stats_all = results[:stats_total][col]
            stats_sel = results[:stats_selected][col]
            stats_not = results[:stats_not_selected][col]

            col_name = string(col)
            if col == :etrk
                col_name = "Track E (keV)"
            elseif col == :trkl
                col_name = "Track L (mm)"
            elseif col == :eblob1
                col_name = "Blob1 E (keV)"
            elseif col == :eblob2
                col_name = "Blob2 E (keV)"
            elseif col == :confidence
                col_name = "Confidence"
            end

            @printf("%-15s | %7.1f ± %6.1f | %7.1f ± %6.1f | %7.1f ± %6.1f\n",
                    col_name,
                    stats_all.mean, stats_all.std,
                    stats_sel.mean, stats_sel.std,
                    stats_not.mean, stats_not.std)
        end
    end

    if verbose
        println("\n" * "="^70)
        println("Detailed Statistics:")
        println("="^70)

        for (label, stats_dict) in [("SELECTED EVENTS", results[:stats_selected]),
                                     ("NOT SELECTED EVENTS", results[:stats_not_selected])]
            println("\n$label:")
            println("-"^35)

            for col in [:etrk, :trkl, :eblob1, :eblob2, :nblob1, :nblob2, :confidence]
                if haskey(stats_dict, col)
                    stats = stats_dict[col]
                    println("\n  $col:")
                    println("    Count:  $(stats.count)")
                    println("    Mean:   $(round(stats.mean, digits=2))")
                    println("    Std:    $(round(stats.std, digits=2))")
                    println("    Min:    $(round(stats.min, digits=2))")
                    println("    Median: $(round(stats.median, digits=2))")
                    println("    Max:    $(round(stats.max, digits=2))")
                end
            end
        end
    end

    println("\n" * "="^70)
end

"""
    create_comparison_plots(results::Dict, title_suffix::String="")

Create comparison plots between selected and non-selected events.
"""
function create_comparison_plots(results::Dict, title_suffix::String="")
    df_sel = results[:df_selected]
    df_not = results[:df_not_selected]

    # Create subplots comparing distributions
    p1 = begin
        histogram(df_sel.etrk, bins=30, alpha=0.5, label="Selected", color=:green,
                 normalize=:probability)
        histogram!(df_not.etrk, bins=30, alpha=0.5, label="Not Selected", color=:red,
                  normalize=:probability)
        xlabel!("Track Energy (keV)")
        ylabel!("Probability")
        title!("Track Energy Distribution")
    end

    p2 = begin
        histogram(df_sel.trkl, bins=30, alpha=0.5, label="Selected", color=:green,
                 normalize=:probability)
        histogram!(df_not.trkl, bins=30, alpha=0.5, label="Not Selected", color=:red,
                  normalize=:probability)
        xlabel!("Track Length (mm)")
        ylabel!("Probability")
        title!("Track Length Distribution")
    end

    p3 = begin
        histogram(df_sel.eblob2, bins=30, alpha=0.5, label="Selected", color=:green,
                 normalize=:probability)
        histogram!(df_not.eblob2, bins=30, alpha=0.5, label="Not Selected", color=:red,
                  normalize=:probability)
        xlabel!("Blob 2 Energy (keV)")
        ylabel!("Probability")
        title!("Blob 2 Energy Distribution")
    end

    p4 = begin
        scatter(df_sel.eblob1, df_sel.eblob2, alpha=0.5, label="Selected",
               color=:green, markersize=3)
        scatter!(df_not.eblob1, df_not.eblob2, alpha=0.5, label="Not Selected",
                color=:red, markersize=3)
        xlabel!("Blob 1 Energy (keV)")
        ylabel!("Blob 2 Energy (keV)")
        title!("Blob Energy Correlation")
    end

    # Combine plots
    plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 900),
         plot_title="Selection Comparison" * (isempty(title_suffix) ? "" : " - $title_suffix"))
end

"""
    save_statistics_summary(all_results::Vector{Tuple{String, Dict}}, output_file::String)

Save a summary of statistics from multiple files to a CSV file.
"""
function save_statistics_summary(all_results::Vector{Tuple{String, Dict}}, output_file::String)
    # Create a summary dataframe
    rows = []

    for (filename, results) in all_results
        row = Dict(
            :file => basename(filename),
            :n_total => results[:n_total],
            :n_selected => results[:n_selected],
            :n_not_selected => results[:n_not_selected],
            :frac_selected => round(results[:frac_selected], digits=4),
            :frac_not_selected => round(results[:frac_not_selected], digits=4)
        )

        # Add mean values for key variables
        for col in [:etrk, :trkl, :eblob1, :eblob2, :confidence]
            if haskey(results[:stats_total], col)
                row[Symbol("mean_$(col)_all")] = round(results[:stats_total][col].mean, digits=2)
                row[Symbol("mean_$(col)_selected")] = round(results[:stats_selected][col].mean, digits=2)
                row[Symbol("mean_$(col)_not_selected")] = round(results[:stats_not_selected][col].mean, digits=2)
            end
        end

        push!(rows, row)
    end

    summary_df = DataFrame(rows)
    CSV.write(output_file, summary_df)
    println("\nSaved statistics summary to: $output_file")
end

"""
    main()

Main function for selection statistics analysis.
"""
function main()
    # Parse command line arguments
    args = parse_commandline()

    println("="^70)
    println("Selection Statistics Analysis")
    println("="^70)

    # Store results from all files
    all_results = Tuple{String, Dict}[]

    # Process each input file
    for input_file in args["input"]
        if !isfile(input_file)
            println("WARNING: File not found: $input_file")
            continue
        end

        println("\nProcessing: $input_file")

        # Read CSV file
        df = CSV.read(input_file, DataFrame)
        println("Loaded $(nrow(df)) rows")

        # Analyze selection
        try
            results = analyze_selection(df; verbose=args["verbose"])
            push!(all_results, (input_file, results))

            # Print statistics
            print_statistics(results, basename(input_file); verbose=args["verbose"])

            # Create comparison plots if requested
            if args["compare-mode"] && !isempty(args["plot-output"])
                p = create_comparison_plots(results, basename(input_file))

                # Determine output filename for this file's plot
                base_name = splitext(args["plot-output"])[1]
                ext = splitext(args["plot-output"])[2]
                if isempty(ext)
                    ext = ".png"
                end

                if length(args["input"]) > 1
                    # Multiple files: add suffix
                    plot_file = "$(base_name)_$(splitext(basename(input_file))[1])$(ext)"
                else
                    plot_file = "$(base_name)$(ext)"
                end

                savefig(p, plot_file)
                println("Saved comparison plot to: $plot_file")
            end

        catch e
            println("ERROR processing $input_file: $e")
        end
    end

    # If multiple files were processed, show combined statistics
    if length(all_results) > 1
        println("\n" * "="^70)
        println("COMBINED SUMMARY")
        println("="^70)

        total_events = sum(r[2][:n_total] for r in all_results)
        total_selected = sum(r[2][:n_selected] for r in all_results)
        total_not_selected = sum(r[2][:n_not_selected] for r in all_results)

        println("Total events across all files: $total_events")
        println("Total selected:     $total_selected ($(round(100*total_selected/total_events, digits=1))%)")
        println("Total not selected: $total_not_selected ($(round(100*total_not_selected/total_events, digits=1))%)")

        # File-by-file summary
        println("\nFile-by-file selection rates:")
        println("-"^50)
        for (filename, results) in all_results
            @printf("%-30s: %5d events, %5.1f%% selected\n",
                    basename(filename)[1:min(30, length(basename(filename)))],
                    results[:n_total],
                    100*results[:frac_selected])
        end
    end

    # Save summary to file if requested
    if !isempty(args["output-stats"]) && !isempty(all_results)
        save_statistics_summary(all_results, args["output-stats"])
    end

    println("\n✓ Analysis complete!")
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end