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
using Distributions
using Statistics
using Printf
using Plots

"""
    parse_commandline()

Parse command line arguments for figure of merit analysis.
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input", "-i"
            help = "Input CSV file containing blob analysis data"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Output CSV file for efficiency data (default: input_fm.csv)"
            arg_type = String
            default = ""
        "--min-cut"
            help = "Minimum eblob2 cut value in keV"
            arg_type = Float64
            default = 500.0
        "--max-cut"
            help = "Maximum eblob2 cut value in keV"
            arg_type = Float64
            default = 700.0
        "--step"
            help = "Step size for eblob2 cut variation in keV"
            arg_type = Float64
            default = 20.0
        "--plot"
            help = "Generate efficiency plot"
            action = :store_true
        "--plot-output"
            help = "Output file for efficiency plot"
            arg_type = String
            default = ""
        "--verbose", "-v"
            help = "Verbose output"
            action = :store_true
        "--confidence-level"
            help = "Confidence level for error calculation (e.g., 0.68 for 1-sigma)"
            arg_type = Float64
            default = 0.6827  # 1-sigma confidence level
    end

    return parse_args(s)
end

"""
    binomial_error(k::Int, n::Int; confidence_level=0.6827)

Calculate the binomial proportion error using the Clopper-Pearson method.

# Arguments
- `k`: Number of successes (events passing cut)
- `n`: Total number of trials (total events)
- `confidence_level`: Confidence level for the interval (default: 0.6827 for 1-sigma)

# Returns
- Tuple of (efficiency, lower_error, upper_error)
"""
function binomial_error(k::Int, n::Int; confidence_level=0.6827)
    if n == 0
        return (0.0, 0.0, 0.0)
    end

    # Point estimate
    p = k / n

    # For the binomial distribution, we use the Clopper-Pearson interval
    # This is the exact confidence interval for a binomial proportion
    alpha = 1.0 - confidence_level
    alpha_2 = alpha / 2.0

    # Lower bound
    if k == 0
        lower = 0.0
    else
        # Use the quantile of the Beta distribution
        # Beta(k, n - k + 1) for lower bound
        lower = quantile(Beta(k, n - k + 1), alpha_2)
    end

    # Upper bound
    if k == n
        upper = 1.0
    else
        # Beta(k + 1, n - k) for upper bound
        upper = quantile(Beta(k + 1, n - k), 1.0 - alpha_2)
    end

    # Return efficiency and symmetric error (average of upper and lower)
    lower_error = p - lower
    upper_error = upper - p

    # For simplicity, return the average error
    avg_error = (lower_error + upper_error) / 2.0

    return (p, avg_error)
end

"""
    simple_binomial_error(k::Int, n::Int)

Calculate the binomial proportion error using simple square root approximation.
This is valid when np and n(1-p) are both > 5.

# Arguments
- `k`: Number of successes (events passing cut)
- `n`: Total number of trials (total events)

# Returns
- Tuple of (efficiency, error)
"""
function simple_binomial_error(k::Int, n::Int)
    if n == 0
        return (0.0, 0.0)
    end

    p = k / n

    # Standard error for binomial proportion
    # σ = sqrt(p(1-p)/n)
    error = sqrt(p * (1 - p) / n)

    return (p, error)
end

"""
    calculate_efficiency_curve(df::DataFrame, min_cut::Float64, max_cut::Float64, step::Float64;
                              verbose=false, use_simple_error=true)

Calculate efficiency as a function of eblob2 cut value.

# Arguments
- `df`: DataFrame containing the data
- `min_cut`: Minimum cut value
- `max_cut`: Maximum cut value
- `step`: Step size for cut variation
- `verbose`: Print detailed information
- `use_simple_error`: Use simple error calculation (true) or Clopper-Pearson (false)

# Returns
- DataFrame with columns: eblob2_cut, efficiency, error
"""
function calculate_efficiency_curve(df::DataFrame, min_cut::Float64, max_cut::Float64,
                                   step::Float64; verbose=false, use_simple_error=true,
                                   confidence_level=0.6827)

    # Check if eblob2 column exists
    if !("eblob2" in names(df))
        error("Column 'eblob2' not found in dataframe. Available columns: $(names(df))")
    end

    # Total number of events
    n_total = nrow(df)

    if verbose
        println("Total events in dataset: $n_total")
        println("eblob2 range in data: [$(minimum(df.eblob2)), $(maximum(df.eblob2))] keV")
    end

    # Initialize results arrays
    cuts = Float64[]
    efficiencies = Float64[]
    errors = Float64[]

    # Iterate over cut values
    cut_value = min_cut
    while cut_value <= max_cut
        # Count events passing the cut
        n_pass = sum(df.eblob2 .>= cut_value)

        # Calculate efficiency and error
        if use_simple_error
            eff, err = simple_binomial_error(n_pass, n_total)
        else
            eff, err = binomial_error(n_pass, n_total; confidence_level=confidence_level)
        end

        push!(cuts, cut_value)
        push!(efficiencies, eff)
        push!(errors, err)

        if verbose
            @printf("Cut: %.1f keV | Pass: %d/%d | Eff: %.4f ± %.4f (%.2f%%)\n",
                    cut_value, n_pass, n_total, eff, err, eff * 100)
        end

        cut_value += step
    end

    # Create output dataframe
    results_df = DataFrame(
        eblob2 = cuts,
        eff = round.(efficiencies, digits=4),
        err = round.(errors, digits=4)
    )

    return results_df
end

"""
    create_efficiency_plot(results_df::DataFrame, title::String="")

Create a plot of efficiency vs eblob2 cut with error bars.
"""
function create_efficiency_plot(results_df::DataFrame, title::String="")
    # Create the plot with error bars
    p = scatter(results_df.eblob2, results_df.eff,
                yerror = results_df.err,
                xlabel = "Blob 2 Energy Cut (keV)",
                ylabel = "Efficiency",
                title = isempty(title) ? "Efficiency vs Energy Cut" : title,
                label = "Data",
                markersize = 5,
                markercolor = :blue,
                markerstrokewidth = 1,
                markerstrokecolor = :black,
                linecolor = :blue,
                linewidth = 1,
                ylims = (0, maximum(results_df.eff) * 1.1),
                grid = true,
                gridstyle = :dot,
                gridalpha = 0.3,
                legend = :topright)

    # Add a smooth line through the points
    plot!(p, results_df.eblob2, results_df.eff,
          line = :solid,
          linewidth = 2,
          label = nothing,
          color = :blue,
          alpha = 0.5)

    # Add percentage labels on secondary y-axis
    plot!(p, yticks = (0:0.1:1.0, ["$(Int(y*100))%" for y in 0:0.1:1.0]))

    return p
end

"""
    main()

Main function for figure of merit analysis.
"""
function main()
    # Parse command line arguments
    args = parse_commandline()

    # Set default output filename if not provided
    if isempty(args["output"])
        base_name = splitext(basename(args["input"]))[1]
        args["output"] = "$(base_name)_fm.csv"
    end

    if isempty(args["plot-output"]) && args["plot"]
        base_name = splitext(basename(args["input"]))[1]
        args["plot-output"] = "$(base_name)_fm.png"
    end

    println("="^70)
    println("Figure of Merit Analysis")
    println("="^70)
    println("Input file:  $(args["input"])")
    println("Output file: $(args["output"])")
    println("Cut range:   $(args["min-cut"]) to $(args["max-cut"]) keV")
    println("Step size:   $(args["step"]) keV")
    println("="^70)

    # Check if input file exists
    if !isfile(args["input"])
        error("Input file not found: $(args["input"])")
    end

    # Read input CSV file
    println("\nReading data from $(args["input"])...")
    df = CSV.read(args["input"], DataFrame)
    println("Loaded $(nrow(df)) events")

    # Print data summary
    if "eblob2" in names(df)
        println("\neblob2 statistics:")
        println("  Mean:   $(round(mean(df.eblob2), digits=1)) keV")
        println("  Median: $(round(median(df.eblob2), digits=1)) keV")
        println("  Std:    $(round(std(df.eblob2), digits=1)) keV")
        println("  Min:    $(round(minimum(df.eblob2), digits=1)) keV")
        println("  Max:    $(round(maximum(df.eblob2), digits=1)) keV")
    end

    # Calculate efficiency curve
    println("\nCalculating efficiency curve...")
    results_df = calculate_efficiency_curve(
        df,
        args["min-cut"],
        args["max-cut"],
        args["step"],
        verbose = args["verbose"],
        use_simple_error = true,
        confidence_level = args["confidence-level"]
    )

    # Save results to CSV
    CSV.write(args["output"], results_df)
    println("\n✓ Efficiency data saved to: $(args["output"])")

    # Print summary table
    println("\n" * "="^70)
    println("Efficiency Summary")
    println("="^70)
    println("Cut (keV) | Efficiency (%) | Error (%)")
    println("-"^40)
    for row in eachrow(results_df)
        @printf("%8.1f  | %13.2f  | %9.2f\n",
                row.eblob2, row.eff * 100, row.err * 100)
    end
    println("="^70)

    # Find and report optimal points
    println("\nKey Points:")
    println("-"^40)

    # Find cut for specific efficiencies
    target_effs = [0.90, 0.95, 0.99]
    for target_eff in target_effs
        # Find the cut value closest to target efficiency
        idx = argmin(abs.(results_df.eff .- target_eff))
        if abs(results_df.eff[idx] - target_eff) < 0.05  # Within 5% of target
            @printf("For ~%.0f%% efficiency: cut = %.1f keV (actual: %.2f%%)\n",
                    target_eff * 100,
                    results_df.eblob2[idx],
                    results_df.eff[idx] * 100)
        end
    end

    # Report efficiency at specific cut values
    println("\nEfficiency at standard cuts:")
    for standard_cut in [550.0, 600.0, 650.0]
        if standard_cut >= args["min-cut"] && standard_cut <= args["max-cut"]
            idx = argmin(abs.(results_df.eblob2 .- standard_cut))
            if abs(results_df.eblob2[idx] - standard_cut) < args["step"]
                @printf("At %.0f keV: %.2f%% ± %.2f%%\n",
                        standard_cut,
                        results_df.eff[idx] * 100,
                        results_df.err[idx] * 100)
            end
        end
    end

    # Create plot if requested
    if args["plot"]
        println("\nGenerating efficiency plot...")

        # Extract base name for title
        title = "Efficiency Curve - $(splitext(basename(args["input"]))[1])"
        p = create_efficiency_plot(results_df, title)

        # Save plot
        savefig(p, args["plot-output"])
        println("✓ Plot saved to: $(args["plot-output"])")
    end

    println("\n✓ Analysis complete!")
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end