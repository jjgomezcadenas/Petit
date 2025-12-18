#!/usr/bin/env julia

"""
Test script for itaca_aux.jl functions
Tests the two Gaussian filtering implementations
"""

using Pkg
Pkg.activate(joinpath(ENV["PROJECTS"], "Petit"))

using Petit
using DataFrames
using Statistics
using Plots
using Images
using ImageFiltering

# Include the auxiliary functions
include("itaca_aux.jl")

println("="^80)
println("Testing Gaussian Filter Implementations")
println("="^80)

# Load data
println("\nLoading data...")
xfile = "/Users/jjgomezcadenas/Data/HD5t/lxe/lxe_0nubb.next.h5"
dfs = get_dataset_dfs(xfile)
hitsdf = dfs["hits"]
println("Loaded $(nrow(hitsdf)) hits")

# Get a test event
nevent = 2
println("\nExtracting event $nevent...")
evtdf = get_event(hitsdf, nevent)
println("Event has $(nrow(evtdf)) hits")

# Run the test
println("\n" * "="^80)
println("Running Gaussian filter comparison...")
println("="^80)

result = test_gaussian_filters(evtdf; sigma_mm=0.6, nbins=100)

# Display results
println("\n" * "="^80)
println("FINAL RESULTS:")
println("="^80)
println("\nStatistics:")
println("  Manual method:")
println("    Total intensity: ", result.stats.manual_total)
println("    Max intensity:   ", result.stats.manual_max)
println("  Images.jl method:")
println("    Total intensity: ", result.stats.efficient_total)
println("    Max intensity:   ", result.stats.efficient_max)
println("\n  Ratios:")
println("    Total (manual/efficient):  ", result.stats.ratio_total)
println("    Max (manual/efficient):    ", result.stats.ratio_max)
println("    Correlation coefficient:   ", result.stats.correlation)

# Save the comparison plot
println("\nSaving comparison plot...")
savefig(result.plot, "gaussian_filter_comparison.png")
println("Plot saved to: gaussian_filter_comparison.png")

# Check if they are compatible (within reasonable tolerance)
println("\n" * "="^80)
println("COMPATIBILITY CHECK:")
println("="^80)

ratio_tol = 0.1  # 10% tolerance
corr_tol = 0.95  # Correlation should be > 0.95

ratio_ok = abs(result.stats.ratio_total - 1.0) < ratio_tol
corr_ok = result.stats.correlation > corr_tol

println("  Ratio test (should be close to 1.0):       ", ratio_ok ? "✓ PASS" : "✗ FAIL")
println("  Correlation test (should be > $corr_tol):  ", corr_ok ? "✓ PASS" : "✗ FAIL")

if ratio_ok && corr_ok
    println("\n✓ Both implementations are COMPATIBLE!")
else
    println("\n✗ Implementations differ significantly!")
    println("\nPossible issues:")
    if !ratio_ok
        println("  - Total intensity differs by more than $(ratio_tol*100)%")
        println("    This suggests different normalization or spreading behavior")
    end
    if !corr_ok
        println("  - Correlation is below $corr_tol")
        println("    This suggests different spatial distributions")
    end
end

println("\n" * "="^80)
println("Test complete!")
println("="^80)
