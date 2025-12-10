#!/usr/bin/env julia

"""
Simple test of KernelDensity package.
"""

const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)

using KernelDensity
using Plots
using Random

# Generate samples from N(0,1)
Random.seed!(42)
data = randn(1000)

# Compute Gaussian KDE with automatic bandwidth
kde_obj = kde(data)

# Evaluate KDE at a grid of points
xs = range(-4, 4, length=200)
ys = pdf(kde_obj, xs)

# Plot
p = plot(xs, ys,
         label="KDE",
         xlabel="x",
         ylabel="Density",
         title="Kernel Density Estimation\n(1000 samples from N(0,1))",
         linewidth=2,
         color=:blue)

# Add theoretical N(0,1) for comparison
using Distributions
theoretical = pdf.(Normal(0, 1), xs)
plot!(p, xs, theoretical,
      label="True N(0,1)",
      linestyle=:dash,
      linewidth=2,
      color=:red)

# Add histogram of data
histogram!(p, data,
           normalize=:pdf,
           alpha=0.3,
           label="Histogram",
           color=:gray)

# Save plot
output_file = joinpath(pdir, "scripts", "kde_test.png")
savefig(p, output_file)
println("Plot saved to: $output_file")

# Print KDE info
println("\nKDE x range: $(kde_obj.x[1]) to $(kde_obj.x[end])")
println("KDE grid points: $(length(kde_obj.x))")
