#!/usr/bin/env julia

"""
Blob asymmetry optimization script.

Reads reconstructed tracks (with central paths) for bb0nu and electrons,
scans sphere radius to find optimal separation based on blob energy asymmetry.

The asymmetry is defined as: A = (E1 - E2) / (E1 + E2)
where E1 >= E2 are the blob energies.

For bb0nu (two Bragg peaks): expect low asymmetry
For electrons (one Bragg peak): expect high asymmetry
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using Printf
using DataFrames
using Statistics
using Plots
using HDF5

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

"""
    BlobResult

Structure to hold blob analysis results for a single track.
"""
struct BlobResult
    Eb1::Float64      # Energy of blob 1 (higher energy) in keV
    Eb2::Float64      # Energy of blob 2 (lower energy) in keV
    asymmetry::Float64
    event_id::Int
end

"""
    compute_blob_results(results, radius)

Compute blob energies and asymmetries for all tracks at given radius.
Returns vector of BlobResult.
"""
function compute_blob_results(results, radius::Float64)
    blob_results = BlobResult[]

    for r in results
        try
            blobs = Petit.find_blob_energies(r.track, r.central_path; radius=radius)
            # Convert to keV
            Eb1_keV = blobs.Eb1 * 1e3
            Eb2_keV = blobs.Eb2 * 1e3
            push!(blob_results, BlobResult(Eb1_keV, Eb2_keV, blobs.asymmetry, r.event_id))
        catch e
            # Skip tracks that fail blob analysis
            continue
        end
    end

    return blob_results
end

"""
    compute_asymmetries(results, radius)

Compute blob energy asymmetries for all tracks at given radius.
Returns vector of asymmetry values.
"""
function compute_asymmetries(results, radius::Float64)
    blob_results = compute_blob_results(results, radius)
    return [br.asymmetry for br in blob_results]
end

"""
    compute_separation(asym_bb, asym_ele)

Compute separation power between two distributions.
Separation = |μ1 - μ2| / sqrt(σ1² + σ2²)
"""
function compute_separation(asym_bb::Vector{Float64}, asym_ele::Vector{Float64})
    if isempty(asym_bb) || isempty(asym_ele)
        return 0.0
    end

    μ_bb = mean(asym_bb)
    μ_ele = mean(asym_ele)
    σ_bb = std(asym_bb)
    σ_ele = std(asym_ele)

    combined_σ = sqrt(σ_bb^2 + σ_ele^2)
    if combined_σ < 1e-10
        return 0.0
    end

    return abs(μ_ele - μ_bb) / combined_σ
end

"""
    save_asymmetry_histogram(asymmetries, radius, label, output_dir)

Save asymmetry histogram to file.
"""
function save_asymmetry_histogram(asymmetries::Vector{Float64}, radius::Float64,
                                   label::String, output_dir::String)
    filename = joinpath(output_dir, "asymmetry_R$(Int(radius))mm_$(label).csv")

    # Save raw values
    df = DataFrame(asymmetry = asymmetries)
    Petit.CSV.write(filename, df)

    return filename
end

"""
    create_radius_plots(bb_blobs, ele_blobs, radius, output_dir)

Create plots for a specific radius:
1. Eb1 vs Eb2 for bb0nu
2. Eb1 vs Eb2 for electrons
3. Asymmetry histograms for both

Returns the combined plot.
"""
function create_radius_plots(bb_blobs::Vector{BlobResult}, ele_blobs::Vector{BlobResult},
                             radius::Float64, output_dir::String)

    # Extract data
    bb_Eb1 = [b.Eb1 for b in bb_blobs]
    bb_Eb2 = [b.Eb2 for b in bb_blobs]
    bb_asym = [b.asymmetry for b in bb_blobs]

    ele_Eb1 = [b.Eb1 for b in ele_blobs]
    ele_Eb2 = [b.Eb2 for b in ele_blobs]
    ele_asym = [b.asymmetry for b in ele_blobs]

    # Determine energy range for consistent axes
    all_E = vcat(bb_Eb1, bb_Eb2, ele_Eb1, ele_Eb2)
    E_max = min(maximum(all_E) * 1.1, 1500)

    # Plot 1: Eb1 vs Eb2 for bb0nu
    p1 = scatter(bb_Eb2, bb_Eb1,
                 xlabel="Eb2 (keV)", ylabel="Eb1 (keV)",
                 title="bb0nu R=$(Int(radius))mm (n=$(length(bb_blobs)))",
                 label="", alpha=0.6, markersize=4, color=:blue,
                 xlims=(0, E_max), ylims=(0, E_max))
    plot!(p1, [0, E_max], [0, E_max], label="Eb1=Eb2", linestyle=:dash, color=:gray)

    # Plot 2: Eb1 vs Eb2 for electrons
    p2 = scatter(ele_Eb2, ele_Eb1,
                 xlabel="Eb2 (keV)", ylabel="Eb1 (keV)",
                 title="electrons R=$(Int(radius))mm (n=$(length(ele_blobs)))",
                 label="", alpha=0.6, markersize=4, color=:red,
                 xlims=(0, E_max), ylims=(0, E_max))
    plot!(p2, [0, E_max], [0, E_max], label="Eb1=Eb2", linestyle=:dash, color=:gray)

    # Plot 3: Asymmetry histograms
    p3 = histogram(bb_asym,
                   bins=range(0, 1, length=25),
                   alpha=0.6, label="bb0nu",
                   xlabel="Asymmetry", ylabel="Counts",
                   title="Asymmetry R=$(Int(radius))mm",
                   color=:blue)
    histogram!(p3, ele_asym,
               bins=range(0, 1, length=25),
               alpha=0.6, label="electrons",
               color=:red)

    # Compute separation for annotation
    sep = compute_separation(bb_asym, ele_asym)
    annotate!(p3, [(0.7, maximum(ylims(p3)[2])*0.9,
                    text("Sep: $(round(sep, digits=2))σ", 10, :left))])

    # Combined plot for this radius
    p_combined = plot(p1, p2, p3, layout=(1, 3), size=(1400, 400))

    # Save
    plot_file = joinpath(output_dir, "blob_analysis_R$(Int(radius))mm.png")
    savefig(p_combined, plot_file)

    return p_combined
end

"""
    run_optimization(bb_dir, bb_pattern, ele_dir, ele_pattern; kwargs...)

Run the asymmetry optimization scan.

# Arguments
- `r_start`: Initial radius in mm (default: 10.0)
- `r_end`: Final radius in mm (default: 30.0)
- `r_step`: Step size in mm (default: 2.0)
- `output_dir`: Directory for output files
"""
function run_optimization(bb_dir::String, bb_pattern::String,
                          ele_dir::String, ele_pattern::String;
                          r_start::Float64=10.0,
                          r_end::Float64=30.0,
                          r_step::Float64=2.0,
                          output_dir::String="asymmetry_results")

    # Create output directory if it doesn't exist
    if !isdir(output_dir)
        mkpath(output_dir)
        println("Created output directory: $output_dir")
    end

    println("="^70)
    println("BLOB ASYMMETRY OPTIMIZATION")
    println("="^70)
    println("Output directory: $output_dir")

    # Load bb0nu results
    println("\nLoading bb0nu data...")
    bb_results, bb_metadata = Petit.chain_reco_results(bb_dir, bb_pattern)
    println("  Loaded $(length(bb_results)) bb0nu tracks")

    # Load electron results
    println("\nLoading electron data...")
    ele_results, ele_metadata = Petit.chain_reco_results(ele_dir, ele_pattern)
    println("  Loaded $(length(ele_results)) electron tracks")

    if isempty(bb_results) || isempty(ele_results)
        error("No data loaded. Check file paths and patterns.")
    end

    # Storage for results
    radii = Float64[]
    separations = Float64[]
    bb_means = Float64[]
    bb_stds = Float64[]
    ele_means = Float64[]
    ele_stds = Float64[]

    # Radius scan
    println("\n" * "-"^70)
    println("RADIUS SCAN: R = $r_start to $r_end mm (step = $r_step mm)")
    println("-"^70)
    @printf("%-10s %10s %10s %10s %10s %12s\n",
            "R (mm)", "bb μ", "bb σ", "ele μ", "ele σ", "Separation")
    println("-"^70)

    best_separation = 0.0
    best_radius = r_start

    r = r_start
    while r <= r_end
        # Compute blob results (energies and asymmetries)
        bb_blobs = compute_blob_results(bb_results, r)
        ele_blobs = compute_blob_results(ele_results, r)

        if isempty(bb_blobs) || isempty(ele_blobs)
            println("Warning: No valid blob results at R = $r mm")
            r += r_step
            continue
        end

        # Extract asymmetries
        asym_bb = [b.asymmetry for b in bb_blobs]
        asym_ele = [b.asymmetry for b in ele_blobs]

        # Compute statistics
        μ_bb = mean(asym_bb)
        σ_bb = std(asym_bb)
        μ_ele = mean(asym_ele)
        σ_ele = std(asym_ele)
        separation = compute_separation(asym_bb, asym_ele)

        @printf("%-10.1f %10.3f %10.3f %10.3f %10.3f %12.3f\n",
                r, μ_bb, σ_bb, μ_ele, σ_ele, separation)

        # Store results
        push!(radii, r)
        push!(separations, separation)
        push!(bb_means, μ_bb)
        push!(bb_stds, σ_bb)
        push!(ele_means, μ_ele)
        push!(ele_stds, σ_ele)

        # Save histograms (CSV)
        save_asymmetry_histogram(asym_bb, r, "bb0nu", output_dir)
        save_asymmetry_histogram(asym_ele, r, "electrons", output_dir)

        # Create and save plots for this radius
        create_radius_plots(bb_blobs, ele_blobs, r, output_dir)
        println("  Saved plots for R = $r mm")

        # Update best
        if separation > best_separation
            best_separation = separation
            best_radius = r
        end

        r += r_step
    end

    # Summary
    println("\n" * "="^70)
    println("OPTIMIZATION RESULTS")
    println("="^70)
    println("\nBest radius: $(best_radius) mm")
    println("Best separation: $(round(best_separation, digits=3)) σ")

    # Find stats at best radius
    best_idx = findfirst(x -> x == best_radius, radii)
    if !isnothing(best_idx)
        println("\nAt R = $(best_radius) mm:")
        println("  bb0nu:     μ = $(round(bb_means[best_idx], digits=3)), σ = $(round(bb_stds[best_idx], digits=3))")
        println("  electrons: μ = $(round(ele_means[best_idx], digits=3)), σ = $(round(ele_stds[best_idx], digits=3))")
    end

    # Create summary plots
    println("\nGenerating summary plots...")

    # Plot 1: Separation vs Radius
    p1 = plot(radii, separations,
              xlabel="Sphere Radius (mm)",
              ylabel="Separation (σ)",
              title="Separation Power vs Radius",
              marker=:circle,
              linewidth=2,
              legend=false,
              grid=true)
    vline!(p1, [best_radius], linestyle=:dash, color=:red, label="Best R=$(Int(best_radius))mm")

    # Plot 2: Mean asymmetry vs Radius
    p2 = plot(radii, bb_means,
              ribbon=bb_stds,
              fillalpha=0.3,
              label="bb0nu",
              xlabel="Sphere Radius (mm)",
              ylabel="Mean Asymmetry",
              title="Asymmetry vs Radius",
              linewidth=2,
              color=:blue)
    plot!(p2, radii, ele_means,
          ribbon=ele_stds,
          fillalpha=0.3,
          label="electrons",
          linewidth=2,
          color=:red)
    vline!(p2, [best_radius], linestyle=:dash, color=:black, label="Best R")

    # Plot 3: Histograms at best radius
    asym_bb_best = compute_asymmetries(bb_results, best_radius)
    asym_ele_best = compute_asymmetries(ele_results, best_radius)

    p3 = histogram(asym_bb_best,
                   bins=range(0, 1, length=30),
                   alpha=0.6,
                   label="bb0nu (n=$(length(asym_bb_best)))",
                   xlabel="Asymmetry",
                   ylabel="Counts",
                   title="Asymmetry at Best R = $(Int(best_radius)) mm",
                   color=:blue)
    histogram!(p3, asym_ele_best,
               bins=range(0, 1, length=30),
               alpha=0.6,
               label="electrons (n=$(length(asym_ele_best)))",
               color=:red)
    annotate!(p3, [(0.7, maximum(ylims(p3)[2])*0.9,
                    text("Sep: $(round(best_separation, digits=2))σ", 10, :left))])

    # Combined summary plot
    p_summary = plot(p1, p2, p3, layout=(1, 3), size=(1500, 400))

    summary_file = joinpath(output_dir, "blob_asymmetry_summary.png")
    savefig(p_summary, summary_file)
    println("Saved: $summary_file")

    # Global asymmetry vs radius plot (standalone)
    p_asym_vs_r = plot(radii, bb_means,
                       yerror=bb_stds,
                       label="bb0nu",
                       xlabel="Sphere Radius (mm)",
                       ylabel="Asymmetry",
                       title="Blob Energy Asymmetry vs Sphere Radius",
                       marker=:circle,
                       linewidth=2,
                       color=:blue,
                       grid=true,
                       size=(800, 500))
    plot!(p_asym_vs_r, radii, ele_means,
          yerror=ele_stds,
          label="electrons",
          marker=:square,
          linewidth=2,
          color=:red)
    vline!(p_asym_vs_r, [best_radius], linestyle=:dash, color=:black,
           label="Best R=$(Int(best_radius))mm")

    asym_vs_r_file = joinpath(output_dir, "asymmetry_vs_radius.png")
    savefig(p_asym_vs_r, asym_vs_r_file)
    println("Saved: $asym_vs_r_file")

    display(p_summary)

    # Save results table
    results_df = DataFrame(
        radius_mm = radii,
        separation = separations,
        bb_mean = bb_means,
        bb_std = bb_stds,
        ele_mean = ele_means,
        ele_std = ele_stds
    )
    results_file = joinpath(output_dir, "blob_asymmetry_results.csv")
    Petit.CSV.write(results_file, results_df)
    println("Saved: $results_file")

    return (best_radius=best_radius, best_separation=best_separation,
            radii=radii, separations=separations, results_df=results_df)
end

"""
    main()

Main function for command-line usage.
"""
function main()
    # Default paths - adjust as needed
    bb_dir = joinpath(ENV["DATA"], "HD5t/itaca/bb0nu")
    bb_pattern = "bb0nu_test_reco_th_*.h5"

    ele_dir = joinpath(ENV["DATA"], "HD5t/itaca/xe137")
    ele_pattern = "electrons_test_reco_th_*.h5"

    # Default parameters
    output_dir = "asymmetry_results"
    r_start = 10.0
    r_end = 30.0
    r_step = 2.0

    # Parse command line arguments
    positional_args = String[]
    for arg in ARGS
        if startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) != 2
                println("Warning: Ignoring malformed argument: $arg")
                continue
            end
            key, value = parts

            try
                if key == "rstart"
                    r_start = parse(Float64, value)
                elseif key == "rend"
                    r_end = parse(Float64, value)
                elseif key == "rstep"
                    r_step = parse(Float64, value)
                elseif key == "outdir"
                    output_dir = String(value)
                else
                    println("Warning: Unknown argument: --$key")
                end
            catch e
                println("Error parsing argument --$key=$value: $e")
                exit(1)
            end
        else
            push!(positional_args, arg)
        end
    end

    # Use positional arguments if provided
    if length(positional_args) >= 4
        bb_dir = positional_args[1]
        bb_pattern = positional_args[2]
        ele_dir = positional_args[3]
        ele_pattern = positional_args[4]
    end

    if length(ARGS) < 4
        println("Usage: julia blob_asymmetry_optimization.jl <bb_dir> <bb_pattern> <ele_dir> <ele_pattern> [options]")
        println("\nOptional arguments:")
        println("  --rstart=X    Initial radius in mm (default: 10.0)")
        println("  --rend=X      Final radius in mm (default: 30.0)")
        println("  --rstep=X     Step size in mm (default: 2.0)")
        println("  --outdir=X    Output directory for plots (default: asymmetry_results)")
        println("\nUsing default paths:")
        println("  bb0nu:     $bb_dir / $bb_pattern")
        println("  electrons: $ele_dir / $ele_pattern")
    end

    # Run optimization
    run_optimization(bb_dir, bb_pattern, ele_dir, ele_pattern;
                     r_start=r_start,
                     r_end=r_end,
                     r_step=r_step,
                     output_dir=output_dir)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
