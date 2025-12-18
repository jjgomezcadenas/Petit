#!/usr/bin/env julia

"""
Signal-to-noise figure of merit script.

Computes the figure of merit for bb0nu vs electron separation using
blob energy cuts. The FOM is defined as:

    FOM = ε_signal / √ε_background

where ε is the efficiency for events passing the Eb2 cut.
Eb1 is the higher energy blob, Eb2 is the lower energy blob.
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
    BlobEnergies

Structure to hold blob energies for a single track.
"""
struct BlobEnergies
    Eb1::Float64      # Energy of blob 1 (higher energy) in keV
    Eb2::Float64      # Energy of blob 2 (lower energy) in keV
    event_id::Int
end

"""
    compute_blob_energies(results, radius)

Compute blob energies for all tracks at given radius.
Returns vector of BlobEnergies.
"""
function compute_blob_energies(results, radius::Float64)
    blob_energies = BlobEnergies[]

    for r in results
        try
            blobs = Petit.find_blob_energies(r.track, r.central_path; radius=radius)
            # Convert to keV - Eb1 is higher energy, Eb2 is lower energy
            Eb1_keV = blobs.Eb1 * 1e3  # Eb1 from find_blob_energies is higher energy
            Eb2_keV = blobs.Eb2 * 1e3  # Eb2 from find_blob_energies is lower energy
            push!(blob_energies, BlobEnergies(Eb1_keV, Eb2_keV, r.event_id))
        catch e
            # Skip tracks that fail blob analysis
            continue
        end
    end

    return blob_energies
end

"""
    compute_efficiency(blob_energies, eb2_cut)

Compute efficiency for events passing Eb2 >= eb2_cut.
Eb2 is the lower energy blob.
"""
function compute_efficiency(blob_energies::Vector{BlobEnergies}, eb2_cut::Float64)
    n_total = length(blob_energies)
    if n_total == 0
        return 0.0
    end
    n_pass = count(b -> b.Eb2 >= eb2_cut, blob_energies)
    return n_pass / n_total
end

"""
    compute_fom(eff_signal, eff_background; min_signal_eff=0.01)

Compute figure of merit: ε_signal / √ε_background
Returns 0 if signal efficiency drops below min_signal_eff or background is zero.
"""
function compute_fom(eff_signal::Float64, eff_background::Float64; min_signal_eff::Float64=0.01)
    if eff_signal < min_signal_eff || eff_background <= 0
        return 0.0
    end
    return eff_signal / sqrt(eff_background)
end

"""
    create_scatter_plots(bb_blobs, ele_blobs, radius, output_dir)

Create Eb1 vs Eb2 scatter plots for signal and background.
"""
function create_scatter_plots(bb_blobs::Vector{BlobEnergies}, ele_blobs::Vector{BlobEnergies},
                               radius::Float64, output_dir::String)

    # Extract data
    bb_Eb1 = [b.Eb1 for b in bb_blobs]
    bb_Eb2 = [b.Eb2 for b in bb_blobs]

    ele_Eb1 = [b.Eb1 for b in ele_blobs]
    ele_Eb2 = [b.Eb2 for b in ele_blobs]

    # Determine energy range for consistent axes
    all_E = vcat(bb_Eb1, bb_Eb2, ele_Eb1, ele_Eb2)
    E_max = min(maximum(all_E) * 1.1, 1500)

    # Plot 1: Eb1 vs Eb2 for bb0nu (signal)
    p1 = scatter(bb_Eb2, bb_Eb1,
                 xlabel="Eb2 [keV]", ylabel="Eb1 [keV]",
                 title="bb0nu (signal) R=$(Int(radius))mm (n=$(length(bb_blobs)))",
                 label="", alpha=0.6, markersize=4, color=:blue,
                 xlims=(0, E_max), ylims=(0, E_max),
                 size=(600, 550),
                 guidefontsize=12, tickfontsize=10, titlefontsize=11,
                 left_margin=8Plots.mm, bottom_margin=5Plots.mm)
    plot!(p1, [0, E_max], [0, E_max], label="Eb1=Eb2", linestyle=:dash, color=:gray)

    bb_plot_file = joinpath(output_dir, "Eb1_vs_Eb2_bb0nu_R$(Int(radius))mm.png")
    savefig(p1, bb_plot_file)

    # Plot 2: Eb1 vs Eb2 for electrons (background)
    p2 = scatter(ele_Eb2, ele_Eb1,
                 xlabel="Eb2 [keV]", ylabel="Eb1 [keV]",
                 title="electrons (background) R=$(Int(radius))mm (n=$(length(ele_blobs)))",
                 label="", alpha=0.6, markersize=4, color=:red,
                 xlims=(0, E_max), ylims=(0, E_max),
                 size=(600, 550),
                 guidefontsize=12, tickfontsize=10, titlefontsize=11,
                 left_margin=8Plots.mm, bottom_margin=5Plots.mm)
    plot!(p2, [0, E_max], [0, E_max], label="Eb1=Eb2", linestyle=:dash, color=:gray)

    ele_plot_file = joinpath(output_dir, "Eb1_vs_Eb2_electrons_R$(Int(radius))mm.png")
    savefig(p2, ele_plot_file)

    # Combined plot
    p_combined = plot(p1, p2, layout=(1, 2), size=(1300, 550))
    combined_file = joinpath(output_dir, "Eb1_vs_Eb2_combined_R$(Int(radius))mm.png")
    savefig(p_combined, combined_file)

    println("Saved: $bb_plot_file")
    println("Saved: $ele_plot_file")
    println("Saved: $combined_file")

    return p_combined
end

"""
    run_fom_analysis(bb_dir, bb_pattern, ele_dir, ele_pattern; kwargs...)

Run the figure of merit analysis.
"""
function run_fom_analysis(bb_dir::String, bb_pattern::String,
                          ele_dir::String, ele_pattern::String;
                          radius::Float64=12.0,
                          eb2_min::Float64=50.0,
                          eb2_max::Float64=500.0,
                          eb2_step::Float64=25.0,
                          output_dir::String="fom_results")

    # Create output directory if it doesn't exist
    if !isdir(output_dir)
        mkpath(output_dir)
        println("Created output directory: $output_dir")
    end

    println("="^70)
    println("SIGNAL-TO-NOISE FIGURE OF MERIT ANALYSIS")
    println("="^70)
    println("Sphere radius: $radius mm")
    println("Eb2 scan: $eb2_min to $eb2_max keV (step $eb2_step keV)")
    println("Output directory: $output_dir")

    # Load bb0nu results (signal)
    println("\nLoading bb0nu (signal) data...")
    bb_results, bb_metadata = Petit.chain_reco_results(bb_dir, bb_pattern)
    println("  Loaded $(length(bb_results)) bb0nu tracks")

    # Load electron results (background)
    println("\nLoading electron (background) data...")
    ele_results, ele_metadata = Petit.chain_reco_results(ele_dir, ele_pattern)
    println("  Loaded $(length(ele_results)) electron tracks")

    if isempty(bb_results) || isempty(ele_results)
        error("No data loaded. Check file paths and patterns.")
    end

    # Compute blob energies
    println("\nComputing blob energies at R = $radius mm...")
    bb_blobs = compute_blob_energies(bb_results, radius)
    ele_blobs = compute_blob_energies(ele_results, radius)
    println("  bb0nu:     $(length(bb_blobs)) valid tracks")
    println("  electrons: $(length(ele_blobs)) valid tracks")

    # Create scatter plots
    println("\nCreating scatter plots...")
    create_scatter_plots(bb_blobs, ele_blobs, radius, output_dir)

    # Scan Eb2 cuts
    println("\n" * "-"^70)
    println("Eb2 CUT SCAN")
    println("-"^70)
    @printf("%-12s %12s %12s %12s\n", "Eb2 (keV)", "ε_bb0nu", "ε_electrons", "FOM")
    println("-"^70)

    eb2_cuts = Float64[]
    eff_bb = Float64[]
    eff_ele = Float64[]
    fom_values = Float64[]

    eb2 = eb2_min
    while eb2 <= eb2_max
        ε_bb = compute_efficiency(bb_blobs, eb2)
        ε_ele = compute_efficiency(ele_blobs, eb2)
        fom = compute_fom(ε_bb, ε_ele)

        @printf("%-12.1f %12.3f %12.3f %12.3f\n", eb2, ε_bb, ε_ele, fom)

        push!(eb2_cuts, eb2)
        push!(eff_bb, ε_bb)
        push!(eff_ele, ε_ele)
        push!(fom_values, fom)

        # Stop scan if signal efficiency drops below threshold
        if ε_bb < 0.01
            println("\nStopping scan: signal efficiency dropped below 1%")
            break
        end

        eb2 += eb2_step
    end

    # Find best FOM
    best_idx = argmax(fom_values)
    best_eb2 = eb2_cuts[best_idx]
    best_fom = fom_values[best_idx]
    best_eff_bb = eff_bb[best_idx]
    best_eff_ele = eff_ele[best_idx]

    println("-"^70)
    println("\nBest FOM:")
    println("  Eb2 cut:      $(best_eb2) keV")
    println("  ε_bb0nu:      $(round(best_eff_bb, digits=3))")
    println("  ε_electrons:  $(round(best_eff_ele, digits=3))")
    println("  FOM:          $(round(best_fom, digits=3))")

    # Create FOM plots
    println("\nGenerating FOM plots...")

    # Plot 1: Efficiency vs Eb2 cut
    p1 = plot(eb2_cuts, eff_bb,
              xlabel="Eb2 cut [keV]", ylabel="Efficiency",
              title="Efficiency vs Eb2 cut",
              label="bb0nu (signal)",
              marker=:circle, linewidth=2, color=:blue,
              grid=true, legend=:right,
              guidefontsize=11, tickfontsize=9, titlefontsize=11,
              left_margin=8Plots.mm, bottom_margin=5Plots.mm)
    plot!(p1, eb2_cuts, eff_ele,
          label="electrons (bkg)",
          marker=:square, linewidth=2, color=:red)
    vline!(p1, [best_eb2], linestyle=:dash, color=:black, label="Best cut")

    # Plot 2: FOM vs Eb2 cut
    p2 = plot(eb2_cuts, fom_values,
              xlabel="Eb2 cut [keV]", ylabel="FOM (ε_s/√ε_b)",
              title="Figure of Merit vs Eb2 cut",
              label="",
              marker=:circle, linewidth=2, color=:green,
              grid=true,
              guidefontsize=11, tickfontsize=9, titlefontsize=11,
              left_margin=8Plots.mm, bottom_margin=5Plots.mm)
    vline!(p2, [best_eb2], linestyle=:dash, color=:black, label="Best cut")
    scatter!(p2, [best_eb2], [best_fom], markersize=8, color=:red, label="Max FOM=$(round(best_fom, digits=2))")

    # Plot 3: Eb2 distributions
    bb_Eb2 = [b.Eb2 for b in bb_blobs]
    ele_Eb2 = [b.Eb2 for b in ele_blobs]

    p3 = histogram(bb_Eb2,
                   bins=range(0, maximum(vcat(bb_Eb2, ele_Eb2))*1.1, length=30),
                   alpha=0.6, label="bb0nu",
                   xlabel="Eb2 [keV]", ylabel="Probability",
                   title="Eb2 Distribution",
                   color=:blue, normalize=:probability,
                   guidefontsize=11, tickfontsize=9, titlefontsize=11,
                   left_margin=8Plots.mm, bottom_margin=5Plots.mm)
    histogram!(p3, ele_Eb2,
               bins=range(0, maximum(vcat(bb_Eb2, ele_Eb2))*1.1, length=30),
               alpha=0.6, label="electrons",
               color=:red, normalize=:probability)
    vline!(p3, [best_eb2], linestyle=:dash, color=:black, label="Best cut")

    # Combined plot
    p_combined = plot(p1, p2, p3, layout=(1, 3), size=(1600, 500))

    fom_plot_file = joinpath(output_dir, "fom_analysis_R$(Int(radius))mm.png")
    savefig(p_combined, fom_plot_file)
    println("Saved: $fom_plot_file")

    display(p_combined)

    # Save results table
    results_df = DataFrame(
        eb2_cut_keV = eb2_cuts,
        eff_bb0nu = eff_bb,
        eff_electrons = eff_ele,
        fom = fom_values
    )
    results_file = joinpath(output_dir, "fom_results_R$(Int(radius))mm.csv")
    Petit.CSV.write(results_file, results_df)
    println("Saved: $results_file")

    return (best_eb2=best_eb2, best_fom=best_fom,
            best_eff_bb=best_eff_bb, best_eff_ele=best_eff_ele,
            results_df=results_df)
end

"""
    main()

Main function for command-line usage.
"""
function main()
    # Default paths
    bb_dir = joinpath(ENV["DATA"], "HD5t/itaca/bb0nu")
    bb_pattern = "bb0nu_test_reco_th_*.h5"

    ele_dir = joinpath(ENV["DATA"], "HD5t/itaca/xe137")
    ele_pattern = "electrons_test_reco_th_*.h5"

    # Default parameters
    output_dir = "fom_results"
    radius = 12.0
    eb2_min = 50.0
    eb2_max = 500.0
    eb2_step = 25.0

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
                if key == "radius"
                    radius = parse(Float64, value)
                elseif key == "eb2min"
                    eb2_min = parse(Float64, value)
                elseif key == "eb2max"
                    eb2_max = parse(Float64, value)
                elseif key == "eb2step"
                    eb2_step = parse(Float64, value)
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
        println("Usage: julia ston_fom.jl <bb_dir> <bb_pattern> <ele_dir> <ele_pattern> [options]")
        println("\nOptional arguments:")
        println("  --radius=X     Sphere radius in mm (default: 12.0)")
        println("  --eb2min=X     Minimum Eb2 cut in keV (default: 50.0)")
        println("  --eb2max=X     Maximum Eb2 cut in keV (default: 500.0)")
        println("  --eb2step=X    Eb2 step size in keV (default: 25.0)")
        println("  --outdir=X     Output directory for plots (default: fom_results)")
        println("\nUsing default paths:")
        println("  bb0nu:     $bb_dir / $bb_pattern")
        println("  electrons: $ele_dir / $ele_pattern")
    end

    # Run analysis
    run_fom_analysis(bb_dir, bb_pattern, ele_dir, ele_pattern;
                     radius=radius,
                     eb2_min=eb2_min,
                     eb2_max=eb2_max,
                     eb2_step=eb2_step,
                     output_dir=output_dir)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
