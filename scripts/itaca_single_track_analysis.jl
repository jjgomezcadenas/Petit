#!/usr/bin/env julia

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)

using Printf
using DataFrames
using CSV
using Glob
using ArgParse
using Plots
using Dates

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
include(joinpath(pdir, "src", "itaca_aux.jl"))

# Import Petit module functions
import .Petit

const cmdir = joinpath(ENV["DATA"], "HD5t/itaca")

#=============================================================================
# Helper Functions
=============================================================================#

"""
    save_plot(plt, outdir, nevent, name)

Save plot to file without displaying.
"""
function save_plot(plt, outdir, nevent, name)
    mkpath(outdir)
    filename = joinpath(outdir, "event_$(nevent)_$(name).png")
    savefig(plt, filename)
    println("  Saved: $filename")
end

"""
    save_and_continue(plt, outdir, nevent, name)

Save plot to file, display, and wait for Enter.
"""
function save_and_continue(plt, outdir, nevent, name)
    mkpath(outdir)
    filename = joinpath(outdir, "event_$(nevent)_$(name).png")
    savefig(plt, filename)
    display(plt)
    println("\nSaved: $filename")
    println("Press Enter to continue...")
    readline()
end

"""
    get_sigma(particle_type, ldrft; dt, dl, tK, edrift, Pbar)

Compute transverse and longitudinal sigma based on particle type.
"""
function get_sigma(particle_type, ldrft;
				   dt = 3.5, dl = 0.9,
				   tK = 297.0, edrift = 500.0, Pbar=15.0)

	if particle_type == "ion"
		σt =  Petit.sigma_t_ion_mm(tK, ldrft, edrift)
		σl =  0.0
	else
		σt =  Petit.sigma_t_mm(ldrft, Pbar; dtmm=dt)
		σl =  Petit.sigma_l_mm(ldrft, Pbar; dlmm=dl)
	end
	σt, σl
end

"""
    get_voxel_size_and_distance(ldrft, σt)

Compute voxel size, MC voxel size, and max distance based on diffusion.
"""
function get_voxel_size_and_distance(ldrft, σt)
	if ldrft > 50.0
		voxel_scale = 2.0
		voxel_distance_scale = 1.5
	else
		voxel_scale = 3.0
		voxel_distance_scale = 2.0
	end

	voxel_size = σt * voxel_scale
	mcvox_size = voxel_size
	max_distance = voxel_size * voxel_distance_scale
	(voxel_size, mcvox_size, max_distance)
end

"""
    get_energy_threshold(particle_type; energy_threshold_ions, energy_threshold_keV)

Get energy threshold based on particle type.
"""
function get_energy_threshold(particle_type;
							  energy_threshold_ions =  10.0,
							  energy_threshold_keV =  10.0)
	f = 1e+5/2.5 # ions per MeV
	fkeV = f*1e-3 # ions per keV

	if particle_type == "ion"
		energy_threshold_keV = energy_threshold_ions/fkeV
	end
	energy_threshold_keV
end

"""
    print_event_results(nevent, track_length, blobs, extreme_dists, peaks, pk)

Pretty print the results for an event (CSV row contents).
Note: blobs.Eb1 and blobs.Eb2 are already in keV.
"""
function print_event_results(nevent, track_length, blobs, extreme_dists, peaks, pk)
    println("┌─────────────────────────────────────────────────────────────────┐")
    println("│  Event $nevent Results                                          ")
    println("├─────────────────────────────────────────────────────────────────┤")
    println("│  Track length:    $(lpad(@sprintf("%.2f", track_length), 8)) mm                          │")
    println("├─────────────────────────────────────────────────────────────────┤")
    println("│  Blob Analysis:                                                 │")
    println("│    Eb1 (high):    $(lpad(@sprintf("%.1f", blobs.Eb1), 8)) keV                        │")
    println("│    Eb2 (low):     $(lpad(@sprintf("%.1f", blobs.Eb2), 8)) keV                        │")
    println("│    Asymmetry:     $(lpad(@sprintf("%.4f", blobs.asymmetry), 8))                            │")
    println("├─────────────────────────────────────────────────────────────────┤")
    println("│  Extreme Distances (RECO vs MC):                                │")
    println("│    d1:            $(lpad(@sprintf("%.2f", extreme_dists.d1), 8)) mm                          │")
    println("│    d2:            $(lpad(@sprintf("%.2f", extreme_dists.d2), 8)) mm                          │")
    println("├─────────────────────────────────────────────────────────────────┤")
    println("│  KDE Peak Analysis:                                             │")
    println("│    N peaks:       $(lpad(length(peaks.indices), 8))                            │")
    println("│    Peak1: L=$(lpad(@sprintf("%.1f", pk.peak1_left), 5)) R=$(lpad(@sprintf("%.1f", pk.peak1_right), 5)) prom=$(lpad(@sprintf("%.3f", pk.peak1_prom), 5))     │")
    println("│    Peak2: L=$(lpad(@sprintf("%.1f", pk.peak2_left), 5)) R=$(lpad(@sprintf("%.1f", pk.peak2_right), 5)) prom=$(lpad(@sprintf("%.3f", pk.peak2_prom), 5))     │")
    println("└─────────────────────────────────────────────────────────────────┘")
end

"""
    print_blobs(blobs)

Print blob analysis results to console.
Note: blobs.Eb1 and blobs.Eb2 are already in keV.
"""
function print_blobs(blobs)
    println("═══════════════════════════════════════")
    println("          Blob Analysis")
    println("═══════════════════════════════════════")
    println("Blob1 (high energy):")
    println("  position = ($(round(blobs.blob1.x, digits=2)), $(round(blobs.blob1.y, digits=2)), $(round(blobs.blob1.z, digits=2))) mm")
    println("  energy   = $(round(blobs.Eb1, digits=1)) keV")
    println("───────────────────────────────────────")
    println("Blob2 (low energy):")
    println("  position = ($(round(blobs.blob2.x, digits=2)), $(round(blobs.blob2.y, digits=2)), $(round(blobs.blob2.z, digits=2))) mm")
    println("  energy   = $(round(blobs.Eb2, digits=1)) keV")
    println("───────────────────────────────────────")
    println("Asymmetry  = $(round(blobs.asymmetry, digits=3))")
    println("═══════════════════════════════════════")
end

"""
    print_peaks(peaks)

Print peak analysis results to console.
"""
function print_peaks(peaks)
    println("═══════════════════════════════════════")
    println("          Peak Analysis")
    println("═══════════════════════════════════════")
    println("Number of peaks: $(length(peaks.indices))")
    for i in eachindex(peaks.indices)
        println("───────────────────────────────────────")
        println("Peak $i:")
        println("  position    = $(round(peaks.positions[i], digits=2)) mm")
        println("  prominence  = $(round(peaks.proms[i], digits=4))")
        println("  width       = $(round(peaks.widths[i], digits=2))")
        println("  left edge   = $(round(peaks.lefts[i], digits=2))")
        println("  right edge  = $(round(peaks.rights[i], digits=2))")
    end
    println("═══════════════════════════════════════")
end

"""
    plot_kde_comparison(reco_kde, mc_kde, nevent, kde_bandwidth, n_kde_eval)

Plot KDE and histogram comparison between RECO and MC.
"""
function plot_kde_comparison(reco_kde, mc_kde, nevent, kde_bandwidth, n_kde_eval)
    # KDE comparison
    p1 = plot(reco_kde.kde_s, reco_kde.kde_f,
              xlabel="Arc length s (mm)",
              ylabel="KDE Energy density f(s)",
              title="Longitudinal Energy Density (Event $nevent)",
              label="RECO KDE (h=$(round(kde_bandwidth, digits=1))mm)",
              linewidth=2,
              color=:blue)
    plot!(p1, mc_kde.kde_s, mc_kde.kde_f,
          label="MC KDE",
          linewidth=2,
          color=:red)

    # Histogram comparison
    p2 = histogram(reco_kde.s, weights=reco_kde.E .* 1e3,
                   bins=range(0, reco_kde.track_length, length=n_kde_eval+1),
                   xlabel="Arc length s (mm)",
                   ylabel="Energy (keV)",
                   title="Energy vs S (Event $nevent)",
                   label="RECO",
                   fillalpha=0.5,
                   color=:blue)
    histogram!(p2, mc_kde.s, weights=mc_kde.E .* 1e3,
               bins=range(0, mc_kde.track_length, length=n_kde_eval+1),
               label="MC",
               fillalpha=0.5,
               color=:red)

    return plot(p1, p2, layout=(1, 2), size=(1200, 500))
end

"""
    save_metadata(outdir, params, n_processed, n_single_track)

Save metadata CSV file with run parameters.
"""
function save_metadata(outdir::String, params::Dict, n_processed::Int, n_single_track::Int)
    mkpath(outdir)

    metadata_df = DataFrame(
        parameter = String[],
        value = String[]
    )

    # Add all parameters
    for (k, v) in params
        push!(metadata_df, (parameter=k, value=string(v)))
    end

    # Add processing stats
    push!(metadata_df, (parameter="n_events_processed", value=string(n_processed)))
    push!(metadata_df, (parameter="n_events_single_track", value=string(n_single_track)))
    push!(metadata_df, (parameter="timestamp", value=string(now())))

    CSV.write(joinpath(outdir, "metadata.csv"), metadata_df)
    println("\nMetadata saved to: $(joinpath(outdir, "metadata.csv"))")
end

#=============================================================================
# Main Script
=============================================================================#

function main(; ievent::Int=1, levent::Int=1,
               input_file::String="bb0nu/bb0nu_15bar_p1.h5",
               particle_type::String="ion",
               ldrft::Float64=100.0,
               dt::Float64=3.5, dl::Float64=0.9,
               tK::Float64=297.0, edrift::Float64=500.0, Pbar::Float64=15.0,
               energy_threshold_ions::Float64=10.0,
               energy_threshold_keV::Float64=10.0,
               nbins::Int=100, nsigma::Float64=3.0,
               n_kde_eval::Int=200, Rb::Float64=10.0,
               outdir::String="results",
               interactive::Bool=false,
               nplot::Int=0,
               print_level::String="verbose")

    println("╔═══════════════════════════════════════════════════════════════╗")
    println("║         ITACA Single Track Analysis                          ║")
    println("╚═══════════════════════════════════════════════════════════════╝")

    # Compute diffusion parameters
    σt, σl = get_sigma(particle_type, ldrft;
                       dt=dt, dl=dl, tK=tK, edrift=edrift, Pbar=Pbar)

    voxel_size, mcvox_size, max_distance = get_voxel_size_and_distance(ldrft, σt)

    eth = get_energy_threshold(particle_type;
                               energy_threshold_ions=energy_threshold_ions,
                               energy_threshold_keV=energy_threshold_keV)

    kde_bandwidth = 2.0 * voxel_size
    mc_kde_bandwidth = 2.0 * mcvox_size

    dfpars = Petit.DiffusionParams(ldrft, σt, σl, voxel_size, max_distance,
                                   eth, nbins, nsigma)

    # Print level flags: "quiet", "verbose", or "very_verbose"
    is_quiet = (print_level == "quiet")
    is_verbose = (print_level == "verbose")
    is_very_verbose = (print_level == "very_verbose")

    if is_verbose || is_very_verbose
        println("\n### Diffusion Parameters")
        diffusion_params_print(dfpars)
        println("\nInteractive mode: $interactive")
        println("Plot every N events: $(nplot == 0 ? "disabled" : nplot)")
        println("Print level: $print_level")
    end

    # Store parameters for metadata
    params = Dict(
        "input_file" => input_file,
        "particle_type" => particle_type,
        "ldrft_cm" => ldrft,
        "sigma_t_mm" => σt,
        "sigma_l_mm" => σl,
        "voxel_size_mm" => voxel_size,
        "max_distance_mm" => max_distance,
        "energy_threshold_keV" => eth,
        "kde_bandwidth_mm" => kde_bandwidth,
        "blob_radius_mm" => Rb,
        "ievent" => ievent,
        "levent" => levent
    )

    # Create output directories
    mkpath(outdir)
    plots_dir = joinpath(outdir, "plots")
    mkpath(plots_dir)

    # Load data
    if is_very_verbose
        println("\n### Loading data from: $input_file")
    end
    hitsdf = Petit.load_data(input_file, cmdir)

    # Initialize results DataFrame
    results_df = DataFrame(
        event = Int[],
        track_length_mm = Float64[],
        Eb1_keV = Float64[],
        Eb2_keV = Float64[],
        asymmetry = Float64[],
        d1_mm = Float64[],
        d2_mm = Float64[],
        n_peaks = Int[],
        peak1_left = Float64[],
        peak1_right = Float64[],
        peak1_prom = Float64[],
        peak2_left = Float64[],
        peak2_right = Float64[],
        peak2_prom = Float64[]
    )

    n_processed = 0
    n_single_track = 0

    # Process events
    for nevent in ievent:levent
        n_processed += 1

        # Determine if we should save plots for this event
        save_plots_this_event = (nplot > 0) && (mod(nevent, nplot) == 0)

        # Print event header based on print_level
        if is_very_verbose
            println("\n")
            println("╔═══════════════════════════════════════════════════════════════╗")
            println("║  Processing Event $nevent                                      ")
            println("╚═══════════════════════════════════════════════════════════════╝")
        elseif is_verbose
            print("Event $nevent: ")
        elseif save_plots_this_event
            # quiet mode: only print on plot events
            print("Event $nevent... ")
        end

        # Get event
        event_df = Petit.get_event(hitsdf, nevent)

        if interactive || save_plots_this_event
            plt = Petit.plot_event(event_df)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "01_mc_hits")
            else
                save_plot(plt, plots_dir, nevent, "01_mc_hits")
            end
        end

        # Compute MC path
        if is_very_verbose
            println("\n### Computing MC path")
        end
        mc_path = Petit.compute_mc_path(event_df, mcvox_size)

        if interactive || save_plots_this_event
            plt = Petit.plot_event(event_df; mc_path=mc_path)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "02_mc_path")
            else
                save_plot(plt, plots_dir, nevent, "02_mc_path")
            end
        end

        # Transform hits
        if is_very_verbose
            println("### Transforming hits")
        end
        event_mc = Petit.transform_hits_df(event_df)

        # Diffuse event
        if is_very_verbose
            println("### Diffusing event")
        end
        diffused_df = Petit.diffuse_xyz_image_mc(event_mc;
                                                 sigma_t_mm=σt,
                                                 sigma_l_mm=σl,
                                                 nbins=nbins,
                                                 nsigma=nsigma)

        if interactive || save_plots_this_event
            plt = Petit.plot_hits(diffused_df, energy_column=:electrons)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "03_diffused")
            else
                save_plot(plt, plots_dir, nevent, "03_diffused")
            end
        end

        # Voxelize event
        if is_very_verbose
            println("### Voxelizing event")
        end
        voxels = Petit.voxelize_event(diffused_df, voxel_size)

        if interactive || save_plots_this_event
            plt = Petit.plot_event(voxels)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "04_voxels")
            else
                save_plot(plt, plots_dir, nevent, "04_voxels")
            end
        end

        # Make tracks
        if is_very_verbose
            println("### Making tracks")
        end
        tracks = Petit.make_tracks(voxels;
                                   max_distance_mm=max_distance,
                                   energy_threshold_kev=eth,
                                   diffusion=dfpars)

        n_tracks = length(tracks)
        if is_very_verbose
            println("Number of tracks found: $n_tracks")
        end

        if n_tracks != 1
            if is_very_verbose
                println("⚠ Number of tracks ≠ 1, skipping event...")
            elseif is_verbose
                println("$n_tracks tracks - rejected")
            elseif save_plots_this_event
                println("rejected ($n_tracks tracks)")
            end
            if save_plots_this_event
                savefig(Petit.plot_event(voxels),
                        joinpath(plots_dir, "event_$(nevent)_multitrack.png"))
            end
            continue
        end

        n_single_track += 1
        track = tracks[1]

        # Walk track
        if is_very_verbose
            println("### Walking track from extremes")
        end
        walk_result = Petit.walk_track_from_extremes(track)
        if interactive
            walk_result_print(walk_result)
        end

        # Get path
        path = Petit.get_raw_path(track, walk_result.path_indices)
        track_length = path.s[end]
        if interactive
            path_print(path)
        end

        # Compute extreme distances
        if is_very_verbose
            println("### Computing extreme distances (RECO vs MC)")
        end
        extreme_dists = Petit.compute_extreme_distances(path, mc_path)
        if interactive
            extreme_distances_print(extreme_dists)
        end

        if interactive || save_plots_this_event
            plt = Petit.plot_track_with_paths(track, path, mc_path; show_distances=true)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "05_track_paths")
            else
                save_plot(plt, plots_dir, nevent, "05_track_paths")
            end
        end

        # KDE analysis
        if is_very_verbose
            println("### KDE Analysis")
        end
        reco_kde = Petit.get_reco_kde(track, path; bandwidth=kde_bandwidth,
                                      n_eval=n_kde_eval)
        mc_kde = Petit.get_mc_kde(mc_path; bandwidth=mc_kde_bandwidth,
                                  n_eval=n_kde_eval)

        if interactive || save_plots_this_event
            plt = plot_kde_comparison(reco_kde, mc_kde, nevent, kde_bandwidth, n_kde_eval)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "06_kde_comparison")
            else
                save_plot(plt, plots_dir, nevent, "06_kde_comparison")
            end
        end

        # Find peaks
        if is_very_verbose
            println("### Finding peaks")
        end
        peaks = Petit.find_peaks(reco_kde.kde_f, reco_kde.kde_s; prom_scale=0.2)
        if interactive
            print_peaks(peaks)
        end

        # Extract peak1 (leftmost) and peak2 (rightmost) info
        pk = kde_peaks(peaks, reco_kde.kde_f)
        if is_very_verbose
            println("    Peak1: left=$(round(pk.peak1_left, digits=1)), right=$(round(pk.peak1_right, digits=1)), prom=$(round(pk.peak1_prom, digits=3))")
            println("    Peak2: left=$(round(pk.peak2_left, digits=1)), right=$(round(pk.peak2_right, digits=1)), prom=$(round(pk.peak2_prom, digits=3))")
        end

        if interactive || save_plots_this_event
            plt = Petit.plot_kde_peaks(reco_kde.kde_s, reco_kde.kde_f, peaks;
                                       title="RECO KDE (Event $nevent)")
            if interactive
                save_and_continue(plt, plots_dir, nevent, "07_kde_peaks")
            else
                save_plot(plt, plots_dir, nevent, "07_kde_peaks")
            end
        end

        # Blob analysis
        if is_very_verbose
            println("### Blob Analysis")
        end
        blobs = Petit.find_blob_energies(track, path; radius=Rb)
        if interactive
            print_blobs(blobs)
        end

        if interactive || save_plots_this_event
            plt = Petit.plot_reco_track_with_voxels_and_spheres(track, path, blobs, Rb)
            if interactive
                save_and_continue(plt, plots_dir, nevent, "08_blobs")
            else
                save_plot(plt, plots_dir, nevent, "08_blobs")
            end
        end

        # Add row to results DataFrame
        push!(results_df, (
            event = nevent,
            track_length_mm = track_length,
            Eb1_keV = blobs.Eb1,  # Already in keV
            Eb2_keV = blobs.Eb2,  # Already in keV
            asymmetry = blobs.asymmetry,
            d1_mm = extreme_dists.d1,
            d2_mm = extreme_dists.d2,
            n_peaks = length(peaks.indices),
            peak1_left = pk.peak1_left,
            peak1_right = pk.peak1_right,
            peak1_prom = pk.peak1_prom,
            peak2_left = pk.peak2_left,
            peak2_right = pk.peak2_right,
            peak2_prom = pk.peak2_prom
        ))

        # Print completion message based on print_level
        if is_very_verbose
            print_event_results(nevent, track_length, blobs, extreme_dists, peaks, pk)
            println("✓ Event $nevent completed")
        elseif is_verbose
            println("1 track")
            print_event_results(nevent, track_length, blobs, extreme_dists, peaks, pk)
        elseif save_plots_this_event
            println("done")
        end
    end

    # Save results CSV
    results_file = joinpath(outdir, "analysis_results.csv")
    CSV.write(results_file, results_df)
    println("\n✓ Results saved to: $results_file")

    # Save metadata
    save_metadata(outdir, params, n_processed, n_single_track)

    println("\n╔═══════════════════════════════════════════════════════════════╗")
    println("║  Analysis Complete                                            ║")
    println("╠═══════════════════════════════════════════════════════════════╣")
    println("║  Events processed:    $(lpad(n_processed, 5))                              ║")
    println("║  Single-track events: $(lpad(n_single_track, 5))                              ║")
    println("║  Output directory:    $(outdir)")
    println("╚═══════════════════════════════════════════════════════════════╝")
end


#=============================================================================
# Command Line Interface
=============================================================================#

if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParseSettings(description="ITACA Single Track Analysis")

    @add_arg_table! s begin
        "--ievent", "-i"
            help = "Initial event number"
            arg_type = Int
            default = 1
        "--levent", "-e"
            help = "Last event number"
            arg_type = Int
            default = 1
        "--input", "-f"
            help = "Input file (relative to DATA/HD5t/itaca)"
            arg_type = String
            default = "bb0nu/bb0nu_15bar_p1.h5"
        "--particle", "-p"
            help = "Particle type: 'ion' or 'electron'"
            arg_type = String
            default = "ion"
        "--ldrft", "-l"
            help = "Drift length in cm"
            arg_type = Float64
            default = 100.0
        "--dt"
            help = "Transverse diffusion coefficient (mm/√cm)"
            arg_type = Float64
            default = 3.5
        "--dl"
            help = "Longitudinal diffusion coefficient (mm/√cm)"
            arg_type = Float64
            default = 0.9
        "--tK"
            help = "Temperature in Kelvin"
            arg_type = Float64
            default = 297.0
        "--edrift"
            help = "Drift field in V/cm"
            arg_type = Float64
            default = 500.0
        "--Pbar"
            help = "Pressure in bar"
            arg_type = Float64
            default = 15.0
        "--eth-ion"
            help = "Energy threshold for ions (keV)"
            arg_type = Float64
            default = 10.0
        "--eth"
            help = "Energy threshold (keV)"
            arg_type = Float64
            default = 10.0
        "--nbins"
            help = "Number of bins for diffusion"
            arg_type = Int
            default = 100
        "--nsigma"
            help = "Number of sigma for diffusion"
            arg_type = Float64
            default = 3.0
        "--nkde"
            help = "Number of KDE evaluation points"
            arg_type = Int
            default = 200
        "--Rb"
            help = "Blob radius in mm"
            arg_type = Float64
            default = 10.0
        "--outdir", "-o"
            help = "Output directory for results"
            arg_type = String
            default = "results"
        "--interactive"
            help = "Enable interactive mode (display plots and wait for Enter)"
            action = :store_true
        "--nplot"
            help = "Save plots every N events (0 = disabled)"
            arg_type = Int
            default = 0
        "--print_level"
            help = "Print level: 'quiet' (minimal), 'verbose' (track count + results), 'very_verbose' (all steps)"
            arg_type = String
            default = "verbose"
    end

    args = parse_args(s)

    main(; ievent=args["ievent"],
           levent=args["levent"],
           input_file=args["input"],
           particle_type=args["particle"],
           ldrft=args["ldrft"],
           dt=args["dt"],
           dl=args["dl"],
           tK=args["tK"],
           edrift=args["edrift"],
           Pbar=args["Pbar"],
           energy_threshold_ions=args["eth-ion"],
           energy_threshold_keV=args["eth"],
           nbins=args["nbins"],
           nsigma=args["nsigma"],
           n_kde_eval=args["nkde"],
           Rb=args["Rb"],
           outdir=args["outdir"],
           interactive=args["interactive"],
           nplot=args["nplot"],
           print_level=args["print_level"])
end
