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
import .Petit
# Helper functions (get_sigma, get_voxel_size_and_distance, get_energy_threshold)
# are now part of Petit module via itaca_functions.jl

const cmdir = joinpath(ENV["DATA"], "HD5t/itaca")

#=============================================================================
# Helper Functions
=============================================================================#

"""
    save_plot(plt, plotdir, nevent, ldrft, name)

Save plot to plotdir with event and drift length in filename.
"""
function save_plot(plt, plotdir, nevent, ldrft, name)
    filename = joinpath(plotdir, "evt$(nevent)_L$(Int(ldrft))cm_$(name).png")
    savefig(plt, filename)
    println("  Saved: $filename")
end

"""
    plot_kde_comparison(reco_kde, mc_kde, nevent, ldrft, kde_bandwidth, n_kde_eval)

Plot KDE and histogram comparison between RECO and MC.
"""
function plot_kde_comparison(reco_kde, mc_kde, nevent, ldrft, kde_bandwidth, n_kde_eval)
    # KDE comparison
    p1 = plot(reco_kde.kde_s, reco_kde.kde_f,
              xlabel="Arc length s (mm)",
              ylabel="KDE Energy density f(s)",
              title="Event $nevent, L=$ldrft cm",
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
                   title="Energy vs S",
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
    save_results_csv(csvdir, nevent, ldrft, params, extreme_dists, peaks, blobs)

Save analysis results to CSV files with drift length in filename.
"""
function save_results_csv(csvdir, nevent, ldrft, params, extreme_dists, peaks, blobs)
    prefix = "evt$(nevent)_L$(Int(ldrft))cm"

    # Save parameters
    params_df = DataFrame(
        parameter = collect(keys(params)),
        value = collect(values(params))
    )
    CSV.write(joinpath(csvdir, "$(prefix)_params.csv"), params_df)

    # Save distances
    dists_df = DataFrame(
        event = nevent,
        ldrft_cm = ldrft,
        d1 = extreme_dists.d1,
        d2 = extreme_dists.d2,
        total = extreme_dists.total,
        pairing = String(extreme_dists.pairing)
    )
    CSV.write(joinpath(csvdir, "$(prefix)_distances.csv"), dists_df)

    # Save peaks (leftmost and rightmost)
    if length(peaks.indices) >= 2
        left_idx = argmin(peaks.lefts)
        right_idx = argmax(peaks.rights)

        peaks_df = DataFrame(
            event = [nevent, nevent],
            ldrft_cm = [ldrft, ldrft],
            peak_type = ["leftmost", "rightmost"],
            position = [peaks.positions[left_idx], peaks.positions[right_idx]],
            prominence = [peaks.proms[left_idx], peaks.proms[right_idx]],
            width = [peaks.widths[left_idx], peaks.widths[right_idx]],
            left_edge = [peaks.lefts[left_idx], peaks.lefts[right_idx]],
            right_edge = [peaks.rights[left_idx], peaks.rights[right_idx]]
        )
    elseif length(peaks.indices) == 1
        peaks_df = DataFrame(
            event = nevent,
            ldrft_cm = ldrft,
            peak_type = "single",
            position = peaks.positions[1],
            prominence = peaks.proms[1],
            width = peaks.widths[1],
            left_edge = peaks.lefts[1],
            right_edge = peaks.rights[1]
        )
    else
        peaks_df = DataFrame()
    end
    CSV.write(joinpath(csvdir, "$(prefix)_peaks.csv"), peaks_df)

    # Save blobs
    blobs_df = DataFrame(
        event = nevent,
        ldrft_cm = ldrft,
        Eb1_keV = blobs.Eb1,
        Eb2_keV = blobs.Eb2,
        asymmetry = blobs.asymmetry,
        blob1_x = blobs.blob1.x,
        blob1_y = blobs.blob1.y,
        blob1_z = blobs.blob1.z,
        blob2_x = blobs.blob2.x,
        blob2_y = blobs.blob2.y,
        blob2_z = blobs.blob2.z
    )
    CSV.write(joinpath(csvdir, "$(prefix)_blobs.csv"), blobs_df)
end

#=============================================================================
# Process single event at given drift length
=============================================================================#

function process_event_at_drift(event_df, nevent, ldrft, particle_type,
                                 plotdir, csvdir;
                                 dt=3.5, dl=0.9, tK=297.0, edrift=500.0, Pbar=15.0,
                                 energy_threshold_ions=10.0, energy_threshold_keV=10.0,
                                 nbins=100, nsigma=3.0, n_kde_eval=200, Rb=10.0)

    println("\n  ─── L = $ldrft cm ───")

    # Compute diffusion parameters
    σt, σl = Petit.get_sigma(particle_type, ldrft;
                       dt=dt, dl=dl, tK=tK, edrift=edrift, Pbar=Pbar)

    voxel_size, mcvox_size, max_distance = Petit.get_voxel_size_and_distance(ldrft, σt)

    eth = Petit.get_energy_threshold(particle_type;
                               energy_threshold_ions=energy_threshold_ions,
                               energy_threshold_keV=energy_threshold_keV)

    kde_bandwidth = 2.0 * voxel_size
    mc_kde_bandwidth = 2.0 * mcvox_size

    dfpars = Petit.DiffusionParams(ldrft, σt, σl, voxel_size, max_distance,
                                   eth, nbins, nsigma)

    println("    σt = $(round(σt, digits=2)) mm, voxel = $(round(voxel_size, digits=2)) mm")

    # Store parameters
    params = Dict(
        "ldrft_cm" => ldrft,
        "particle_type" => particle_type,
        "sigma_t_mm" => σt,
        "sigma_l_mm" => σl,
        "voxel_size_mm" => voxel_size,
        "max_distance_mm" => max_distance,
        "energy_threshold_keV" => eth,
        "kde_bandwidth_mm" => kde_bandwidth,
        "blob_radius_mm" => Rb
    )

    # Compute MC path
    mc_path = Petit.compute_mc_path(event_df, mcvox_size)

    # Transform and diffuse
    event_mc = Petit.transform_hits_df(event_df)
    diffused_df = Petit.diffuse_xyz_image_mc(event_mc;
                                             sigma_t_mm=σt,
                                             sigma_l_mm=σl,
                                             nbins=nbins,
                                             nsigma=nsigma)

    # Voxelize
    voxels = Petit.voxelize_event(diffused_df, voxel_size)

    # Make tracks
    tracks = Petit.make_tracks(voxels;
                               max_distance_mm=max_distance,
                               energy_threshold_kev=eth,
                               diffusion=dfpars)

    println("    Tracks found: $(length(tracks))")

    if length(tracks) != 1
        println("    ⚠ Skipping (tracks ≠ 1)")
        plt = Petit.plot_event(voxels)
        save_plot(plt, plotdir, nevent, ldrft, "multitrack")
        return nothing
    end

    track = tracks[1]

    # Walk track and get path
    walk_result = Petit.walk_track_from_extremes(track)
    path = Petit.get_raw_path(track, walk_result.path_indices)

    # Extreme distances
    extreme_dists = Petit.compute_extreme_distances(path, mc_path)
    println("    d1=$(round(extreme_dists.d1, digits=2)) mm, d2=$(round(extreme_dists.d2, digits=2)) mm")

    # Save track paths plot
    plt = Petit.plot_track_with_paths(track, path, mc_path; show_distances=true)
    save_plot(plt, plotdir, nevent, ldrft, "track_paths")

    # KDE analysis
    reco_kde = Petit.get_reco_kde(track, path; bandwidth=kde_bandwidth, n_eval=n_kde_eval)
    mc_kde = Petit.get_mc_kde(mc_path; bandwidth=mc_kde_bandwidth, n_eval=n_kde_eval)

    plt = plot_kde_comparison(reco_kde, mc_kde, nevent, ldrft, kde_bandwidth, n_kde_eval)
    save_plot(plt, plotdir, nevent, ldrft, "kde")

    # Find peaks
    peaks = Petit.find_peaks(reco_kde.kde_f, reco_kde.kde_s; prom_scale=0.2)

    plt = Petit.plot_kde_peaks(reco_kde.kde_s, reco_kde.kde_f, peaks;
                               title="RECO KDE (Event $nevent, L=$ldrft cm)")
    save_plot(plt, plotdir, nevent, ldrft, "peaks")

    # Extract peak1 (leftmost) and peak2 (rightmost) info
    pk = Petit.kde_peaks(peaks, reco_kde.kde_f)

    println("    Peak1: left=$(round(pk.peak1_left, digits=1)), right=$(round(pk.peak1_right, digits=1)), prom=$(round(pk.peak1_prom, digits=3))")
    println("    Peak2: left=$(round(pk.peak2_left, digits=1)), right=$(round(pk.peak2_right, digits=1)), prom=$(round(pk.peak2_prom, digits=3))")

    # Blob analysis
    blobs = Petit.find_blob_energies(track, path; radius=Rb)
    println("    Eb1=$(round(blobs.Eb1, digits=1)) keV, Eb2=$(round(blobs.Eb2, digits=1)) keV, asym=$(round(blobs.asymmetry, digits=3))")

    plt = Petit.plot_reco_track_with_voxels_and_spheres(track, path, blobs, Rb)
    save_plot(plt, plotdir, nevent, ldrft, "blobs")

    # Save CSV results
    save_results_csv(csvdir, nevent, ldrft, params, extreme_dists, peaks, blobs)

    return (d1=extreme_dists.d1, d2=extreme_dists.d2,
            Eb1=blobs.Eb1, Eb2=blobs.Eb2, asymmetry=blobs.asymmetry,
            peak1_left=pk.peak1_left, peak1_right=pk.peak1_right, peak1_prom=pk.peak1_prom,
            peak2_left=pk.peak2_left, peak2_right=pk.peak2_right, peak2_prom=pk.peak2_prom)
end

#=============================================================================
# Main Script
=============================================================================#

function main(; ievent::Int=1, levent::Int=1,
               input_file::String="bb0nu/bb0nu_15bar_p1.h5",
               particle_type::String="ion",
               ldrft_values::Vector{Float64}=collect(20.0:20.0:200.0),
               dt::Float64=3.5, dl::Float64=0.9,
               tK::Float64=297.0, edrift::Float64=500.0, Pbar::Float64=15.0,
               energy_threshold_ions::Float64=10.0,
               energy_threshold_keV::Float64=10.0,
               nbins::Int=100, nsigma::Float64=3.0,
               n_kde_eval::Int=200, Rb::Float64=10.0,
               outdir::String="resolution_scan")

    println("╔═══════════════════════════════════════════════════════════════╗")
    println("║      ITACA Resolution Scan - Single Track Analysis           ║")
    println("╚═══════════════════════════════════════════════════════════════╝")

    # Create output directories
    plotdir = joinpath(outdir, "plots")
    csvdir = joinpath(outdir, "csv")
    mkpath(plotdir)
    mkpath(csvdir)

    println("\nOutput directories:")
    println("  Plots: $plotdir")
    println("  CSV:   $csvdir")
    println("\nDrift lengths to scan: $ldrft_values cm")
    println("Events: $ievent to $levent")

    # Load data
    println("\n### Loading data from: $input_file")
    hitsdf = Petit.load_data(input_file, cmdir)

    # Collect summary results
    summary_results = DataFrame(
        event = Int[],
        ldrft_cm = Float64[],
        d1_mm = Float64[],
        d2_mm = Float64[],
        Eb1_keV = Float64[],
        Eb2_keV = Float64[],
        asymmetry = Float64[],
        peak1_left = Float64[],
        peak1_right = Float64[],
        peak1_prom = Float64[],
        peak2_left = Float64[],
        peak2_right = Float64[],
        peak2_prom = Float64[]
    )

    # Process events
    for nevent in ievent:levent
        println("\n╔═══════════════════════════════════════════════════════════════╗")
        println("║  Event $nevent                                                  ")
        println("╚═══════════════════════════════════════════════════════════════╝")

        event_df = Petit.get_event(hitsdf, nevent)

        # Save MC hits plot (once per event)
        plt = Petit.plot_event(event_df)
        savefig(plt, joinpath(plotdir, "evt$(nevent)_mc_hits.png"))

        # Scan over drift lengths
        for ldrft in ldrft_values
            result = process_event_at_drift(event_df, nevent, ldrft, particle_type,
                                            plotdir, csvdir;
                                            dt=dt, dl=dl, tK=tK, edrift=edrift, Pbar=Pbar,
                                            energy_threshold_ions=energy_threshold_ions,
                                            energy_threshold_keV=energy_threshold_keV,
                                            nbins=nbins, nsigma=nsigma,
                                            n_kde_eval=n_kde_eval, Rb=Rb)

            if result !== nothing
                push!(summary_results, (nevent, ldrft, result.d1, result.d2,
                                        result.Eb1, result.Eb2, result.asymmetry,
                                        result.peak1_left, result.peak1_right, result.peak1_prom,
                                        result.peak2_left, result.peak2_right, result.peak2_prom))
            end
        end
    end

    # Save summary CSV
    summary_file = joinpath(csvdir, "summary_all_events.csv")
    CSV.write(summary_file, summary_results)
    println("\n✓ Summary saved to: $summary_file")

    # Plot summary: d1+d2 vs drift length
    if nrow(summary_results) > 0
        summary_results.d_total = summary_results.d1_mm .+ summary_results.d2_mm

        plt = plot(xlabel="Drift length (cm)", ylabel="Distance (mm)",
                   title="Extreme Distance vs Drift Length", legend=:topleft)

        for evt in unique(summary_results.event)
            evt_data = filter(row -> row.event == evt, summary_results)
            plot!(plt, evt_data.ldrft_cm, evt_data.d_total,
                  marker=:circle, label="Event $evt (d1+d2)")
        end

        savefig(plt, joinpath(plotdir, "summary_distance_vs_drift.png"))
        println("✓ Summary plot saved")
    end

    println("\n╔═══════════════════════════════════════════════════════════════╗")
    println("║  Resolution Scan Complete                                     ║")
    println("╚═══════════════════════════════════════════════════════════════╝")
end

#=============================================================================
# Command Line Interface
=============================================================================#

if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParseSettings(description="ITACA Resolution Scan - Single Track Analysis")

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
            help = "Energy threshold for ions"
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
            default = "resolution_scan"
    end

    args = parse_args(s)

    main(; ievent=args["ievent"],
           levent=args["levent"],
           input_file=args["input"],
           particle_type=args["particle"],
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
           outdir=args["outdir"])
end
