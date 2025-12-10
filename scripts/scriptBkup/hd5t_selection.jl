#!/usr/bin/env julia

"""
HD5t Selection Script

This script performs event selection for signal and background events
and writes the results to JSON summary files.
"""

# Activate the project environment
using Pkg
Pkg.activate(dirname(@__DIR__))

using CSV
using DataFrames
using Plots
using Printf
using HDF5
using Statistics
using StatsBase
using Distributions
using StatsPlots

using Petit

# Configuration
cmdir = joinpath(ENV["DATA"], "HD5t")
pdir = joinpath(ENV["PROJECTS"], "Petit")
results_base_dir = joinpath(pdir, "AnalysisResults")
results_znubb = joinpath(results_base_dir, "0nubb")

# Analysis Parameters
@printf("=== HD5t Event Selection Analysis ===\n")
@printf("Analysis Parameters:\n")

# Radial cut (mm): All hits radial position √(x²+y²) below this value
xr = 1320.0
@printf("  Radial cut: %.1f mm\n", xr)

# Track length cut (mm)
trklm = 50
@printf("  Track length cut: %d mm\n", trklm)

# Z cuts: Define fiducial volume excluding central region [zir, zil]
# Allow hits in ranges (-1900, 0) or (0, 1900) - excluding central dead region
zil = 0.0      # Inner left boundary  
zir = 0.0      # Inner right boundary
zol = -1900.0  # Outer left boundary
zor = 1900.0   # Outer right boundary
@printf("  Z cuts: inner left=%.1f, inner right=%.1f, outer left=%.1f, outer right=%.1f mm\n", 
        zil, zir, zol, zor)

# Energy windows (keV)
energy_window_low = 2300.0  # keV 
energy_window_up = 2700.0
@printf("  Energy window: (%.1f, %.1f) keV\n", energy_window_low, energy_window_up)

# Energy resolution (keV): 
erex = 12.5  # keV
@printf("  Energy resolution: %.1f keV FWHM\n", erex)

# ROI bounds (keV): 
roi_xi = 2440.0
roi_xu = 2480.0
@printf("  ROI bounds: (%.1f, %.1f) keV\n", roi_xi, roi_xu)

# Data sets
dfstats = Dict(
    "0nubb" => 1e+4,
    "bi214_copper_endcaps" => 5e+6,
    "bi214_copper_shell" => 5e+6,
    "bi214_cathode_volume" => 1e+6,
    "bi214_ptfe_barrel" => 5e+6,
    "bi214_ptfe_endcap" => 5e+6,
    "tl208_copper_shell" => 5e+5,
    "tl208_copper_endcaps" => 5e+5,
    "tl208_cathode_volume" => 1e+5,
    "tl208_ptfe_barrel" => 1e+5,
    "tl208_ptfe_endcap" => 1e+5,
    "electron_0nubb_energy" => 1e+4
)

background_files = [
    "bi214_copper_shell",
    "bi214_copper_endcaps",
    "bi214_cathode_volume",
    "bi214_cathode_surface",
    "bi214_ptfe_barrel",
    "bi214_ptfe_endcap",
    "tl208_copper_shell",
    "tl208_copper_endcaps",
    "tl208_cathode_volume",
    "tl208_ptfe_barrel",
    "tl208_ptfe_endcap"
]

@printf("\nDataset Information:\n")
@printf("| %-25s | %10s |\n", "File", "Size")
@printf("|%s|%s|\n", "-"^27, "-"^12)
@printf("| %-25s | %10.0f |\n", "0nubb.h5", dfstats["0nubb"])
for file in background_files
    if haskey(dfstats, file)
        @printf("| %-25s | %10.0f |\n", file * ".h5", dfstats[file])
    end
end

# Event selection function
function event_selection(analysis_results, xr, trklm, 
                         zil, zir, zol, zor,
                         energy_window_low, energy_window_up, 
                         roi_xi, roi_xu, erex; 
                         signal_2e = 0.8, 
                         path_signal = "AnalysisSummary/znubb.js",
                         signal_file = "znubb")
    
    @printf("\n--- Event Selection for %s ---\n", signal_file)
    
    # Check if we have any single track events to begin with
    if analysis_results.n_single_track == 0
        @printf("No single track events available for analysis\n")
        signal_teff = 0.0
        json_analysis_summary(path_signal, signal_file, trklm, xr,
                              zil, zir, zol, zor,
                              erex, roi_xi, roi_xu,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, signal_2e, signal_teff)
        @printf("  Summary saved to: %s\n", path_signal)
        return signal_teff
    end
    
    # Apply radial cut
    filtered_st = filter_radial(analysis_results.single_track, xr, xr)
    n_filtered_st = number_of_events(filtered_st)
    
    # Check if we have any events left after radial cut
    if n_filtered_st > 0
        @printf("Single Track Events: mean energy = %.2f keV\n", mean(filtered_st.energy))
        @printf("  Radial cut > %.1f mm: efficiency = %.4f\n", 
                xr, n_filtered_st/analysis_results.n_single_track)
    else
        @printf("Single Track Events: No events after radial cut\n")
        @printf("  Radial cut > %.1f mm: efficiency = %.4f\n", 
                xr, n_filtered_st/analysis_results.n_single_track)
        
        # Return early with zero efficiencies if no events survive
        signal_teff = 0.0
        json_analysis_summary(path_signal, signal_file, trklm, xr,
                              zil, zir, zol, zor,
                              erex, roi_xi, roi_xu,
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, signal_2e, signal_teff)
        @printf("  Summary saved to: %s\n", path_signal)
        return signal_teff
    end

    # Apply track length cut
    filtered_xs = filter_short_tracks(filtered_st, trklm)
    n_filtered_xs = number_of_events(filtered_xs)
    
    @printf("  Track length > %d mm: efficiency = %.4f\n", 
            trklm, n_filtered_xs/n_filtered_st)

    # Apply Z cuts
    filtered_xz = filter_z(filtered_xs, zil, zir, zol, zor)
    n_filtered_xz = number_of_events(filtered_xz)
    
    @printf("  Z fiducial cuts: efficiency = %.4f\n", 
            n_filtered_xz/n_filtered_xs)

    # Energy Resolution Analysis
    @printf("Energy Resolution Analysis:\n")
    
    energy_df = combine(groupby(filtered_xz, :event_id), :energy => first => :energy)
    eht1, _ = step_hist(energy_df.energy;
         nbins = 40,
         xlim = (energy_window_low, energy_window_up),
         xlabel = "E (keV)",
         ylabel = "Frequency",
         title = "Single Track Energy Distribution")

    eres = smear_histogram(eht1, erex)
    ehrx, _ = step_hist(eres;
         nbins = 80,
         xlim = (energy_window_low, energy_window_up),
         xlabel = "E (keV)",
         ylabel = "Frequency",
         title = "Energy with $(erex) keV resolution")

    counts_roi = counts_in_range(ehrx, roi_xi, roi_xu)

    # Fix to take into account low statistics that may result in zero counts.
    # take an upper limit of 0.5 if 0 counts found. 

    if counts_roi == 0.0
        counts_roi = 0.5
    end
    
    eff_roi = counts_roi/n_filtered_xz
    
    @printf("  Counts in ROI [%.1f, %.1f] keV: %d\n", roi_xi, roi_xu, counts_roi)
    @printf("  ROI efficiency: %.4f\n", eff_roi)

    # Calculate all efficiencies
    signal_contained = analysis_results.n_events_processed/dfstats[signal_file]
    signal_1trk = analysis_results.n_single_track/analysis_results.n_events_processed
    signal_radial = n_filtered_st/analysis_results.n_single_track
    signal_trkl = n_filtered_xs/n_filtered_st
    signal_zfid = n_filtered_xz/n_filtered_xs
    signal_roi = eff_roi
    
    signal_teff = signal_contained * signal_1trk * signal_radial * signal_trkl * signal_zfid * signal_roi * signal_2e 

    @printf("\nSignal Efficiency Summary (%.1f keV FWHM resolution):\n", erex)
    @printf("  Fraction of events contained: %.4f\n", signal_contained)
    @printf("  Fraction of events 1 track: %.4f\n", signal_1trk)
    @printf("  Fraction of events radial cut: %.4f\n", signal_radial)
    @printf("  Fraction of events track length cut: %.4f\n", signal_trkl)
    @printf("  Fraction of events Z fiducial cut: %.4f\n", signal_zfid)
    @printf("  Fraction of events ROI: %.4f\n", signal_roi)
    @printf("  Fraction of events 2 electron ID: %.4f\n", signal_2e)
    @printf("  ** Total signal efficiency: %.6f **\n", signal_teff)
    
    # Save to JSON
    json_analysis_summary(path_signal, signal_file, trklm, xr,
                          zil, zir, zol, zor,
                          erex, roi_xi, roi_xu,
                          signal_contained,
                          signal_1trk,
                          signal_radial,
                          signal_trkl,
                          signal_zfid,
                          signal_roi,
                          signal_2e,
                          signal_teff)
    
    @printf("  Summary saved to: %s\n", path_signal)
    
    return signal_teff
end

# Main execution
function main()
    # Create AnalysisSummary directory if it doesn't exist (in parent directory)
    analysis_summary_dir = joinpath(dirname(@__DIR__), "AnalysisSummary")
    if !isdir(analysis_summary_dir)
        @printf("Creating output directory: %s\n", analysis_summary_dir)
        mkpath(analysis_summary_dir)
    end
    
    @printf("\n=== Processing Signal (0νββ) ===\n")
    @printf("Loading results from: %s\n", results_znubb)
    
    if !isdir(results_znubb)
        @printf("ERROR: Results directory not found: %s\n", results_znubb)
        return
    end
    
    rznubb = load_analysis_results(results_znubb)
    
    @printf("0νββ Results Loaded:\n")
    @printf("  Events processed: %d\n", rznubb.n_events_processed)
    @printf("  Single track events: %d\n", rznubb.n_single_track)
    @printf("  Two track events: %d\n", rznubb.n_two_track)
    @printf("  Three+ track events: %d\n", rznubb.n_three_plus_track)
    @printf("  Failed events: %d\n", rznubb.n_failed)

    @printf("\nCalling event selection for signal...\n")
    
    signal_efficiency = event_selection(rznubb, xr, trklm, 
                         zil, zir, zol, zor,
                         energy_window_low, energy_window_up, 
                         roi_xi, roi_xu, erex; 
                         signal_2e = 0.8, 
                         path_signal = joinpath(analysis_summary_dir, "znubb.js"),
                         signal_file = "0nubb")

    @printf("\n=== Processing Background Files ===\n")
    
    # Loop over all background files
    for background_file in background_files
        # Construct the results directory path
        results_bkg = joinpath(results_base_dir, background_file)
        # Extract display name from the file name
        bkg_display_name = replace(background_file, "_" => " ") |> titlecase

        if isdir(results_bkg)
            rbi214 = load_analysis_results(results_bkg)
        
            @printf("\n%s Results Loaded:\n", bkg_display_name)
            @printf("  Events processed: %d\n", rbi214.n_events_processed)
            @printf("  Single track events: %d\n", rbi214.n_single_track)
            @printf("  Two track events: %d\n", rbi214.n_two_track)
            @printf("  Three+ track events: %d\n", rbi214.n_three_plus_track)
            @printf("  Failed events: %d\n", rbi214.n_failed)
            
            @printf("\nEvent Classification Summary (%s):\n", bkg_display_name)
            @printf("| %-15s | %8s | %10s |\n", "Event Type", "Count", "Percentage")
            @printf("|%s|%s|%s|\n", "-"^17, "-"^10, "-"^12)
            @printf("| %-15s | %8d | %9.1f%% |\n", "Single Track", rbi214.n_single_track, 
                    100*rbi214.n_single_track/rbi214.n_events_processed)
            @printf("| %-15s | %8d | %9.1f%% |\n", "Two Tracks", rbi214.n_two_track, 
                    100*rbi214.n_two_track/rbi214.n_events_processed)
            @printf("| %-15s | %8d | %9.1f%% |\n", "Three+ Tracks", rbi214.n_three_plus_track, 
                    100*rbi214.n_three_plus_track/rbi214.n_events_processed)
            @printf("| %-15s | %8d | %9.1f%% |\n", "Failed", rbi214.n_failed, 
                    100*rbi214.n_failed/rbi214.n_events_processed)
            @printf("| %-15s | %8d | %9.1f%% |\n", "**Total**", rbi214.n_events_processed, 100.0)

            # Analysis
            path_bkg = joinpath(analysis_summary_dir, background_file * ".js")

            event_selection(rbi214, xr, trklm, 
                           zil, zir, zol, zor,
                           energy_window_low, energy_window_up, 
                           roi_xi, roi_xu, erex; 
                           signal_2e = 0.04, 
                           path_signal = path_bkg,
                           signal_file = background_file)
        else
            @printf("\nWARNING: Results directory not found: %s\n", results_bkg)
            @printf("Please ensure the analysis has been run for %s.\n", bkg_display_name)
        end
    end
    
    @printf("\n=== Analysis Complete ===\n")
    @printf("Check AnalysisSummary/ directory for JSON output files.\n")
end

# Run the script if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end