
module Petit

using DataFrames
using CSV
using Graphs
using StatsBase
using Statistics
using SparseArrays
using NearestNeighbors
using Plots
using HDF5
using ProgressMeter
using Serialization
using JSON
using Dates
using Glob
using Interpolations
using KernelDensity
using Peaks
using LinearAlgebra

include("histos.jl")
using .histos

include("util.jl")
include("itaca_functions.jl")
include("event_and_hits.jl")
include("voxels.jl")
include("track_building.jl")
include("track_extreme_finding.jl")
include("track_analysis.jl")
include("track_blob_analysis.jl")
include("track_path_reconstruction.jl")
include("hdt5_analysis.jl")
include("hd5t_blob_analysis.jl")
include("hdt5_mt.jl")
include("plots.jl")
#include("migration_helper.jl")
include("json_io.jl")
include("track_io.jl")
include("fits.jl")
include("hd5st_functions.jl")
include("kde.jl")
include("mctrack.jl")
include("fom_and_roc.jl")
include("blob_functions.jl")


# Track building exports
export Tracks, DiffusionParams, AnalysisResults
export build_tracks, build_tracks_knn, build_tracks_kdtree, make_tracks

# Event and voxel exports
export get_event, voxelize_hits, voxelize_event, euclidean_distance, select_events, analysis_loop
export voxel_distances, voxel_closest_distance, voxel_energy
export hits_per_event, hits_per_all_events, number_of_events, energy_primary, energy_deposited
export find_events_with_alphas, filter_fiducial_events, filter_radial, filter_z, filter_short_tracks

# Track extreme finding exports
export find_track_extremes, find_track_extremes_internal, walk_track_from_extremes
export TrackCoords, extract_coords, find_path_bfs, calculate_path_length, calculate_path_length_from_coords
export euclidean_distance_coords, calculate_vertex_curvatures
export find_extremes_topology, find_extremes_curvature, find_extremes_spatial, find_extremes_combined
export find_extremes_energy_weighted, find_extremes_energy_weighted_relaxed
export find_extremes_edge_energy_weighted_opt, find_extremes_distance_fallback_opt
export find_extremes_mst_diameter, compute_mst_graph
export energy_weight_matrix, dijkstra_path

# Track analysis exports
export track_positions, track_energies_keV, track_stats
export walk_result_print, extreme_distances_print, path_print, blobs_print, diffusion_params_print
export kde_peaks_print
export diagnose_path_efficiency, diagnose_endpoint_degrees, diagnose_skeleton_coverage
export diagnose_endpoint_stability, diagnose_track_extent
# Markdown display functions (for Pluto notebooks)
export walk_result_md, path_md, diffusion_params_md, extreme_distances_md
export kde_peaks_md, blobs_md

# Track blob analysis exports
export find_blob_energies, energy_in_spheres_around_extremes, energy_blobs_from_path
export blob_asymmetry, blob_asymmetry_from_path, blob_analysis_vs_radius
export energy_in_variable_spheres_around_extremes

# Track path reconstruction exports
export reconstruct_central_path, reconstruct_end_point_energy
export get_raw_path, project_voxels_to_path, compute_extreme_distances

export plot_hits_evt, plot_event, plot_hits_trk, plot_hits, plot_track_with_extremes, plot_track_blobs, plot_track_walk, plot_paths, plot_track_with_paths
export event_loop, event_loop_single_track, analysis_loop_single_track, determine_events_to_process, validate_event_loop_parameters, load_and_prepare_data
export event_loop_single_track2, analysis_loop_single_track2
export event_loop_single_track_mt, analysis_loop_single_track_mt, get_optimal_threads, split_events_for_threads
export get_dataset_dfs, histogram_results, getdirs, smear_histogram, example_smearing, counts_in_range, save_analysis_results, load_analysis_results
export extract_unique_energies, get_track_energies, get_track_positions
export json_analysis_summary, read_json_analysis_summary
export save_track_csv, load_track_csv, save_track_binary, load_track_binary, save_track_with_analysis, print_track_summary
export nof_events, count_events
export save_tracks_to_hdf5, read_tracks_from_hdf5, chain_track_files, merge_csv_files
export read_reco_results_from_hdf5, chain_reco_results

# Export blob analysis functions
export Eff, TRK, NTRKS, Blobs, BlobsVariable
export get_bkgnd_tracks, get_bb0nu_tracks, get_xe137_tracks
export get_itaca_tracks, read_itaca_metadata, get_itaca_metadata_value
export match_itaca_events, match_itaca_events_df
export blob_analysis, blob_analysis_variable, write_blobs_csv, fom_blobs
export read_cnn_efficiencies, compute_roc, compute_all_rocs, plot_roc, plot_efficiency_vs_threshold

# Re-export histos functions
export hist1d, hist2d, p1df, step_hist, get_histo1d, Histo1d
export save_histo1d, load_histo1d
export save_histos, load_histos
export centers, edges

# Export fitting functions
export lfit, RFit, FGauss
export gpol1, gpol2, gpol3
export fit_pol1, fit_pol2, fit_pol3
export gausg, gausg2, gaussfm, gauss2fm
export fit_gauss, fit_gauss_fm, gfit_gauss2_cmean, fit_2gauss_cmean
export fitg1, fitg2, plot_fit_gauss
export func1dfit, fit_profile

# KDE functions for track energy analysis
export compute_arc_length, project_voxel_to_path, project_voxels_to_path
export energy_weighted_kde, find_maxima, analyze_track_energy_profile
export get_reco_kde, get_mc_kde, find_peaks, plot_kde_peaks

# MC track functions
export get_mc_extremes, compute_mc_path

# Diffusion functions
export diffuse_xyz_image_mc #, diffuse_xyz_image_kernel  # commented out - requires ImageFiltering

# FOM and ROC analysis functions
export efficiency_vs_cut, compute_efficiencies, figure_of_merit, roc_curve
export signal_eff, background_eff, background_rejection
export plot_efficiency_vs_cut, plot_fom_vs_cut, plot_roc, plot_all_fom_analysis
export print_fom_summary

# ITACA helper functions (from itaca_functions.jl)
export get_sigma, get_energy_threshold, get_voxel_size_and_distance
export length_and_energy_of_tracks, kde_peaks
export print_diagnostics

# Blob analysis plotting functions (from blob_functions.jl)
export apply_blob_cut, select_blob
export plot_eb2_cut_eff_and_fom, plot_eb1_vs_eb2
export plot_track_length, plot_asymmetry, plot_d1_vs_d2

# Utility functions
export load_csv_results

end
