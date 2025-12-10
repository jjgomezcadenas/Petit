
module Petit

using DataFrames
using CSV
using Graphs
using StatsBase
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

include("histos.jl")
using .histos

include("util.jl")
include("itaca_functions.jl")
include("event_and_hits.jl")
include("voxels.jl")
include("tracks.jl")
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


export Tracks, DiffusionParams, AnalysisResults
export get_event, voxelize_hits, voxelize_event, euclidean_distance, build_tracks, select_events, analysis_loop, make_tracks
export voxel_distances, voxel_closest_distance, voxel_energy
export hits_per_event, hits_per_all_events, number_of_events, energy_primary, energy_deposited, find_events_with_alphas, filter_fiducial_events, filter_radial, filter_z, filter_short_tracks
export find_track_extremes, find_track_extremes_legacy, find_track_extremes_improved, find_extremes_topology_based, find_extremes_curvature_based, find_extremes_combined, find_extremes_spatial_based
export get_raw_path, compute_extreme_distances
export calculate_vertex_curvatures, calculate_path_length, are_vertices_too_close, find_extremes_distance_fallback
export find_path_bfs, walk_track_from_extremes, energy_in_spheres_around_extremes, energy_in_variable_spheres_around_extremes, energy_blobs_from_path, blob_asymmetry, blob_asymmetry_from_path, blob_analysis_vs_radius
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
export track_energies_keV, blob_analysis, blob_analysis_variable, write_blobs_csv, fom_blobs
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

# FOM and ROC analysis functions
export efficiency_vs_cut, compute_efficiencies, figure_of_merit, roc_curve
export plot_efficiency_vs_cut, plot_fom_vs_cut, plot_roc, plot_all_fom_analysis
export print_fom_summary

end
