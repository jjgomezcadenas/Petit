module Petit

using DataFrames
using Graphs
using StatsBase
using Plots
using HDF5
using ProgressMeter
using Serialization

include("histos.jl")
using .histos

include("util.jl")
include("event_and_hits.jl")
include("voxels.jl")
include("tracks.jl")
include("hdt5_analysis.jl")
include("plots.jl")


export Tracks, AnalysisResults, TracksSummary
export get_event, voxelize_hits, euclidean_distance, build_tracks, select_events, analysis_loop, make_tracks
export voxel_distances, voxel_closest_distance, voxel_energy
export hits_per_event, hits_per_all_events, number_of_events, energy_primary, energy_deposited, find_events_with_alphas, filter_fiducial_events
export plot_hits_evt, plot_hits_trk, plot_hits
export event_loop, get_dataset_dfs, histogram_results, getdirs, smear_histogram, example_smearing, counts_in_range, save_analysis_results, load_analysis_results

# Re-export histos functions
export hist1d, hist2d, p1df, step_hist, get_histo1d, Histo1d
export save_histo1d, load_histo1d
export save_histos, load_histos

end
