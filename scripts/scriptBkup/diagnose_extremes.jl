#!/usr/bin/env julia

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using Printf
using DataFrames
using CSV
using Glob
using Statistics
using Graphs

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

cmdir = joinpath(ENV["DATA"], "HD5t/itaca")

function load_data(input_file)
    input_path = joinpath(cmdir, input_file)
    dfs = Petit.get_dataset_dfs(input_path)
    return dfs["hits"]
end

function transform_hits_df(df::DataFrame; energy_to_electrons::Float64=1e5/2.5)
    df2 = select(df, Not([:time, :label, :particle_id, :hit_id]))
    df2.electrons = round.(Int, df2.energy .* energy_to_electrons)
    return df2
end

# Constants
f = 1e+5/2.5
fkeV = f*1e-3
nbins = 100
nsigma = 3.0
tK = 297.0
edrift = 500.0
σl_ion = 0.0
energy_threshold_ions = 10.0

function analyze_track_graph(track::Petit.Tracks)
    g = track.graph
    n_verts = nv(g)
    n_edges = ne(g)

    # Degree distribution
    degrees = [degree(g, v) for v in vertices(g)]
    min_deg = minimum(degrees)
    max_deg = maximum(degrees)
    avg_deg = mean(degrees)

    # Count degree-1 vertices (true endpoints)
    deg1_verts = findall(d -> d == 1, degrees)
    n_deg1 = length(deg1_verts)

    # Connected components
    comps = connected_components(g)
    n_comps = length(comps)

    return (n_vertices=n_verts, n_edges=n_edges,
            min_degree=min_deg, max_degree=max_deg, avg_degree=avg_deg,
            n_degree1=n_deg1, degree1_vertices=deg1_verts,
            n_components=n_comps)
end

function analyze_extremes(track::Petit.Tracks, walk_result)
    vox = track.voxels
    n_vox = nrow(vox)

    # Get extremes info
    ext1, ext2 = walk_result.extremes
    path = walk_result.path_indices
    path_len = walk_result.total_length
    confidence = walk_result.confidence

    if isnothing(ext1)
        return (valid=false, n_voxels=n_vox)
    end

    # Positions of extremes
    pos1 = (ext1.x, ext1.y, ext1.z)
    pos2 = (ext2.x, ext2.y, ext2.z)

    # Straight-line distance between extremes
    straight_dist = sqrt((pos2[1]-pos1[1])^2 + (pos2[2]-pos1[2])^2 + (pos2[3]-pos1[3])^2)

    # Path efficiency (straight/path)
    efficiency = path_len > 0 ? straight_dist / path_len : 0.0

    # Energy at extremes
    e1 = ext1.energy
    e2 = ext2.energy

    # Find spatial extremes of track (bounding box corners)
    x_min, x_max = extrema(vox.x)
    y_min, y_max = extrema(vox.y)
    z_min, z_max = extrema(vox.z)

    # Check if found extremes are near spatial extremes
    is_ext1_spatial = (ext1.x == x_min || ext1.x == x_max ||
                       ext1.y == y_min || ext1.y == y_max ||
                       ext1.z == z_min || ext1.z == z_max)
    is_ext2_spatial = (ext2.x == x_min || ext2.x == x_max ||
                       ext2.y == y_min || ext2.y == y_max ||
                       ext2.z == z_min || ext2.z == z_max)

    # Check for high energy at extremes (Bragg peak indicator)
    mean_energy = mean(vox.energy)
    is_ext1_high_energy = e1 > 1.5 * mean_energy
    is_ext2_high_energy = e2 > 1.5 * mean_energy

    return (valid=true, n_voxels=n_vox,
            pos1=pos1, pos2=pos2,
            straight_dist=straight_dist, path_length=path_len,
            efficiency=efficiency, confidence=confidence,
            n_path_voxels=length(path),
            energy1=e1, energy2=e2,
            is_ext1_spatial=is_ext1_spatial, is_ext2_spatial=is_ext2_spatial,
            is_ext1_high_energy=is_ext1_high_energy, is_ext2_high_energy=is_ext2_high_energy)
end

function diagnose_event(ievent::Int, bbdf::DataFrame, ion_dfpars::Petit.DiffusionParams,
                        σt_ion::Float64, ion_voxel::Float64, max_distance::Float64,
                        energy_threshold_keV::Float64)

    bbevt = Petit.get_event(bbdf, ievent)
    bbevtmc = transform_hits_df(bbevt)

    # Original MC hits info
    mc_x_range = extrema(bbevtmc.x)
    mc_y_range = extrema(bbevtmc.y)
    mc_z_range = extrema(bbevtmc.z)
    mc_extent = sqrt((mc_x_range[2]-mc_x_range[1])^2 +
                     (mc_y_range[2]-mc_y_range[1])^2 +
                     (mc_z_range[2]-mc_z_range[1])^2)

    # Diffuse and voxelize
    bbion_df = Petit.diffuse_xyz_image_mc(bbevtmc;
                                          sigma_t_mm=σt_ion,
                                          sigma_l_mm=σl_ion,
                                          nbins=nbins, nsigma=nsigma)
    bbion_vx = Petit.voxelize_event(bbion_df, ion_voxel)

    # Make tracks
    tracks = Petit.make_tracks(bbion_vx;
                               max_distance_mm=max_distance,
                               energy_threshold_kev=energy_threshold_keV,
                               diffusion=ion_dfpars)

    n_tracks = length(tracks)

    if n_tracks != 1
        return (event=ievent, n_tracks=n_tracks, mc_extent=mc_extent,
                success=false, reason="multiple_tracks")
    end

    track = tracks[1]
    graph_info = analyze_track_graph(track)

    # Walk track to find extremes
    walk_result = Petit.walk_track_from_extremes(track)
    extremes_info = analyze_extremes(track, walk_result)

    if !extremes_info.valid
        return (event=ievent, n_tracks=n_tracks, mc_extent=mc_extent,
                success=false, reason="invalid_extremes", graph_info=graph_info)
    end

    # Quality assessment
    # Good extremes should:
    # 1. Have high path efficiency (> 0.5 for curved tracks)
    # 2. Cover most of the track extent
    # 3. Have at least one high-energy endpoint (Bragg peak)

    extent_coverage = extremes_info.straight_dist / mc_extent
    has_bragg_peak = extremes_info.is_ext1_high_energy || extremes_info.is_ext2_high_energy

    quality = "good"
    issues = String[]

    if extremes_info.efficiency < 0.3
        quality = "poor"
        push!(issues, "low_efficiency=$(round(extremes_info.efficiency, digits=2))")
    end

    if extent_coverage < 0.5
        quality = "poor"
        push!(issues, "low_coverage=$(round(extent_coverage, digits=2))")
    end

    if !has_bragg_peak
        push!(issues, "no_bragg_peak")
    end

    if extremes_info.confidence < 0.7
        push!(issues, "low_confidence=$(round(extremes_info.confidence, digits=2))")
    end

    return (event=ievent, n_tracks=n_tracks, mc_extent=mc_extent,
            success=true, quality=quality, issues=issues,
            graph_info=graph_info, extremes_info=extremes_info,
            extent_coverage=extent_coverage)
end

function main(; ievent::Int=1, levent::Int=10, ldrft::Float64=100.0,
               voxel_scale::Float64=1.0, voxel_distance_scale::Float64=3.0)

    σt_ion = Petit.sigma_t_ion_mm(tK, ldrft, edrift)
    ion_voxel = σt_ion * voxel_scale
    max_distance = ion_voxel * voxel_distance_scale

    energy_threshold_keV = energy_threshold_ions / fkeV
    ion_dfpars = Petit.DiffusionParams(ldrft, σt_ion, σl_ion,
                                       ion_voxel, max_distance, energy_threshold_keV,
                                       nbins, nsigma)

    println("=" ^ 80)
    println("TRACK EXTREMES DIAGNOSTIC ANALYSIS")
    println("=" ^ 80)
    println("Parameters:")
    println("  σ_t = $(round(σt_ion, digits=3)) mm")
    println("  voxel_size = $(round(ion_voxel, digits=3)) mm")
    println("  max_distance = $(round(max_distance, digits=3)) mm")
    println("  Events: $ievent to $levent")
    println()

    bbdf = load_data("bb0nu/bb0nu_15bar_p1.h5")

    results = []
    good_events = Int[]
    poor_events = Int[]
    skipped_events = Int[]

    for evt in ievent:levent
        result = diagnose_event(evt, bbdf, ion_dfpars, σt_ion, ion_voxel,
                                max_distance, energy_threshold_keV)
        push!(results, result)

        if !result.success
            push!(skipped_events, evt)
        elseif result.quality == "good"
            push!(good_events, evt)
        else
            push!(poor_events, evt)
        end
    end

    # Print summary
    println("=" ^ 80)
    println("SUMMARY")
    println("=" ^ 80)
    println("Total events: $(levent - ievent + 1)")
    println("Good quality: $(length(good_events)) - events: $good_events")
    println("Poor quality: $(length(poor_events)) - events: $poor_events")
    println("Skipped (multi-track): $(length(skipped_events)) - events: $skipped_events")
    println()

    # Detailed analysis of poor quality events
    if !isempty(poor_events)
        println("=" ^ 80)
        println("POOR QUALITY EVENT DETAILS")
        println("=" ^ 80)
        for evt in poor_events
            r = results[evt - ievent + 1]
            println("\nEvent $evt:")
            println("  Issues: $(join(r.issues, ", "))")
            println("  Graph: $(r.graph_info.n_vertices) voxels, $(r.graph_info.n_edges) edges")
            println("  Degree-1 vertices: $(r.graph_info.n_degree1)")
            println("  Avg degree: $(round(r.graph_info.avg_degree, digits=2))")
            println("  Path: $(r.extremes_info.n_path_voxels) voxels, $(round(r.extremes_info.path_length, digits=1)) mm")
            println("  Straight dist: $(round(r.extremes_info.straight_dist, digits=1)) mm")
            println("  MC extent: $(round(r.mc_extent, digits=1)) mm")
            println("  Efficiency: $(round(r.extremes_info.efficiency, digits=3))")
            println("  Confidence: $(round(r.extremes_info.confidence, digits=3))")
            println("  Ext1 energy: $(round(r.extremes_info.energy1*1e3, digits=1)) keV, high=$(r.extremes_info.is_ext1_high_energy)")
            println("  Ext2 energy: $(round(r.extremes_info.energy2*1e3, digits=1)) keV, high=$(r.extremes_info.is_ext2_high_energy)")
        end
    end

    # Detailed analysis of good quality events for comparison
    if !isempty(good_events) && length(good_events) <= 5
        println("\n" * "=" ^ 80)
        println("GOOD QUALITY EVENT DETAILS (for comparison)")
        println("=" ^ 80)
        for evt in good_events
            r = results[evt - ievent + 1]
            println("\nEvent $evt:")
            println("  Graph: $(r.graph_info.n_vertices) voxels, $(r.graph_info.n_edges) edges")
            println("  Degree-1 vertices: $(r.graph_info.n_degree1)")
            println("  Avg degree: $(round(r.graph_info.avg_degree, digits=2))")
            println("  Path: $(r.extremes_info.n_path_voxels) voxels, $(round(r.extremes_info.path_length, digits=1)) mm")
            println("  Straight dist: $(round(r.extremes_info.straight_dist, digits=1)) mm")
            println("  MC extent: $(round(r.mc_extent, digits=1)) mm")
            println("  Efficiency: $(round(r.extremes_info.efficiency, digits=3))")
            println("  Confidence: $(round(r.extremes_info.confidence, digits=3))")
            println("  Ext1 energy: $(round(r.extremes_info.energy1*1e3, digits=1)) keV, high=$(r.extremes_info.is_ext1_high_energy)")
            println("  Ext2 energy: $(round(r.extremes_info.energy2*1e3, digits=1)) keV, high=$(r.extremes_info.is_ext2_high_energy)")
        end
    end

    return results
end

# Run
if abspath(PROGRAM_FILE) == @__FILE__
    main(ievent=1, levent=20, voxel_distance_scale=3.0)
end
