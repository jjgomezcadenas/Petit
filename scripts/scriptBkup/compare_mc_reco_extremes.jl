#!/usr/bin/env julia

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using Printf
using DataFrames
using Statistics
using Plots
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

function get_extremes_positions(walk_result)
    ext1, ext2 = walk_result.extremes
    if isnothing(ext1)
        return nothing, nothing
    end
    pos1 = (ext1.x, ext1.y, ext1.z)
    pos2 = (ext2.x, ext2.y, ext2.z)
    return pos1, pos2
end

function distance_3d(p1, p2)
    return sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2 + (p1[3]-p2[3])^2)
end

function match_extremes(mc_pos1, mc_pos2, reco_pos1, reco_pos2)
    # Try both pairings and choose the one with minimum total distance
    d11 = distance_3d(mc_pos1, reco_pos1)
    d22 = distance_3d(mc_pos2, reco_pos2)
    pairing1_total = d11 + d22

    d12 = distance_3d(mc_pos1, reco_pos2)
    d21 = distance_3d(mc_pos2, reco_pos1)
    pairing2_total = d12 + d21

    if pairing1_total <= pairing2_total
        return (d1=d11, d2=d22, pairing=1)
    else
        return (d1=d12, d2=d21, pairing=2)
    end
end

function process_event(ievent::Int, bbdf::DataFrame,
                       mc_voxel_size::Float64, mc_max_distance::Float64,
                       reco_σt::Float64, reco_voxel::Float64, reco_max_distance::Float64,
                       energy_threshold_keV::Float64, reco_dfpars::Petit.DiffusionParams)

    bbevt = Petit.get_event(bbdf, ievent)
    bbevtmc = transform_hits_df(bbevt)

    # 1) Monte Carlo track: voxelize without diffusion
    mc_vx = Petit.voxelize_event(bbevtmc, mc_voxel_size)
    mc_dfpars = Petit.DiffusionParams()  # default params for MC

    mc_tracks = Petit.make_tracks(mc_vx;
                                  max_distance_mm=mc_max_distance,
                                  energy_threshold_kev=energy_threshold_keV,
                                  diffusion=mc_dfpars)

    if length(mc_tracks) != 1
        return (success=false, reason="mc_multi_track", n_mc_tracks=length(mc_tracks))
    end

    mc_track = mc_tracks[1]
    mc_walk = Petit.walk_track_from_extremes(mc_track)
    mc_pos1, mc_pos2 = get_extremes_positions(mc_walk)

    if isnothing(mc_pos1)
        return (success=false, reason="mc_no_extremes")
    end

    # 2) Reconstructed track: diffuse, voxelize, find extremes
    reco_df = Petit.diffuse_xyz_image_mc(bbevtmc;
                                         sigma_t_mm=reco_σt,
                                         sigma_l_mm=σl_ion,
                                         nbins=nbins, nsigma=nsigma)

    reco_vx = Petit.voxelize_event(reco_df, reco_voxel)

    reco_tracks = Petit.make_tracks(reco_vx;
                                    max_distance_mm=reco_max_distance,
                                    energy_threshold_kev=energy_threshold_keV,
                                    diffusion=reco_dfpars)

    if length(reco_tracks) != 1
        return (success=false, reason="reco_multi_track", n_reco_tracks=length(reco_tracks))
    end

    reco_track = reco_tracks[1]
    reco_walk = Petit.walk_track_from_extremes(reco_track)
    reco_pos1, reco_pos2 = get_extremes_positions(reco_walk)

    if isnothing(reco_pos1)
        return (success=false, reason="reco_no_extremes")
    end

    # 3) Match and compute distances
    match = match_extremes(mc_pos1, mc_pos2, reco_pos1, reco_pos2)

    # 4) Collect diagnostic info
    mc_graph = mc_track.graph
    reco_graph = reco_track.graph

    mc_degrees = [Graphs.degree(mc_graph, v) for v in Graphs.vertices(mc_graph)]
    reco_degrees = [Graphs.degree(reco_graph, v) for v in Graphs.vertices(reco_graph)]

    mc_deg1_count = count(d -> d == 1, mc_degrees)
    reco_deg1_count = count(d -> d == 1, reco_degrees)

    # Track extent
    mc_extent = sqrt((extrema(mc_track.voxels.x)[2] - extrema(mc_track.voxels.x)[1])^2 +
                     (extrema(mc_track.voxels.y)[2] - extrema(mc_track.voxels.y)[1])^2 +
                     (extrema(mc_track.voxels.z)[2] - extrema(mc_track.voxels.z)[1])^2)

    return (success=true,
            mc_pos1=mc_pos1, mc_pos2=mc_pos2,
            reco_pos1=reco_pos1, reco_pos2=reco_pos2,
            d1=match.d1, d2=match.d2,
            pairing=match.pairing,
            mc_nvox=nrow(mc_track.voxels), reco_nvox=nrow(reco_track.voxels),
            mc_deg1=mc_deg1_count, reco_deg1=reco_deg1_count,
            mc_avg_deg=mean(mc_degrees), reco_avg_deg=mean(reco_degrees),
            mc_confidence=mc_walk.confidence, reco_confidence=reco_walk.confidence,
            mc_path_len=mc_walk.total_length, reco_path_len=reco_walk.total_length,
            mc_extent=mc_extent)
end

function main(; ievent::Int=1, levent::Int=100, ldrft::Float64=100.0,
               voxel_scale::Float64=1.0, voxel_distance_scale::Float64=3.0,
               confidence_threshold::Float64=0.0,
               input_file::String="bb0nu/bb0nu_15bar_p1.h5")

    # MC parameters (no diffusion, 1mm voxels)
    mc_voxel_size = 1.0  # mm
    mc_max_distance = 3.0  # mm

    # Reco parameters (with diffusion)
    reco_σt = Petit.sigma_t_ion_mm(tK, ldrft, edrift)
    reco_voxel = reco_σt * voxel_scale
    reco_max_distance = reco_voxel * voxel_distance_scale

    energy_threshold_keV = energy_threshold_ions / fkeV
    reco_dfpars = Petit.DiffusionParams(ldrft, reco_σt, σl_ion,
                                        reco_voxel, reco_max_distance, energy_threshold_keV,
                                        nbins, nsigma)

    bbdf = load_data(input_file)

    d1_values = Float64[]
    d2_values = Float64[]
    conf_values = Float64[]
    event_ids = Int[]
    results_all = []

    n_success = 0
    n_mc_multi = 0
    n_reco_multi = 0
    n_low_confidence = 0

    for evt in ievent:levent
        result = process_event(evt, bbdf, mc_voxel_size, mc_max_distance,
                               reco_σt, reco_voxel, reco_max_distance,
                               energy_threshold_keV, reco_dfpars)

        push!(results_all, (event=evt, result=result))

        if result.success
            # Apply confidence threshold
            if confidence_threshold > 0.0 && result.reco_confidence < confidence_threshold
                n_low_confidence += 1
                continue
            end

            n_success += 1
            push!(d1_values, result.d1)
            push!(d2_values, result.d2)
            push!(conf_values, result.reco_confidence)
            push!(event_ids, evt)
        else
            if result.reason == "mc_multi_track"
                n_mc_multi += 1
            elseif result.reason == "reco_multi_track"
                n_reco_multi += 1
            end
        end
    end

    n_total = levent - ievent + 1
    n_single_track_reco = n_success + n_low_confidence  # events that passed single-track before conf cut

    println("=" ^ 70)
    println("MC vs RECO EXTREMES COMPARISON")
    println("=" ^ 70)

    # Efficiency table
    println("\n--- EFFICIENCY TABLE ---")
    println("Total events:                  $n_total")
    println("MC single-track:               $(n_total - n_mc_multi) ($(round(100*(n_total-n_mc_multi)/n_total, digits=1))%)")
    println("Reco single-track:             $n_single_track_reco ($(round(100*n_single_track_reco/n_total, digits=1))%)")
    if confidence_threshold > 0.0
        println("Confidence threshold:          $confidence_threshold")
        println("Rejected (low conf):           $n_low_confidence")
        println("Final reconstructed:           $n_success ($(round(100*n_success/n_total, digits=1))%)")
    else
        println("Final reconstructed:           $n_success ($(round(100*n_success/n_total, digits=1))%)")
    end

    println("\n--- RECONSTRUCTION QUALITY ---")
    println("d1: $(round(mean(d1_values), digits=2)) ± $(round(std(d1_values), digits=2)) mm")
    println("d2: $(round(mean(d2_values), digits=2)) ± $(round(std(d2_values), digits=2)) mm")

    if n_success > 0
        # Find outliers (distance > 10 mm)
        outlier_threshold = 10.0
        outlier_indices = findall(i -> d1_values[i] > outlier_threshold || d2_values[i] > outlier_threshold,
                                   1:length(d1_values))
        outlier_events = event_ids[outlier_indices]

        println("\nOutliers (d > $(outlier_threshold) mm): $(length(outlier_events)) events")
        println("Outlier event IDs: $outlier_events")

        # Histogram 1: d1
        p1 = histogram(d1_values,
                       bins=range(0, min(60, maximum(d1_values)*1.1), length=25),
                       xlabel="d1 (mm)", ylabel="Counts",
                       title="Extreme 1 Error (mean=$(round(mean(d1_values), digits=1)) mm)",
                       label="", fillalpha=0.7, color=:blue)

        # Histogram 2: d2
        p2 = histogram(d2_values,
                       bins=range(0, min(60, maximum(d2_values)*1.1), length=25),
                       xlabel="d2 (mm)", ylabel="Counts",
                       title="Extreme 2 Error (mean=$(round(mean(d2_values), digits=1)) mm)",
                       label="", fillalpha=0.7, color=:red)

        # Scatter: d1 vs d2
        p3 = scatter(d1_values, d2_values,
                     xlabel="d1 (mm)", ylabel="d2 (mm)",
                     title="d1 vs d2",
                     label="", alpha=0.6, markersize=4)
        plot!(p3, [0, 60], [0, 60], label="d1=d2", linestyle=:dash, color=:gray)

        # Combined layout
        p = plot(p1, p2, p3, layout=(1, 3), size=(1200, 350))
        display(p)
        savefig(p, joinpath(pdir, "scripts", "mc_reco_extremes_comparison.png"))

        # Analyze outliers vs good events
        good_indices = findall(i -> d1_values[i] <= outlier_threshold && d2_values[i] <= outlier_threshold,
                               1:length(d1_values))

        println("\n" * "=" ^ 70)
        println("OUTLIER ANALYSIS")
        println("=" ^ 70)

        # Collect stats for good vs outlier events
        good_reco_deg1 = Float64[]
        good_reco_avg_deg = Float64[]
        good_reco_conf = Float64[]
        good_reco_nvox = Float64[]

        outlier_reco_deg1 = Float64[]
        outlier_reco_avg_deg = Float64[]
        outlier_reco_conf = Float64[]
        outlier_reco_nvox = Float64[]

        for (i, evt) in enumerate(event_ids)
            r = results_all[evt].result
            if i in good_indices
                push!(good_reco_deg1, r.reco_deg1)
                push!(good_reco_avg_deg, r.reco_avg_deg)
                push!(good_reco_conf, r.reco_confidence)
                push!(good_reco_nvox, r.reco_nvox)
            else
                push!(outlier_reco_deg1, r.reco_deg1)
                push!(outlier_reco_avg_deg, r.reco_avg_deg)
                push!(outlier_reco_conf, r.reco_confidence)
                push!(outlier_reco_nvox, r.reco_nvox)
            end
        end

        println("\nGood events (n=$(length(good_indices))):")
        println("  Reco deg-1 vertices: $(round(mean(good_reco_deg1), digits=1)) ± $(round(std(good_reco_deg1), digits=1))")
        println("  Reco avg degree:     $(round(mean(good_reco_avg_deg), digits=2)) ± $(round(std(good_reco_avg_deg), digits=2))")
        println("  Reco confidence:     $(round(mean(good_reco_conf), digits=2)) ± $(round(std(good_reco_conf), digits=2))")
        println("  Reco n_voxels:       $(round(mean(good_reco_nvox), digits=0)) ± $(round(std(good_reco_nvox), digits=0))")

        println("\nOutlier events (n=$(length(outlier_indices))):")
        println("  Reco deg-1 vertices: $(round(mean(outlier_reco_deg1), digits=1)) ± $(round(std(outlier_reco_deg1), digits=1))")
        println("  Reco avg degree:     $(round(mean(outlier_reco_avg_deg), digits=2)) ± $(round(std(outlier_reco_avg_deg), digits=2))")
        println("  Reco confidence:     $(round(mean(outlier_reco_conf), digits=2)) ± $(round(std(outlier_reco_conf), digits=2))")
        println("  Reco n_voxels:       $(round(mean(outlier_reco_nvox), digits=0)) ± $(round(std(outlier_reco_nvox), digits=0))")

        # Detailed look at worst outliers
        worst_indices = sortperm(max.(d1_values, d2_values), rev=true)[1:min(5, length(outlier_indices))]
        println("\n" * "-" ^ 70)
        println("WORST 5 OUTLIERS:")
        println("-" ^ 70)
        for idx in worst_indices
            evt = event_ids[idx]
            r = results_all[evt].result
            println("\nEvent $evt: d1=$(round(d1_values[idx], digits=1)) mm, d2=$(round(d2_values[idx], digits=1)) mm")
            println("  MC:   pos1=$(round.(r.mc_pos1, digits=1)), pos2=$(round.(r.mc_pos2, digits=1))")
            println("  Reco: pos1=$(round.(r.reco_pos1, digits=1)), pos2=$(round.(r.reco_pos2, digits=1))")
            println("  Reco graph: $(r.reco_nvox) voxels, $(r.reco_deg1) deg-1, avg_deg=$(round(r.reco_avg_deg, digits=2))")
            println("  Reco confidence: $(round(r.reco_confidence, digits=2)), path_len=$(round(r.reco_path_len, digits=1)) mm")
            println("  MC extent: $(round(r.mc_extent, digits=1)) mm")
        end

        # Summary and recommendations
        println("\n" * "=" ^ 70)
        println("ANALYSIS SUMMARY & RECOMMENDATIONS")
        println("=" ^ 70)

        # Check correlations
        has_few_deg1 = mean(outlier_reco_deg1) < mean(good_reco_deg1) - 0.5
        has_high_degree = mean(outlier_reco_avg_deg) > mean(good_reco_avg_deg) + 0.5
        has_low_conf = mean(outlier_reco_conf) < mean(good_reco_conf) - 0.05

        println("\nKey findings:")
        if has_few_deg1
            println("  - Outliers have FEWER degree-1 vertices ($(round(mean(outlier_reco_deg1), digits=1)) vs $(round(mean(good_reco_deg1), digits=1)))")
            println("    -> Tracks are more densely connected, true endpoints not isolated")
        end
        if has_high_degree
            println("  - Outliers have HIGHER average degree ($(round(mean(outlier_reco_avg_deg), digits=2)) vs $(round(mean(good_reco_avg_deg), digits=2)))")
            println("    -> Dense graphs make topology-based extreme finding unreliable")
        end
        if has_low_conf
            println("  - Outliers have LOWER confidence ($(round(mean(outlier_reco_conf), digits=2)) vs $(round(mean(good_reco_conf), digits=2)))")
            println("    -> Algorithm is uncertain, could use this to flag bad events")
        end

        println("\nRecommended improvements:")
        println("  1. Use energy-weighted extreme selection: Bragg peaks have high energy")
        println("  2. For dense tracks (avg_degree > 4), prioritize spatial extremes over topology")
        println("  3. Add confidence threshold: reject events with confidence < 0.7")
        println("  4. Consider path coverage: extremes should span most of track extent")
    end

    return (d1=d1_values, d2=d2_values, events=event_ids, results=results_all,
            outlier_events=outlier_events)
end

# Run
if abspath(PROGRAM_FILE) == @__FILE__
    # Electrons (single Bragg peak)
    println("\n" * "=" ^ 70)
    println("ELECTRONS (2400-2500 keV)")
    println("=" ^ 70)
    main(ievent=1, levent=200, voxel_scale=1.0, voxel_distance_scale=1.5,
         confidence_threshold=0.0,
         input_file="xe137/electrons_2400_2500_15bar_100mum.next.h5")
end
