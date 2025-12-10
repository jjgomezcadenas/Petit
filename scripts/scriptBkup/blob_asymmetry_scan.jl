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

const cmdir = joinpath(ENV["DATA"], "HD5t/itaca")

# Constants
const f = 1e+5/2.5
const fkeV = f*1e-3
const nbins = 100
const nsigma = 3.0
const tK = 297.0
const edrift = 500.0
const σl_ion = 0.0
const energy_threshold_ions = 10.0

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

"""
    compute_asymmetry(E1, E2)

Compute blob energy asymmetry: (E1 - E2) / (E1 + E2)
where E1 >= E2.
"""
function compute_asymmetry(E1::Float64, E2::Float64)
    if E1 + E2 <= 0
        return NaN
    end
    return (E1 - E2) / (E1 + E2)
end

"""
    process_event(ievent, bbdf, σt, voxel_size, max_distance, energy_threshold_keV, dfpars)

Process a single event: voxelize, make tracks, find extremes using central path reconstruction.
Returns (track, walk_result, central_path) if single-track event, nothing otherwise.
"""
function process_event(ievent::Int, bbdf::DataFrame,
                       σt::Float64, voxel_size::Float64, max_distance::Float64,
                       energy_threshold_keV::Float64, dfpars::Petit.DiffusionParams)

    bbevt = Petit.get_event(bbdf, ievent)
    bbevtmc = transform_hits_df(bbevt)

    # Diffuse and voxelize
    reco_df = Petit.diffuse_xyz_image_mc(bbevtmc;
                                         sigma_t_mm=σt,
                                         sigma_l_mm=σl_ion,
                                         nbins=nbins, nsigma=nsigma)

    reco_vx = Petit.voxelize_event(reco_df, voxel_size)

    tracks = Petit.make_tracks(reco_vx;
                               max_distance_mm=max_distance,
                               energy_threshold_kev=energy_threshold_keV,
                               diffusion=dfpars)

    if length(tracks) != 1
        return nothing
    end

    track = tracks[1]

    # Walk track to get extremes and path
    walk_result = Petit.walk_track_from_extremes(track)

    if isnothing(walk_result.extremes[1]) || isempty(walk_result.path_indices)
        return nothing
    end

    # Reconstruct central path (smoothed) - same as test_blob_analysis.jl
    central_path = Petit.reconstruct_central_path(track,
                                                   walk_result.path_indices;
                                                   filter_radius = nsigma * σt)

    if nrow(central_path) < 2
        return nothing
    end

    return (track=track, walk_result=walk_result, central_path=central_path)
end

function run_asymmetry_scan(; ievent::Int=1, levent::Int=200, ldrft::Float64=100.0,
                            voxel_scale::Float64=1.0, voxel_distance_scale::Float64=1.5,
                            radii::Vector{Float64}=[5.0, 10.0, 15.0],
                            input_file::String, label::String)

    # Reco parameters
    σt = Petit.sigma_t_ion_mm(tK, ldrft, edrift)
    voxel_size = σt * voxel_scale
    max_distance = voxel_size * voxel_distance_scale
    energy_threshold_keV = energy_threshold_ions / fkeV
    dfpars = Petit.DiffusionParams(ldrft, σt, σl_ion,
                                   voxel_size, max_distance, energy_threshold_keV,
                                   nbins, nsigma)

    println("Loading data from: $input_file")
    bbdf = load_data(input_file)

    # Storage for results per radius
    results = Dict{Float64, NamedTuple}()

    for r in radii
        results[r] = (E1=Float64[], E2=Float64[], asym=Float64[],
                      track_len=Float64[], confidence=Float64[], event_ids=Int[])
    end

    n_total = levent - ievent + 1
    n_single_track = 0

    println("Processing events $ievent to $levent...")
    for evt in ievent:levent
        result = process_event(evt, bbdf, σt, voxel_size, max_distance,
                               energy_threshold_keV, dfpars)

        if isnothing(result)
            continue
        end

        n_single_track += 1

        track = result.track
        walk_result = result.walk_result
        central_path = result.central_path

        track_len = walk_result.total_length
        confidence = walk_result.confidence

        # Compute blob energies for each radius using find_blob_energies (with central path)
        for r in radii
            blobs = Petit.find_blob_energies(track, central_path; radius=r)

            # Eb1 is already the higher energy blob, Eb2 is lower
            E1 = blobs.Eb1 * 1e3  # Convert to keV
            E2 = blobs.Eb2 * 1e3
            asym = blobs.asymmetry

            push!(results[r].E1, E1)
            push!(results[r].E2, E2)
            push!(results[r].asym, asym)
            push!(results[r].track_len, track_len)
            push!(results[r].confidence, confidence)
            push!(results[r].event_ids, evt)
        end

        if n_single_track % 50 == 0
            println("  Processed $n_single_track single-track events...")
        end
    end

    println("\n" * "=" ^ 70)
    println("BLOB ASYMMETRY SCAN (using central path): $label")
    println("=" ^ 70)
    println("Total events:        $n_total")
    println("Single-track events: $n_single_track ($(round(100*n_single_track/n_total, digits=1))%)")

    # Print table for each radius
    println("\n" * "-" ^ 70)
    println("ASYMMETRY STATISTICS BY RADIUS")
    println("-" ^ 70)
    @printf("%-8s %10s %10s %10s %10s %10s\n", "R (mm)", "E1 (keV)", "E2 (keV)", "Asym μ", "Asym σ", "N events")
    println("-" ^ 70)

    for r in radii
        n = length(results[r].asym)
        if n > 0
            valid_asym = filter(!isnan, results[r].asym)
            @printf("%-8.1f %10.1f %10.1f %10.3f %10.3f %10d\n",
                    r,
                    mean(results[r].E1),
                    mean(results[r].E2),
                    mean(valid_asym),
                    std(valid_asym),
                    n)
        end
    end

    return results, n_single_track
end

function main()
    radii = [20.0]

    # Run for bb0nu (double beta)
    println("\n" * "=" ^ 70)
    println("DOUBLE BETA (bb0nu) - Expected: LOW asymmetry (two Bragg peaks)")
    println("=" ^ 70)
    bb_results, bb_n = run_asymmetry_scan(
        ievent=1, levent=200,
        radii=radii,
        input_file="bb0nu/0nubb_15bar_100mum.next.h5",
        label="bb0nu"
    )

    # Run for single electrons
    println("\n" * "=" ^ 70)
    println("SINGLE ELECTRONS (xe137) - Expected: HIGH asymmetry (one Bragg peak)")
    println("=" ^ 70)
    ele_results, ele_n = run_asymmetry_scan(
        ievent=1, levent=200,
        radii=radii,
        input_file="xe137/electrons_2400_2500_15bar_100mum.next.h5",
        label="electrons"
    )

    # Create comparison histograms
    println("\n" * "=" ^ 70)
    println("GENERATING HISTOGRAMS")
    println("=" ^ 70)

    plots_array = []

    for r in radii
        # Get valid asymmetries
        bb_asym = filter(!isnan, bb_results[r].asym)
        ele_asym = filter(!isnan, ele_results[r].asym)

        # Create histogram
        p = histogram(bb_asym, bins=range(0, 1, length=25),
                     alpha=0.6, label="bb0nu (n=$(length(bb_asym)))",
                     xlabel="Asymmetry (E1-E2)/(E1+E2)", ylabel="Counts",
                     title="R = $r mm", color=:blue)
        histogram!(p, ele_asym, bins=range(0, 1, length=25),
                  alpha=0.6, label="electrons (n=$(length(ele_asym)))",
                  color=:red)

        push!(plots_array, p)
    end

    # Combined plot
    p_combined = plot(plots_array..., layout=(1, 3), size=(1400, 400))
    display(p_combined)
    savefig(p_combined, joinpath(pdir, "scripts", "blob_asymmetry_scan.png"))
    println("Saved: blob_asymmetry_scan.png")

    # Summary comparison table
    println("\n" * "=" ^ 70)
    println("SUMMARY COMPARISON TABLE")
    println("=" ^ 70)
    println("\nbb0nu (double beta - 2 Bragg peaks):")
    @printf("%-8s %12s %12s %12s\n", "R (mm)", "Asym mean", "Asym std", "N")
    for r in radii
        valid = filter(!isnan, bb_results[r].asym)
        @printf("%-8.1f %12.3f %12.3f %12d\n", r, mean(valid), std(valid), length(valid))
    end

    println("\nElectrons (single - 1 Bragg peak):")
    @printf("%-8s %12s %12s %12s\n", "R (mm)", "Asym mean", "Asym std", "N")
    for r in radii
        valid = filter(!isnan, ele_results[r].asym)
        @printf("%-8.1f %12.3f %12.3f %12d\n", r, mean(valid), std(valid), length(valid))
    end

    # Separation metric
    println("\n" * "-" ^ 70)
    println("SEPARATION POWER (difference in means / combined std)")
    println("-" ^ 70)
    for r in radii
        bb_valid = filter(!isnan, bb_results[r].asym)
        ele_valid = filter(!isnan, ele_results[r].asym)

        if length(bb_valid) > 0 && length(ele_valid) > 0
            mean_diff = abs(mean(ele_valid) - mean(bb_valid))
            combined_std = sqrt(std(bb_valid)^2 + std(ele_valid)^2)
            separation = mean_diff / combined_std
            @printf("R = %5.1f mm:  separation = %.2f σ\n", r, separation)
        end
    end

    return bb_results, ele_results
end

# Run
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
