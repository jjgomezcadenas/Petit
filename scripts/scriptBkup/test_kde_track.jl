#!/usr/bin/env julia

"""
Test energy-weighted KDE along the central path of a track.

Projects voxel energies onto the 1D path coordinate s (arc length),
then uses KDE to estimate the longitudinal energy density f(s).
Maxima of f(s) correspond to Bragg peaks.
"""

const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)

using DataFrames
using Statistics
using Plots
using StatsPlots  # for groupedbar

include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

function main(reco_file::String; event_idx::Int=1, bandwidth::Float64=5.0, rblob::Float64=15.0, binwidth::Float64=3.0)
    println("Loading reconstructed tracks from: $reco_file")

    results, metadata = Petit.read_reco_results_from_hdf5(reco_file)
    println("Loaded $(length(results)) tracks")

    if event_idx > length(results)
        error("Event index $event_idx out of range (max $(length(results)))")
    end

    r = results[event_idx]
    track = r.track
    central_path = r.central_path

    println("\n--- Event $(r.event_id) ---")
    println("Voxels: $(nrow(track.voxels))")
    println("Path points: $(nrow(central_path))")
    println("Track length: $(round(r.track_length, digits=1)) mm")

    # Analyze energy profile using KDE
    profile = Petit.analyze_track_energy_profile(track, central_path; bandwidth=bandwidth)

    println("Arc length range: 0 to $(round(profile.track_length, digits=1)) mm")
    println("Total energy: $(round(sum(profile.energies)*1e3, digits=1)) keV")
    println("KDE bandwidth: $(round(profile.bandwidth, digits=2)) mm")

    println("\n--- Detected maxima (Bragg peak candidates) ---")
    y_range = maximum(profile.f_kde) - minimum(profile.f_kde)
    threshold = 0.1 * y_range
    println("  (min_prominence = 10% of y_range = $(round(threshold, digits=4)))")
    for (rank, i) in enumerate(1:min(3, length(profile.maxima_idx)))
        idx = profile.maxima_idx[i]
        s_max = profile.s_eval[idx]
        f_max = profile.f_kde[idx]
        prom = profile.prominences[i]
        println("  Peak $rank: s = $(round(s_max, digits=1)) mm, f = $(round(f_max, digits=4)), prominence = $(round(prom, digits=4))")
    end

    # Compare integrated energy near each endpoint
    println("\n--- Integrated energy near endpoints (rblob=$rblob mm) ---")
    s_vox = profile.s_voxels
    E_vox = profile.energies

    # Energy near start (s ≈ 0)
    mask_start = s_vox .< rblob
    E_start = sum(E_vox[mask_start]) * 1e3  # keV

    # Energy near end (s ≈ track_length)
    mask_end = s_vox .> (profile.track_length - rblob)
    E_end = sum(E_vox[mask_end]) * 1e3  # keV

    println("  Energy near start (s<$(rblob)mm): $(round(E_start, digits=1)) keV")
    println("  Energy near end (s>$(round(profile.track_length - rblob, digits=1))mm): $(round(E_end, digits=1)) keV")
    println("  Ratio (start/end): $(round(E_start/E_end, digits=2))")

    # Analyze energy concentration at each endpoint
    println("\n--- Energy concentration analysis ---")

    bin_width = binwidth
    n_bins_total = max(1, Int(ceil(profile.track_length / bin_width)))
    println("  Bin width: $(round(bin_width, digits=2)) mm (n_bins=$n_bins_total)")

    # Build histogram over full track (same as p2)
    hist_full = zeros(n_bins_total)
    for (s, e) in zip(s_vox, E_vox)
        bin_idx = min(Int(floor(s / bin_width)) + 1, n_bins_total)
        if bin_idx >= 1
            hist_full[bin_idx] += e
        end
    end

    # How many bins in the endpoint regions?
    n_bins_endpoint = max(1, Int(ceil(rblob / bin_width)))
    println("  Bins in endpoint region (rblob=$(rblob)mm): $n_bins_endpoint")

    # Start region: first n_bins_endpoint bins
    hist_start = hist_full[1:n_bins_endpoint]
    E_start_hist = sum(hist_start) * 1e3
    if maximum(hist_start) > 0
        max_bin_start = maximum(hist_start) * 1e3
        n_significant_start = count(h -> h > 0.1 * maximum(hist_start), hist_start)
        concentration_start = max_bin_start / E_start_hist
    else
        max_bin_start = 0.0
        n_significant_start = 0
        concentration_start = 0.0
    end

    # End region: last n_bins_endpoint bins
    hist_end = hist_full[end-n_bins_endpoint+1:end]
    E_end_hist = sum(hist_end) * 1e3
    if maximum(hist_end) > 0
        max_bin_end = maximum(hist_end) * 1e3
        n_significant_end = count(h -> h > 0.1 * maximum(hist_end), hist_end)
        concentration_end = max_bin_end / E_end_hist
    else
        max_bin_end = 0.0
        n_significant_end = 0
        concentration_end = 0.0
    end

    println("  Start region:")
    println("    Max bin energy: $(round(max_bin_start, digits=1)) keV")
    println("    N significant bins (>10% of max): $n_significant_start")
    println("    Concentration (max_bin/total): $(round(concentration_start*100, digits=1))%")

    println("  End region:")
    println("    Max bin energy: $(round(max_bin_end, digits=1)) keV")
    println("    N significant bins (>10% of max): $n_significant_end")
    println("    Concentration (max_bin/total): $(round(concentration_end*100, digits=1))%")

    # Compute max bin and max/next-to-max ratio for each endpoint
    println("\n--- Max bin analysis ---")

    # Helper function to get max and next-to-max values from histogram
    function get_max_and_next(hist)
        if length(hist) == 0
            return 0.0, 0.0, 0.0
        elseif length(hist) == 1
            return hist[1], 0.0, Inf  # Only one bin, ratio is infinite
        else
            sorted_hist = sort(hist, rev=true)
            max_val = sorted_hist[1]
            next_val = sorted_hist[2]
            ratio = next_val > 0 ? max_val / next_val : Inf
            return max_val, next_val, ratio
        end
    end

    # Start region
    max1_start, max2_start, ratio_start = get_max_and_next(hist_start)
    max1_start_keV = max1_start * 1e3
    max2_start_keV = max2_start * 1e3

    # End region
    max1_end, max2_end, ratio_end = get_max_and_next(hist_end)
    max1_end_keV = max1_end * 1e3
    max2_end_keV = max2_end * 1e3

    println("  Start region:")
    println("    Max bin energy: $(round(max1_start_keV, digits=1)) keV")
    println("    Next-to-max bin energy: $(round(max2_start_keV, digits=1)) keV")
    println("    Ratio (max/next): $(isinf(ratio_start) ? "∞" : round(ratio_start, digits=2))")

    println("  End region:")
    println("    Max bin energy: $(round(max1_end_keV, digits=1)) keV")
    println("    Next-to-max bin energy: $(round(max2_end_keV, digits=1)) keV")
    println("    Ratio (max/next): $(isinf(ratio_end) ? "∞" : round(ratio_end, digits=2))")

    # Interpretation: Higher ratio = more concentrated = more Bragg-like
    println("  Interpretation: Higher ratio → more concentrated → Bragg peak")

    # Compute FWHM of KDE peaks at each endpoint
    println("\n--- FWHM analysis ---")

    # Helper function to compute FWHM around a peak
    function compute_fwhm(s_eval, f_kde, peak_idx)
        peak_val = f_kde[peak_idx]
        half_max = peak_val / 2.0

        # Find left crossing point
        left_idx = peak_idx
        for i in (peak_idx-1):-1:1
            if f_kde[i] < half_max
                left_idx = i
                break
            end
            left_idx = i
        end

        # Find right crossing point
        right_idx = peak_idx
        for i in (peak_idx+1):length(f_kde)
            if f_kde[i] < half_max
                right_idx = i
                break
            end
            right_idx = i
        end

        # Interpolate for more accurate crossing points
        if left_idx > 1 && f_kde[left_idx] < half_max
            # Linear interpolation
            t = (half_max - f_kde[left_idx]) / (f_kde[left_idx+1] - f_kde[left_idx])
            s_left = s_eval[left_idx] + t * (s_eval[left_idx+1] - s_eval[left_idx])
        else
            s_left = s_eval[left_idx]
        end

        if right_idx < length(f_kde) && f_kde[right_idx] < half_max
            t = (half_max - f_kde[right_idx-1]) / (f_kde[right_idx] - f_kde[right_idx-1])
            s_right = s_eval[right_idx-1] + t * (s_eval[right_idx] - s_eval[right_idx-1])
        else
            s_right = s_eval[right_idx]
        end

        return s_right - s_left
    end

    # Find peak in start region (s < rblob)
    start_region_mask = profile.s_eval .< rblob
    if any(start_region_mask)
        start_region_indices = findall(start_region_mask)
        local_max_idx = start_region_indices[argmax(profile.f_kde[start_region_indices])]
        fwhm_start = compute_fwhm(profile.s_eval, profile.f_kde, local_max_idx)
        peak_height_start = profile.f_kde[local_max_idx]
    else
        fwhm_start = 0.0
        peak_height_start = 0.0
    end

    # Find peak in end region (s > track_length - rblob)
    end_region_mask = profile.s_eval .> (profile.track_length - rblob)
    if any(end_region_mask)
        end_region_indices = findall(end_region_mask)
        local_max_idx = end_region_indices[argmax(profile.f_kde[end_region_indices])]
        fwhm_end = compute_fwhm(profile.s_eval, profile.f_kde, local_max_idx)
        peak_height_end = profile.f_kde[local_max_idx]
    else
        fwhm_end = 0.0
        peak_height_end = 0.0
    end

    println("  Start region: FWHM = $(round(fwhm_start, digits=2)) mm, peak height = $(round(peak_height_start, digits=4))")
    println("  End region: FWHM = $(round(fwhm_end, digits=2)) mm, peak height = $(round(peak_height_end, digits=4))")
    println("  FWHM ratio (start/end): $(round(fwhm_start/max(fwhm_end, 0.001), digits=2))")

    # Compute asymmetry metrics
    println("\n--- Asymmetry metrics ---")

    # Energy asymmetry: A_E = (E_max - E_min) / (E_max + E_min)
    # Always positive: 0 = symmetric, 1 = completely asymmetric
    E_max = max(E_end, E_start)
    E_min = min(E_end, E_start)
    E_sum = E_end + E_start
    A_E = E_sum > 0 ? (E_max - E_min) / E_sum : 0.0

    # FWHM asymmetry: A_FWHM = (FWHM_max - FWHM_min) / (FWHM_max + FWHM_min)
    fwhm_max = max(fwhm_end, fwhm_start)
    fwhm_min = min(fwhm_end, fwhm_start)
    fwhm_sum = fwhm_end + fwhm_start
    A_FWHM = fwhm_sum > 0 ? (fwhm_max - fwhm_min) / fwhm_sum : 0.0

    # Peak height asymmetry: A_h = (h_max - h_min) / (h_max + h_min)
    h_max = max(peak_height_end, peak_height_start)
    h_min = min(peak_height_end, peak_height_start)
    h_sum = peak_height_end + peak_height_start
    A_h = h_sum > 0 ? (h_max - h_min) / h_sum : 0.0

    println("  Energy asymmetry A_E = $(round(A_E, digits=3))")
    println("  FWHM asymmetry A_FWHM = $(round(A_FWHM, digits=3))")
    println("  Peak height asymmetry A_h = $(round(A_h, digits=3))")
    println("  (Values: 0 → symmetric/double Bragg, 1 → single Bragg)")

    # Compute Peak-to-Valley ratio
    println("\n--- Peak-to-Valley analysis ---")

    # Find the valley region: s ∈ [rblob, L - rblob]
    valley_mask = (profile.s_eval .>= rblob) .& (profile.s_eval .<= (profile.track_length - rblob))

    if any(valley_mask)
        valley_indices = findall(valley_mask)
        f_valley = profile.f_kde[valley_indices]
        s_valley = profile.s_eval[valley_indices]

        # Valley depth: minimum of f(s) in the middle region
        valley_min_idx = argmin(f_valley)
        V_depth = f_valley[valley_min_idx]
        s_valley_min = s_valley[valley_min_idx]

        # Central flatness: std deviation of f(s) in middle 50% of track
        middle_start = profile.track_length * 0.25
        middle_end = profile.track_length * 0.75
        middle_mask = (profile.s_eval .>= middle_start) .& (profile.s_eval .<= middle_end)
        if any(middle_mask)
            central_flatness = std(profile.f_kde[middle_mask])
        else
            central_flatness = 0.0
        end

        # Peak-to-Valley ratio: PV = (h_start + h_end) / (2 × V)
        PV_ratio = V_depth > 0 ? (peak_height_start + peak_height_end) / (2.0 * V_depth) : Inf

        println("  Valley region: s ∈ [$(round(rblob, digits=1)), $(round(profile.track_length - rblob, digits=1))] mm")
        println("  Valley minimum: V = $(round(V_depth, digits=4)) at s = $(round(s_valley_min, digits=1)) mm")
        println("  Peak heights: h_start = $(round(peak_height_start, digits=4)), h_end = $(round(peak_height_end, digits=4))")
        println("  Peak-to-Valley ratio: PV = $(isinf(PV_ratio) ? "∞" : round(PV_ratio, digits=2))")
        println("  Central flatness (σ): $(round(central_flatness, digits=4))")
        println("  (High PV → deep valley → double Bragg signature)")
    else
        V_depth = 0.0
        s_valley_min = profile.track_length / 2
        PV_ratio = 1.0
        central_flatness = 0.0
        println("  Valley region too small to analyze (track length < 2×rblob)")
    end

    # Compute Skewness and Kurtosis at each endpoint
    println("\n--- Skewness and Kurtosis analysis (Phase 2) ---")

    # Helper function to compute skewness and kurtosis from weighted samples
    function compute_skew_kurt(positions, energies)
        if length(positions) < 3 || sum(energies) == 0
            return 0.0, 0.0
        end

        # Weighted mean
        total_E = sum(energies)
        μ = sum(positions .* energies) / total_E

        # Weighted variance
        σ² = sum(energies .* (positions .- μ).^2) / total_E
        σ = sqrt(σ²)

        if σ < 1e-10
            return 0.0, 0.0
        end

        # Weighted skewness: γ₁ = E[(X-μ)³] / σ³
        m3 = sum(energies .* (positions .- μ).^3) / total_E
        skewness = m3 / σ^3

        # Weighted kurtosis (excess): γ₂ = E[(X-μ)⁴] / σ⁴ - 3
        m4 = sum(energies .* (positions .- μ).^4) / total_E
        kurtosis = m4 / σ^4 - 3.0  # Excess kurtosis (normal = 0)

        return skewness, kurtosis
    end

    # Start region: s < rblob
    s_start_region = s_vox[mask_start]
    E_start_region = E_vox[mask_start]
    skew_start, kurt_start = compute_skew_kurt(s_start_region, E_start_region)

    # End region: s > L - rblob (use distance from end for consistent interpretation)
    s_end_region = profile.track_length .- s_vox[mask_end]  # Distance from end
    E_end_region = E_vox[mask_end]
    skew_end, kurt_end = compute_skew_kurt(s_end_region, E_end_region)

    println("  Start region (s < $(rblob) mm):")
    println("    Skewness γ₁ = $(round(skew_start, digits=3)) (>0: tail toward track center)")
    println("    Kurtosis γ₂ = $(round(kurt_start, digits=3)) (>0: sharper than normal)")

    println("  End region (s > $(round(profile.track_length - rblob, digits=1)) mm):")
    println("    Skewness γ₁ = $(round(skew_end, digits=3)) (>0: tail toward track center)")
    println("    Kurtosis γ₂ = $(round(kurt_end, digits=3)) (>0: sharper than normal)")

    # Interpretation
    println("  Interpretation:")
    println("    High kurtosis → sharp, concentrated peak (Bragg-like)")
    println("    Positive skewness → energy tails toward track center")

    # Compute Entropy (Phase 3)
    println("\n--- Entropy analysis (Phase 3) ---")

    # Helper function to compute entropy from energy distribution
    function compute_entropy(energies)
        if length(energies) == 0 || sum(energies) == 0
            return 0.0
        end
        # Normalize to probabilities
        total_E = sum(energies)
        p = energies ./ total_E
        # Filter out zero probabilities (log(0) is undefined)
        p_nonzero = p[p .> 0]
        # Entropy: H = -Σ pᵢ log(pᵢ)
        H = -sum(p_nonzero .* log.(p_nonzero))
        return H
    end

    # Compute entropy for different regions
    # Full track entropy (using histogram bins)
    H_full = compute_entropy(hist_full)
    H_max_full = log(n_bins_total)  # Maximum possible entropy (uniform distribution)
    H_norm_full = H_max_full > 0 ? H_full / H_max_full : 0.0  # Normalized entropy [0,1]

    # Start region entropy
    H_start = compute_entropy(hist_start)
    H_max_start = log(length(hist_start))
    H_norm_start = H_max_start > 0 ? H_start / H_max_start : 0.0

    # End region entropy
    H_end = compute_entropy(hist_end)
    H_max_end = log(length(hist_end))
    H_norm_end = H_max_end > 0 ? H_end / H_max_end : 0.0

    println("  Full track:")
    println("    Entropy H = $(round(H_full, digits=3)) (max possible: $(round(H_max_full, digits=3)))")
    println("    Normalized H = $(round(H_norm_full, digits=3)) (0=concentrated, 1=uniform)")

    println("  Start region:")
    println("    Entropy H = $(round(H_start, digits=3)), Normalized = $(round(H_norm_start, digits=3))")

    println("  End region:")
    println("    Entropy H = $(round(H_end, digits=3)), Normalized = $(round(H_norm_end, digits=3))")

    # Entropy asymmetry
    H_sum = H_norm_start + H_norm_end
    A_H = H_sum > 0 ? (H_norm_end - H_norm_start) / H_sum : 0.0
    println("  Entropy asymmetry A_H = $(round(A_H, digits=3))")

    println("  Interpretation:")
    println("    Low entropy → energy concentrated (Bragg-like)")
    println("    High entropy → energy spread uniformly")

    # Decision based on concentration (fewer bins = more concentrated = Bragg peak)
    println("\n--- Bragg peak identification ---")
    if concentration_start > concentration_end && n_significant_start <= n_significant_end
        println("  -> Bragg peak at START (more concentrated)")
    elseif concentration_end > concentration_start && n_significant_end <= n_significant_start
        println("  -> Bragg peak at END (more concentrated)")
    elseif n_significant_start < n_significant_end
        println("  -> Bragg peak likely at START (fewer significant bins)")
    elseif n_significant_end < n_significant_start
        println("  -> Bragg peak likely at END (fewer significant bins)")
    else
        println("  -> Ambiguous (similar concentration at both ends)")
    end

    # Plot
    p1 = plot(profile.s_eval, profile.f_kde,
              xlabel="Arc length s (mm)",
              ylabel="Energy density f(s)",
              title="Longitudinal Energy Density (Event $(r.event_id))",
              label="KDE (h=$(round(profile.bandwidth, digits=1))mm)",
              linewidth=2,
              color=:blue)

    # Mark maxima with prominence labels
    for (i, idx) in enumerate(profile.maxima_idx[1:min(3, length(profile.maxima_idx))])
        prom = profile.prominences[i]
        scatter!(p1, [profile.s_eval[idx]], [profile.f_kde[idx]],
                 markersize=8, color=:red, label=(i==1 ? "peaks" : ""))
        annotate!(p1, profile.s_eval[idx], profile.f_kde[idx] + 0.001,
                  text("p=$(round(prom, digits=3))", 8, :center))
    end

    # Also plot histogram of energy vs s
    p2 = histogram(profile.s_voxels, weights=profile.energies .* 1e3,
                   bins=range(0, profile.track_length, length=n_bins_total+1),
                   xlabel="Arc length s (mm)",
                   ylabel="Energy (keV)",
                   title="Energy Histogram",
                   label="",
                   fillalpha=0.7,
                   color=:gray)

    # Bar chart of integrated energy at endpoints
    p3 = bar(["Start (s<$(rblob)mm)", "End (s>$(round(profile.track_length-rblob, digits=0))mm)"],
             [E_start, E_end],
             xlabel="Endpoint",
             ylabel="Integrated Energy (keV)",
             title="Blob Energy (r=$(rblob)mm)",
             label="",
             color=[:blue, :red],
             fillalpha=0.7)

    # Histogram of start region (binned energy)
    if length(hist_start) > 1
        bins_centers_start = range(bin_width/2, rblob - bin_width/2, length=length(hist_start))
    else
        bins_centers_start = [bin_width/2]
    end
    p4 = bar(collect(bins_centers_start), hist_start .* 1e3,
             xlabel="s (mm)",
             ylabel="Energy (keV)",
             title="Start region\nn_bins=$n_significant_start\nconc=$(round(concentration_start*100, digits=0))%",
             label="",
             color=:blue,
             fillalpha=0.7,
             bar_width=bin_width*0.8)

    # Histogram of end region (binned energy)
    if length(hist_end) > 1
        bins_centers_end = range(bin_width/2, rblob - bin_width/2, length=length(hist_end))
    else
        bins_centers_end = [bin_width/2]
    end
    p5 = bar(collect(bins_centers_end), hist_end .* 1e3,
             xlabel="s from end (mm)",
             ylabel="Energy (keV)",
             title="End region\nn_bins=$n_significant_end\nconc=$(round(concentration_end*100, digits=0))%",
             label="",
             color=:red,
             fillalpha=0.7,
             bar_width=bin_width*0.8)

    # FWHM comparison bar chart
    p6 = bar(["Start", "End"],
             [fwhm_start, fwhm_end],
             xlabel="Endpoint",
             ylabel="FWHM (mm)",
             title="KDE Peak Width (FWHM)",
             label="",
             color=[:blue, :red],
             fillalpha=0.7)

    # Asymmetry metrics bar chart
    # Color based on asymmetry level: green (symmetric) to red (asymmetric)
    asymmetry_colors = [A_E < 0.3 ? :green : (A_E < 0.6 ? :orange : :red),
                        A_FWHM < 0.3 ? :green : (A_FWHM < 0.6 ? :orange : :red),
                        A_h < 0.3 ? :green : (A_h < 0.6 ? :orange : :red)]
    p7 = bar(["A_E\n(Energy)", "A_FWHM\n(Width)", "A_h\n(Height)"],
             [A_E, A_FWHM, A_h],
             xlabel="Asymmetry Type",
             ylabel="Asymmetry Value",
             title="Asymmetry Metrics\n(0=symmetric, 1=single Bragg)",
             label="",
             color=asymmetry_colors,
             fillalpha=0.7,
             ylims=(0, 1.1))
    hline!(p7, [0.3], color=:green, linestyle=:dash, label="", linewidth=1)
    hline!(p7, [0.6], color=:orange, linestyle=:dash, label="", linewidth=1)

    # Max bin energy comparison
    p8 = groupedbar(["Start", "End"],
                    [max1_start_keV max2_start_keV; max1_end_keV max2_end_keV],
                    xlabel="Endpoint",
                    ylabel="Energy (keV)",
                    title="Max Bin Energies",
                    label=["Max bin" "Next-to-max"],
                    color=[:darkblue :lightblue],
                    fillalpha=0.7)

    # Max/next-to-max ratio comparison
    # Handle infinite ratios for plotting
    ratio_start_plot = isinf(ratio_start) ? 10.0 : min(ratio_start, 10.0)
    ratio_end_plot = isinf(ratio_end) ? 10.0 : min(ratio_end, 10.0)
    p9 = bar(["Start", "End"],
             [ratio_start_plot, ratio_end_plot],
             xlabel="Endpoint",
             ylabel="Ratio (max/next)",
             title="Max/Next-to-Max Ratio\n(higher → more concentrated)",
             label="",
             color=[:blue, :red],
             fillalpha=0.7)
    # Add reference line at ratio=2 (significant concentration)
    hline!(p9, [2.0], color=:green, linestyle=:dash, label="", linewidth=1)

    # Peak-to-Valley visualization
    p10 = plot(profile.s_eval, profile.f_kde,
               xlabel="Arc length s (mm)",
               ylabel="Energy density f(s)",
               title="Peak-to-Valley Analysis",
               label="KDE",
               linewidth=2,
               color=:blue)
    # Mark the valley region
    vspan!(p10, [rblob, profile.track_length - rblob],
           color=:gray, alpha=0.2, label="Valley region")
    # Mark the valley minimum
    scatter!(p10, [s_valley_min], [V_depth],
             markersize=10, color=:green, markershape=:diamond, label="Valley min")
    # Mark endpoint peaks
    scatter!(p10, [0, profile.track_length], [peak_height_start, peak_height_end],
             markersize=8, color=[:blue, :red], label="")
    # Draw horizontal line at valley depth
    hline!(p10, [V_depth], color=:green, linestyle=:dash, alpha=0.5, label="")

    # PV ratio bar (single value display)
    PV_plot = isinf(PV_ratio) ? 10.0 : min(PV_ratio, 10.0)
    p11 = bar(["PV Ratio"],
              [PV_plot],
              xlabel="",
              ylabel="Peak-to-Valley Ratio",
              title="PV = $(isinf(PV_ratio) ? "∞" : round(PV_ratio, digits=2))\n(high → double Bragg)",
              label="",
              color=PV_ratio > 2.0 ? :green : :orange,
              fillalpha=0.7,
              ylims=(0, max(PV_plot * 1.2, 3.0)))
    # Reference lines
    hline!(p11, [1.5], color=:orange, linestyle=:dash, label="", linewidth=1)
    hline!(p11, [2.5], color=:green, linestyle=:dash, label="", linewidth=1)

    # Skewness comparison
    skew_colors = [skew_start >= 0 ? :blue : :lightblue,
                   skew_end >= 0 ? :red : :lightcoral]
    p12 = bar(["Start", "End"],
              [skew_start, skew_end],
              xlabel="Endpoint",
              ylabel="Skewness γ₁",
              title="Skewness\n(>0: tail toward center)",
              label="",
              color=skew_colors,
              fillalpha=0.7)
    hline!(p12, [0], color=:black, linestyle=:dash, label="", linewidth=1)

    # Kurtosis comparison
    kurt_colors = [kurt_start >= 0 ? :darkgreen : :lightgreen,
                   kurt_end >= 0 ? :darkgreen : :lightgreen]
    p13 = bar(["Start", "End"],
              [kurt_start, kurt_end],
              xlabel="Endpoint",
              ylabel="Kurtosis γ₂",
              title="Excess Kurtosis\n(>0: sharper than normal)",
              label="",
              color=kurt_colors,
              fillalpha=0.7)
    hline!(p13, [0], color=:black, linestyle=:dash, label="", linewidth=1)
    hline!(p13, [1], color=:green, linestyle=:dot, label="", linewidth=1, alpha=0.5)

    # Entropy comparison (normalized)
    p14 = bar(["Full Track", "Start", "End"],
              [H_norm_full, H_norm_start, H_norm_end],
              xlabel="Region",
              ylabel="Normalized Entropy",
              title="Entropy (normalized)\n(0=concentrated, 1=uniform)",
              label="",
              color=[:gray, :blue, :red],
              fillalpha=0.7,
              ylims=(0, 1.1))
    hline!(p14, [0.5], color=:orange, linestyle=:dash, label="", linewidth=1)

    # Entropy asymmetry (single bar like PV)
    A_H_color = abs(A_H) < 0.3 ? :green : :orange
    p15 = bar(["A_H"],
              [A_H],
              xlabel="",
              ylabel="Entropy Asymmetry",
              title="A_H = $(round(A_H, digits=3))\n(0=symmetric)",
              label="",
              color=A_H_color,
              fillalpha=0.7,
              ylims=(-1.1, 1.1))
    hline!(p15, [0], color=:black, linestyle=:dash, label="", linewidth=1)
    hline!(p15, [-0.3, 0.3], color=:green, linestyle=:dot, label="", linewidth=1, alpha=0.5)

    # Combined plot - 9 rows, each feature in its own row
    # Row 1: KDE (full width)
    # Row 2: Energy histogram (full width)
    # Row 3: Integrated energy (left) + FWHM (right)
    # Row 4: Start concentration (left) + End concentration (right)
    # Row 5: Asymmetry metrics (full width)
    # Row 6: Max bin energies (left) + Max/next ratio (right)
    # Row 7: Peak-to-Valley visualization (left) + PV ratio (right)
    # Row 8: Skewness (left) + Kurtosis (right)
    # Row 9: Entropy (left) + Entropy asymmetry (right)
    p = plot(p1, p2, p3, p6, p4, p5, p7, p8, p9, p10, p11, p12, p13, p14, p15,
             layout=@layout([a; b; c d; e f; g; h i; j k; l m; n o]),
             size=(1000, 2100))

    output_file = joinpath(dirname(reco_file), "kde_track_test.png")
    savefig(p, output_file)
    println("\nPlot saved to: $output_file")

    return profile
end

# Parse command line
if length(ARGS) < 1
    println("Usage: julia test_kde_track.jl <reco_file.h5> [--event=N] [--bandwidth=H] [--rblob=R] [--binwidth=W]")
    println("  --event=N      Event number in file, 1-based (default: 1 = first event)")
    println("  --bandwidth=H  KDE bandwidth in mm (default: 5.0)")
    println("  --rblob=R      Blob integration radius in mm (default: 15.0)")
    println("  --binwidth=W   Histogram bin width in mm (default: 3.0)")
    exit(0)
end

reco_file = ARGS[1]
event_idx = 1
bandwidth = 5.0
rblob = 15.0
binwidth = 3.0

for arg in ARGS[2:end]
    if startswith(arg, "--event=")
        global event_idx = parse(Int, split(arg, "=")[2])
    elseif startswith(arg, "--bandwidth=")
        global bandwidth = parse(Float64, split(arg, "=")[2])
    elseif startswith(arg, "--rblob=")
        global rblob = parse(Float64, split(arg, "=")[2])
    elseif startswith(arg, "--binwidth=")
        global binwidth = parse(Float64, split(arg, "=")[2])
    end
end

main(reco_file; event_idx=event_idx, bandwidth=bandwidth, rblob=rblob, binwidth=binwidth)
