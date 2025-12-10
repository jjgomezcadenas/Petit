"""
Figure of Merit and ROC Analysis Module

Provides functions for evaluating signal/background discrimination:
- Efficiency curves
- Figure of merit optimization
- ROC curves with AUC
"""


"""
    efficiency_vs_cut(data; cuts) -> NamedTuple

Compute efficiency as a function of cut threshold.

# Arguments
- `data::AbstractVector{<:Real}`: Vector of values to cut on
- `cuts::AbstractVector{<:Real}`: Vector of cut thresholds (default: 0 to 1000 in 101 steps)

# Returns
NamedTuple with fields:
- `cuts`: Vector of cut values
- `efficiency`: Vector of efficiencies (fraction passing data < cut)
- `n_total`: Total number of events
- `n_passing`: Vector of events passing each cut
"""
function efficiency_vs_cut(data::AbstractVector{<:Real};
                           cuts::AbstractVector{<:Real} = range(0, 1000, length=101))
    n_total = length(data)
    n_passing = [count(x -> x > cut, data) for cut in cuts]
    efficiency = n_passing ./ n_total

    return (cuts = collect(cuts),
            efficiency = efficiency,
            n_total = n_total,
            n_passing = n_passing)
end


"""
    compute_efficiencies(signal, background; cuts) -> NamedTuple

Compute signal and background efficiencies for a range of cuts.

# Arguments
- `signal::AbstractVector{<:Real}`: Vector of signal values
- `background::AbstractVector{<:Real}`: Vector of background values
- `cuts::AbstractVector{<:Real}`: Vector of cut thresholds

# Returns
NamedTuple with fields:
- `cuts`: Vector of cut values
- `signal_eff`: Vector of signal efficiencies
- `background_eff`: Vector of background efficiencies
- `n_signal`: Total signal events
- `n_background`: Total background events
"""
function compute_efficiencies(signal::AbstractVector{<:Real},
                              background::AbstractVector{<:Real};
                              cuts::AbstractVector{<:Real} = range(0, 1000, length=101))
    sig_result = efficiency_vs_cut(signal; cuts=cuts)
    bkg_result = efficiency_vs_cut(background; cuts=cuts)

    return (cuts = sig_result.cuts,
            signal_eff = sig_result.efficiency,
            background_eff = bkg_result.efficiency,
            n_signal = sig_result.n_total,
            n_background = bkg_result.n_total)
end


"""
    figure_of_merit(signal, background; cuts, epsilon) -> NamedTuple

Compute figure of merit (FOM = signal_eff / sqrt(background_eff)) vs cut.

# Arguments
- `signal::AbstractVector{<:Real}`: Vector of signal values
- `background::AbstractVector{<:Real}`: Vector of background values
- `cuts::AbstractVector{<:Real}`: Vector of cut thresholds
- `epsilon::Float64`: Small value to avoid division by zero (default: 1e-10)

# Returns
NamedTuple with fields:
- `cuts`: Vector of cut values
- `fom`: Vector of FOM values
- `signal_eff`: Vector of signal efficiencies
- `background_eff`: Vector of background efficiencies
- `optimal_cut`: Cut value that maximizes FOM
- `optimal_fom`: Maximum FOM value
- `optimal_signal_eff`: Signal efficiency at optimal cut
- `optimal_background_eff`: Background efficiency at optimal cut
"""
function figure_of_merit(signal::AbstractVector{<:Real},
                         background::AbstractVector{<:Real};
                         cuts::AbstractVector{<:Real} = range(0, 1000, length=101),
                         epsilon::Float64 = 1e-10)
    eff = compute_efficiencies(signal, background; cuts=cuts)

    fom = eff.signal_eff ./ sqrt.(eff.background_eff .+ epsilon)

    # Find optimal cut
    idx_optimal = argmax(fom)

    return (cuts = eff.cuts,
            fom = fom,
            signal_eff = eff.signal_eff,
            background_eff = eff.background_eff,
            optimal_cut = eff.cuts[idx_optimal],
            optimal_fom = fom[idx_optimal],
            optimal_signal_eff = eff.signal_eff[idx_optimal],
            optimal_background_eff = eff.background_eff[idx_optimal])
end


"""
    roc_curve(signal, background; cuts) -> NamedTuple

Compute ROC curve: background rejection vs signal efficiency.

# Arguments
- `signal::AbstractVector{<:Real}`: Vector of signal values
- `background::AbstractVector{<:Real}`: Vector of background values
- `cuts::AbstractVector{<:Real}`: Vector of cut thresholds

# Returns
NamedTuple with fields:
- `signal_eff`: Vector of signal efficiencies (x-axis for ROC)
- `background_rejection`: Vector of background rejection (1 - background_eff)
- `background_eff`: Vector of background efficiencies
- `cuts`: Vector of cut values
- `auc`: Area Under the ROC Curve
"""
function roc_curve(signal::AbstractVector{<:Real},
                   background::AbstractVector{<:Real};
                   cuts::AbstractVector{<:Real} = range(0, 1000, length=101))
    eff = compute_efficiencies(signal, background; cuts=cuts)

    background_rejection = 1.0 .- eff.background_eff

    # Compute AUC using trapezoidal integration
    # Sort by signal_eff for proper integration
    sorted_idx = sortperm(eff.signal_eff)
    x_sorted = eff.signal_eff[sorted_idx]
    y_sorted = background_rejection[sorted_idx]

    # Trapezoidal integration
    auc = 0.0
    for i in 2:length(x_sorted)
        dx = x_sorted[i] - x_sorted[i-1]
        y_avg = (y_sorted[i] + y_sorted[i-1]) / 2
        auc += dx * y_avg
    end

    return (signal_eff = eff.signal_eff,
            background_rejection = background_rejection,
            background_eff = eff.background_eff,
            cuts = eff.cuts,
            auc = auc)
end


"""
    plot_efficiency_vs_cut(signal, background; kwargs...) -> Plot

Plot signal and background efficiency vs cut threshold.

# Arguments
- `signal::AbstractVector{<:Real}`: Vector of signal values
- `background::AbstractVector{<:Real}`: Vector of background values
- `cuts::AbstractVector{<:Real}`: Vector of cut thresholds
- `xlabel::String`: X-axis label (default: "Cut Value")
- `title::String`: Plot title (default: "Efficiency vs Cut")
- `signal_label::String`: Legend label for signal (default: "Signal")
- `background_label::String`: Legend label for background (default: "Background")
"""
function plot_efficiency_vs_cut(signal::AbstractVector{<:Real},
                                background::AbstractVector{<:Real};
                                cuts::AbstractVector{<:Real} = range(0, 1000, length=101),
                                xlabel::String = "Cut Value",
                                title::String = "Efficiency vs Cut",
                                signal_label::String = "Signal",
                                background_label::String = "Background")
    eff = compute_efficiencies(signal, background; cuts=cuts)

    p = plot(eff.cuts, eff.signal_eff,
             label=signal_label,
             color=:blue,
             linewidth=2,
             xlabel=xlabel,
             ylabel="Efficiency",
             title=title,
             legend=:right,
             grid=true)

    plot!(p, eff.cuts, eff.background_eff,
          label=background_label,
          color=:red,
          linewidth=2)

    return p
end


"""
    plot_fom_vs_cut(signal, background; kwargs...) -> Plot

Plot figure of merit vs cut threshold with optimal point marked.

# Arguments
- `signal::AbstractVector{<:Real}`: Vector of signal values
- `background::AbstractVector{<:Real}`: Vector of background values
- `cuts::AbstractVector{<:Real}`: Vector of cut thresholds
- `xlabel::String`: X-axis label (default: "Cut Value")
- `title::String`: Plot title (default: "Figure of Merit vs Cut")
- `mark_optimal::Bool`: Whether to mark optimal point (default: true)
"""
function plot_fom_vs_cut(signal::AbstractVector{<:Real},
                         background::AbstractVector{<:Real};
                         cuts::AbstractVector{<:Real} = range(0, 1000, length=101),
                         xlabel::String = "Cut Value",
                         title::String = "Figure of Merit vs Cut",
                         mark_optimal::Bool = true)
    fom_result = figure_of_merit(signal, background; cuts=cuts)

    p = plot(fom_result.cuts, fom_result.fom,
             label="FOM",
             color=:green,
             linewidth=2,
             xlabel=xlabel,
             ylabel="FOM = sig_eff / sqrt(bkg_eff)",
             title=title,
             legend=:topright,
             grid=true)

    if mark_optimal
        # Mark optimal point
        scatter!(p, [fom_result.optimal_cut], [fom_result.optimal_fom],
                 label="Optimal",
                 color=:red,
                 markersize=8,
                 markershape=:star5)

        # Add vertical line at optimal cut
        vline!(p, [fom_result.optimal_cut],
               label="",
               color=:red,
               linestyle=:dash,
               alpha=0.5)

        # Annotation
        annotate!(p, fom_result.optimal_cut, fom_result.optimal_fom * 0.9,
                  text("cut=$(round(fom_result.optimal_cut, digits=1))\nFOM=$(round(fom_result.optimal_fom, digits=3))",
                       8, :left))
    end

    return p
end


"""
    plot_roc(signal, background; kwargs...) -> Plot

Plot ROC curve (background rejection vs signal efficiency).

# Arguments
- `signal::AbstractVector{<:Real}`: Vector of signal values
- `background::AbstractVector{<:Real}`: Vector of background values
- `cuts::AbstractVector{<:Real}`: Vector of cut thresholds
- `title::String`: Plot title (default: "ROC Curve")
- `show_auc::Bool`: Whether to show AUC annotation (default: true)
- `show_diagonal::Bool`: Whether to show diagonal reference line (default: true)
"""
function plot_roc(signal::AbstractVector{<:Real},
                  background::AbstractVector{<:Real};
                  cuts::AbstractVector{<:Real} = range(0, 1000, length=101),
                  title::String = "ROC Curve",
                  show_auc::Bool = true,
                  show_diagonal::Bool = true)
    roc = roc_curve(signal, background; cuts=cuts)

    p = plot(roc.signal_eff, roc.background_rejection,
             label="ROC",
             color=:purple,
             linewidth=2,
             xlabel="Signal Efficiency",
             ylabel="Background Rejection",
             title=title,
             legend=:bottomleft,
             grid=true,
             xlims=(0, 1),
             ylims=(0, 1),
             aspect_ratio=:equal)

    if show_diagonal
        plot!(p, [0, 1], [0, 1],
              label="Random",
              color=:gray,
              linestyle=:dash,
              alpha=0.5)
    end

    if show_auc
        annotate!(p, 0.6, 0.2,
                  text("AUC = $(round(roc.auc, digits=4))", 10, :left))
    end

    return p
end


"""
    plot_all_fom_analysis(signal, background; kwargs...) -> Plot

Create 2x2 layout with efficiency, FOM, ROC, and distribution plots.

# Arguments
- `signal::AbstractVector{<:Real}`: Vector of signal values
- `background::AbstractVector{<:Real}`: Vector of background values
- `cuts::AbstractVector{<:Real}`: Vector of cut thresholds
- `size::Tuple{Int,Int}`: Figure size (default: (1000, 800))
- `xlabel::String`: X-axis label for cut plots (default: "Cut Value")
- `signal_label::String`: Legend label for signal (default: "Signal")
- `background_label::String`: Legend label for background (default: "Background")
"""
function plot_all_fom_analysis(signal::AbstractVector{<:Real},
                               background::AbstractVector{<:Real};
                               cuts::AbstractVector{<:Real} = range(0, 1000, length=101),
                               size::Tuple{Int,Int} = (1000, 800),
                               xlabel::String = "Cut Value",
                               signal_label::String = "Signal",
                               background_label::String = "Background")

    # Efficiency plot
    p1 = plot_efficiency_vs_cut(signal, background;
                                cuts=cuts,
                                xlabel=xlabel,
                                signal_label=signal_label,
                                background_label=background_label)

    # FOM plot
    p2 = plot_fom_vs_cut(signal, background;
                         cuts=cuts,
                         xlabel=xlabel)

    # ROC curve
    p3 = plot_roc(signal, background; cuts=cuts)

    # Distribution histogram
    p4 = histogram(signal,
                   label=signal_label,
                   alpha=0.6,
                   color=:blue,
                   normalize=:probability,
                   xlabel=xlabel,
                   ylabel="Probability",
                   title="Distributions")
    histogram!(p4, background,
               label=background_label,
               alpha=0.6,
               color=:red,
               normalize=:probability)

    # Combine into 2x2 layout
    p = plot(p1, p2, p3, p4,
             layout=(2, 2),
             size=size,
             margin=5Plots.mm)

    return p
end


"""
    print_fom_summary(signal, background; kwargs...)

Print summary statistics for FOM analysis.

# Arguments
- `signal::AbstractVector{<:Real}`: Vector of signal values
- `background::AbstractVector{<:Real}`: Vector of background values
- `cuts::AbstractVector{<:Real}`: Vector of cut thresholds
- `signal_label::String`: Label for signal (default: "Signal")
- `background_label::String`: Label for background (default: "Background")
"""
function print_fom_summary(signal::AbstractVector{<:Real},
                           background::AbstractVector{<:Real};
                           cuts::AbstractVector{<:Real} = range(0, 1000, length=101),
                           signal_label::String = "Signal",
                           background_label::String = "Background")

    fom_result = figure_of_merit(signal, background; cuts=cuts)
    roc = roc_curve(signal, background; cuts=cuts)

    bkg_rejection = 1.0 - fom_result.optimal_background_eff

    println("="^63)
    println("                    FOM Analysis Summary")
    println("="^63)
    println("$signal_label events:".*(repeat(" ", max(0, 20-length(signal_label)-8)))*"$(length(signal))")
    println("$background_label events:".*(repeat(" ", max(0, 20-length(background_label)-8)))*"$(length(background))")
    println("-"^63)
    println("Optimal cut:        $(round(fom_result.optimal_cut, digits=2))")
    println("-"^63)
    println("At optimal cut:")
    println("  Signal efficiency:      $(round(fom_result.optimal_signal_eff, digits=4)) ($(round(100*fom_result.optimal_signal_eff, digits=2))%)")
    println("  Background efficiency:  $(round(fom_result.optimal_background_eff, digits=4)) ($(round(100*fom_result.optimal_background_eff, digits=2))%)")
    println("  Background rejection:   $(round(bkg_rejection, digits=4)) ($(round(100*bkg_rejection, digits=2))%)")
    println("  Figure of Merit:        $(round(fom_result.optimal_fom, digits=4))")
    println("-"^63)
    println("ROC AUC:            $(round(roc.auc, digits=4))")
    println("="^63)
end
