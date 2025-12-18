"""
Blob analysis plotting functions for ITACA track analysis.

Provides visualization functions for blob energies, track lengths,
asymmetry, and extreme distances.
"""


"""
    apply_blob_cut(data; eb2cut=400.0)

Apply blob energy cut to signal and background data.

# Arguments
- `data`: NamedTuple with fields `bb` (signal) and `xe` (background) DataFrames
- `eb2cut`: Energy threshold in keV for Eb2 (default: 400.0)

# Returns
NamedTuple with:
- `bb`: Filtered signal DataFrame
- `xe`: Filtered background DataFrame
- `effbb`: Signal efficiency (fraction passing cut)
- `effxe`: Background efficiency (fraction passing cut)
"""
function apply_blob_cut(data; eb2cut=400.0)
    bb, xe = select_blob(data.bb, data.xe; eb2cut=eb2cut)
    eff_bb = size(bb)[1] / size(data.bb)[1]
    eff_xe = size(xe)[1] / size(data.xe)[1]
    (bb=bb, xe=xe, effbb=eff_bb, effxe=eff_xe)
end


"""
    select_blob(bb, xe; eb2cut=400.0)

Filter signal and background DataFrames by Eb2 threshold.

# Arguments
- `bb`: Signal DataFrame with Eb2_keV column
- `xe`: Background DataFrame with Eb2_keV column
- `eb2cut`: Energy threshold in keV (default: 400.0)

# Returns
- `(bbsel, xesel)`: Tuple of filtered DataFrames
"""
function select_blob(bb, xe; eb2cut=400.0)
    bbsel = filter(row -> row.Eb2_keV > eb2cut, bb)
    xesel = filter(row -> row.Eb2_keV > eb2cut, xe)
    return bbsel, xesel
end


"""
    plot_eb2_cut_eff_and_fom(bbdf, xedf; kwargs...)

Plot efficiency and figure of merit vs Eb2 cut threshold.

Creates a 1x2 layout with:
- Left: Signal and background efficiency vs cut
- Right: FOM with optimal cut marked

# Arguments
- `bbdf`: Signal DataFrame (bb0nu)
- `xedf`: Background DataFrame (Xe137)
- `Eb2`: Column name for blob2 energy (default: "Eb2_keV")
- `length`: Number of cut points (default: 100)
- `rinf`: Lower cut range (default: 100)
- `rsup`: Upper cut range (default: 800)
- `fom_bias`: Bias to avoid division by zero in FOM (default: 0.1)
- `suptitle`: Plot super title (default: "Blob2 Energy Cut Analysis")
"""
function plot_eb2_cut_eff_and_fom(bbdf, xedf;
                                  Eb2="Eb2_keV",
                                  length=100,
                                  rinf=100,
                                  rsup=800,
                                  fom_bias=0.1,
                                  suptitle="Blob2 Energy Cut Analysis")
    eb2_sig = bbdf[:, Eb2]
    eb2_bkg = xedf[:, Eb2]
    p1 = plot_efficiency_vs_cut(eb2_sig, eb2_bkg;
                                cuts=range(rinf, rsup, length=length),
                                xlabel="Eblob2 (keV)",
                                title="Efficiency vs Cut",
                                signal_label="bb0nu",
                                background_label="Xe137")
    p2 = plot_fom_vs_cut(eb2_sig, eb2_bkg;
                         cuts=range(rinf, rsup, length=length),
                         xlabel="Eblob2 (keV)",
                         title="FOM vs Cut",
                         mark_optimal=true,
                         bias=fom_bias)
    plot(p1, p2;
         layout=(1, 2),
         size=(1000, 450),
         margin=8Plots.mm,
         left_margin=12Plots.mm,
         bottom_margin=10Plots.mm,
         dpi=150,
         plot_title=suptitle,
         plot_titlefontsize=12)
end


"""
    plot_eb1_vs_eb2(df::DataFrame; kwargs...)

Scatter plot of Eb1 vs Eb2 (blob energies in keV).

# Arguments
- `df`: DataFrame with blob energy columns
- `Eb1`: Column name for blob1 energy (default: "Eb1_keV")
- `Eb2`: Column name for blob2 energy (default: "Eb2_keV")
- `title`: Plot title (default: "Eb1 vs Eb2")
- `xlabel`: X-axis label (default: "Eb1 (keV)")
- `ylabel`: Y-axis label (default: "Eb2 (keV)")
- `kwargs...`: Additional arguments passed to scatter
"""
function plot_eb1_vs_eb2(df::DataFrame;
                         Eb1="Eb1_keV",
                         Eb2="Eb2_keV",
                         title="Eb1 vs Eb2",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)",
                         kwargs...)
    scatter(df[:, Eb1], df[:, Eb2];
            xlabel=xlabel,
            ylabel=ylabel,
            title=title,
            titlefontsize=10,
            legend=false,
            markersize=1,
            markercolor=:black,
            alpha=0.6,
            kwargs...)
end


"""
    plot_track_length(df::DataFrame; kwargs...)

Histogram of track lengths.

# Arguments
- `df`: DataFrame with track_length_mm column
- `title`: Plot title (default: "Track Length Distribution")
- `xlabel`: X-axis label (default: "Track Length (mm)")
- `ylabel`: Y-axis label (default: "Counts")
- `bins`: Number of histogram bins (default: 50)
- `kwargs...`: Additional arguments passed to histogram
"""
function plot_track_length(df::DataFrame;
                           title="Track Length Distribution",
                           xlabel="Track Length (mm)",
                           ylabel="Counts",
                           bins=50,
                           kwargs...)
    histogram(df.track_length_mm;
              xlabel=xlabel,
              ylabel=ylabel,
              title=title,
              legend=false,
              bins=bins,
              kwargs...)
end


"""
    plot_asymmetry(df::DataFrame; kwargs...)

Histogram of blob energy asymmetry.

# Arguments
- `df`: DataFrame with asymmetry column
- `title`: Plot title (default: "Blob Energy Asymmetry")
- `xlabel`: X-axis label (default: "Asymmetry (Eb1-Eb2)/(Eb1+Eb2)")
- `ylabel`: Y-axis label (default: "Counts")
- `bins`: Number of histogram bins (default: 50)
- `kwargs...`: Additional arguments passed to histogram
"""
function plot_asymmetry(df::DataFrame;
                        title="Blob Energy Asymmetry",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50,
                        kwargs...)
    histogram(df.asymmetry;
              xlabel=xlabel,
              ylabel=ylabel,
              title=title,
              legend=false,
              bins=bins,
              kwargs...)
end


"""
    plot_d1_vs_d2(df::DataFrame; kwargs...)

Scatter plot of d1 vs d2 (RECO-MC extreme distances).

# Arguments
- `df`: DataFrame with d1_mm and d2_mm columns
- `title`: Plot title (default: "Extreme Distances: d1 vs d2")
- `xlabel`: X-axis label (default: "d1 (mm)")
- `ylabel`: Y-axis label (default: "d2 (mm)")
- `kwargs...`: Additional arguments passed to scatter
"""
function plot_d1_vs_d2(df::DataFrame;
                       title="Extreme Distances: d1 vs d2",
                       xlabel="d1 (mm)",
                       ylabel="d2 (mm)",
                       kwargs...)
    scatter(df.d1_mm, df.d2_mm;
            xlabel=xlabel,
            ylabel=ylabel,
            title=title,
            legend=false,
            markersize=2,
            alpha=0.6,
            kwargs...)
end
