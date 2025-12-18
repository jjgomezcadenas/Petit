# Plot functions for ITACA analysis results
# Usage: include("plot_analysis.jl") in Pluto notebook

using DataFrames
using CSV
using Plots
using Statistics
using Glob
using Petit

"""
    load_analysis_results(filepath::String) -> DataFrame

Load analysis_results.csv into a DataFrame.
"""
function load_analysis_results(filepath::String)
    CSV.read(filepath, DataFrame)
end


function get_data(path; data_type="ion")
	has(patterns...) = f -> all(p -> occursin(p, f), patterns)
	
	pattern = "*$(data_type)*.csv"
	files = glob(pattern, path)

	(bb = load_analysis_results(filter(has("analysis", "bb"), files)[1]),
  	 xe = load_analysis_results(filter(has("analysis", "xe"), files)[1]),
  	 bbm = load_analysis_results(filter(has("metadata", "bb"), files)[1]),
  	 xem = load_analysis_results(filter(has("metadata", "xe"), files)[1]))
	
end

function select_kde(bbdf, xedf; pmin=0.5, lmax = 50, rmin=150)
	
	function is_signal_like(ev; p_min, L_max, R_min)
	    p1 = ev.peak1_prom
	    p2 = ev.peak2_prom
	    L1 = ev.peak1_left
	    R2 = ev.peak2_right

    	left_ok  = (L1 ≤ lmax)
    	right_ok = (p2 ≥ p_min) && (R2 ≥ rmin)

    	return left_ok && right_ok
	end

	# --- Apply to BB0ν (signal) ---
	bbsel = filter(ev -> is_signal_like(ev; p_min=pmin,
										 L_max=lmax,
										 R_min=rmin),
					                     bbdf)
	xesel = filter(ev -> is_signal_like(ev; p_min=pmin,
										 L_max=lmax,
										 R_min=rmin),
				   						 xedf)
	return bbsel, xesel
end


function kde_eff(bbdf, xedf; lmax = 50, rmin=150)
	EFBB = []
	EFXE = []
	range = 0.2:0.1:1.0
	for p in 0.2:0.1:1.0
		bbsel, xesel =select_kde(bbdf, xedf; pmin=p, lmax = lmax, rmin=rmin)
		push!(EFBB, size(bbsel)[1]/size(bbdf)[1])
		push!(EFXE, size(xesel)[1]/size(xedf)[1])
	end
	collect(range), EFBB, EFXE
	
end


function apply_ked_cut(data; pmin = 0.4, lmax = 50, rmin = 150)
	
	bb, xe =select_kde(data.bb, data.xe; pmin=pmin, lmax = lmax,  
					   rmin=rmin)
	eff_kde_bb = size(bb)[1]/size(data.bb)[1]
	eff_kde_xe = size(xe)[1]/size(data.xe)[1]
	(bb=bb, xe=xe, effbb = eff_kde_bb, effxe = eff_kde_xe)
	
end



function apply_blob_cut(data; eb2cut=400.0)
	bb, xe =select_blob(data.bb, data.xe, eb2cut=eb2cut)
										
    eff_bb = size(bb)[1]/size(data.bb)[1]
	eff_xe = size(xe)[1]/size(data.xe)[1]
	(bb=bb, xe=xe, effbb = eff_bb, effxe = eff_xe)
end
	

function apply_kde_and_blob_cut(data; pmin = 0.4, lmax = 50, rmin = 150, eb2cut=400.0)
	
    kde =apply_ked_cut(data; pmin = pmin, lmax = lmax, rmin = rmin)
    blob =apply_blob_cut(kde; eb2cut=eb2cut)
										
	(bb=blob.bb, xe=blob.xe, effbb = blob.effbb * kde.effbb, effxe = blob.effxe * kde.effxe)
end



function select_blob(bb, xe; pmin = 0.4, lmax = 50, rmin = 150, eb2cut=400.0)
		bbsel = filter(row -> row.Eb2_keV > eb2cut, bb )
		xesel = filter(row -> row.Eb2_keV > eb2cut, xe )
		return bbsel, xesel
end


####
## plots
###

function plot_peak_pos_vs_prom(bb, xe; ylims=(0, 1.5))
	sbb = scatter(bb.peak1_left, bb.peak1_prom;
            xlabel="peak1 left",
            ylabel="peak1 prom",
            title="bb0nu peak1 vs prom1",
            legend=false,
            markersize=2,
            alpha=0.6,  ylims=ylims)
	sxe = scatter(xe.peak1_left, xe.peak1_prom;
            xlabel="peak1 left",
            ylabel="peak1 prom",
            title="xe137 peak1 vs prom1",
            legend=false,
            markersize=2,
            alpha=0.6, ylims=ylims)
	sbb2 = scatter(bb.peak2_right, bb.peak2_prom;
            xlabel="peak2 right",
            ylabel="peak2 prom",
            title="bb0nu peak2 vs prom2",
            legend=false,
            markersize=2,
            alpha=0.6, ylims=ylims)
	sxe2 = scatter(xe.peak2_right, xe.peak2_prom;
            xlabel="peak2 right",
            ylabel="peak2 prom",
            title="xe137 peak2 vs prom2",
            legend=false,
            markersize=2,
            alpha=0.6, ylims=ylims)
	plot(sbb, sxe, sbb2, sxe2, layout=(2,2), size=(1200, 1000), margin=6Plots.mm)
end


function plot_kdef_eff_and_fom(bb, xe; lmax = 50, rmin=150)
	cuts, effbb, effxe = kde_eff(bb, xe; lmax = lmax, rmin=rmin)
	p = plot(cuts, effbb,
             label="bb0nu",
             color=:blue,
             linewidth=2,
             xlabel="prom",
             ylabel="Efficiency",
             title="Efficiency KDE",
             legend=:best,
             grid=true)

    p2 = plot!(p, cuts, effxe,
          label="Xe137",
          color=:red,
          linewidth=2)
	p3 = plot(cuts, effbb./sqrt.(effxe),
             label="fomd",
             color=:blue,
             linewidth=2,
             xlabel="prom",
             ylabel="Efficiency",
             title="fom KDE",
             legend=:best,
             grid=true)
	plot(p2, p3)
end


function plot_eb2_cut_eff_and_fom(bb,xe; length=100)
	eb2_sig = bb.Eb2_keV
	eb2_bkg = xe.Eb2_keV
	p1 = Petit.plot_efficiency_vs_cut(eb2_sig,
                                eb2_bkg;
                                cuts = range(200, 600, length=length),
                                xlabel = "Eblob2",
                                title = "Efficiency vs Cut",
                                signal_label = "bb0nu",
                                background_label = "Xe137")
	p2 = Petit.plot_fom_vs_cut(eb2_sig,
                         eb2_bkg;
                         cuts = range(200, 600, length=length),
                         xlabel = "Eblob2",
                         title= "FOM vs Cut",
                         mark_optimal = true)
	plot(p1, p2)

end


function plot_peak1_prom_vs_peak2_prom(bb,xe)
	
	sbb1 = scatter(bb.peak1_prom,bb.peak2_prom;
	            xlabel="peak1 prom",
	            ylabel="peak2 prom",
	            title="bb0nu p1 vs p2 proms",
	            legend=false,
	            markersize=2,
	            alpha=0.6, ylims=(0, 1.5), xlims=(0, 1.4))
	sxe1 = scatter(xe.peak1_prom,xe.peak2_prom;
	            xlabel="peak1 prom",
	            ylabel="peak2 prom",
	            title="xe137 p1 vs p2 proms",
	            legend=false,
	            markersize=2,
	            alpha=0.6, ylims=(0, 1.5), xlims=(0, 1.4))
	plot(sbb1, sxe1)
end



"""
    plot_eb1_vs_eb2(df::DataFrame; kwargs...)

Scatter plot of Eb1 vs Eb2 (blob energies in keV).
"""
function plot_eb1_vs_eb2(df::DataFrame;
                          title="Blob Energies: Eb1 vs Eb2",
                          xlabel="Eb1 (keV)",
                          ylabel="Eb2 (keV)",
                          kwargs...)
    scatter(df.Eb1_keV, df.Eb2_keV;
            xlabel=xlabel,
            ylabel=ylabel,
            title=title,
            legend=false,
            markersize=2,
            alpha=0.6,
            kwargs...)
end

"""
    plot_track_length(df::DataFrame; kwargs...)

Histogram of track lengths.
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

"""
    plot_n_peaks(df::DataFrame; kwargs...)

Histogram of number of KDE peaks found.
"""
function plot_n_peaks(df::DataFrame;
                      title="Number of KDE Peaks",
                      xlabel="N peaks",
                      ylabel="Counts",
                      kwargs...)
    histogram(df.n_peaks;
              xlabel=xlabel,
              ylabel=ylabel,
              title=title,
              legend=false,
              bins=range(0.5, maximum(df.n_peaks) + 0.5, step=1),
              kwargs...)
end

"""
    plot_peak_position_vs_prom(df::DataFrame; kwargs...)

Scatter plot of peak position vs prominence for both peak1 and peak2.
Peak1 shown as circles, Peak2 shown as squares.
Position is computed as midpoint of left and right edges.
"""
function plot_peak_position_vs_prom(df::DataFrame;
                                    title="Peak Position vs Prominence",
                                    xlabel="Position (mm)",
                                    ylabel="Prominence (normalized)",
                                    kwargs...)
    # Compute peak positions as midpoint of left and right edges
    pos1 = (df.peak1_left .+ df.peak1_right) ./ 2
    pos2 = (df.peak2_left .+ df.peak2_right) ./ 2

    # Filter out zero values (peaks not found)
    mask1 = df.peak1_prom .> 0
    mask2 = df.peak2_prom .> 0

    p = scatter(pos1[mask1], df.peak1_prom[mask1];
                xlabel=xlabel,
                ylabel=ylabel,
                title=title,
                label="Peak 1 (left)",
                marker=:circle,
                markersize=5,
                alpha=0.6,
                kwargs...)

    scatter!(p, pos2[mask2], df.peak2_prom[mask2];
             label="Peak 2 (right)",
             marker=:square,
             markersize=2,
             alpha=0.6)

    return p
end

"""
    plot_all_analysis(df::DataFrame; size=(1200, 800))

Create a 2x3 layout with all analysis plots.
"""
function plot_all_analysis(df::DataFrame; size=(1200, 800))
    p1 = plot_eb1_vs_eb2(df)
    p2 = plot_track_length(df)
    p3 = plot_asymmetry(df)
    p4 = plot_d1_vs_d2(df)
    p5 = plot_n_peaks(df)
    p6 = plot_peak_position_vs_prom(df)

    plot(p1, p2, p3, p4, p5, p6;
         layout=(2, 3),
         size=size,
         margin=5Plots.mm)
end

# Summary statistics
"""
    print_analysis_summary(df::DataFrame)

Print summary statistics for the analysis results.
"""
function print_analysis_summary(df::DataFrame)
    println("═══════════════════════════════════════════════════════════════")
    println("              Analysis Results Summary")
    println("═══════════════════════════════════════════════════════════════")
    println("Number of events: $(nrow(df))")
    println("───────────────────────────────────────────────────────────────")
    println("Track Length (mm):")
    println("  mean = $(round(mean(df.track_length_mm), digits=2))")
    println("  std  = $(round(std(df.track_length_mm), digits=2))")
    println("  min  = $(round(minimum(df.track_length_mm), digits=2))")
    println("  max  = $(round(maximum(df.track_length_mm), digits=2))")
    println("───────────────────────────────────────────────────────────────")
    println("Blob Energies (keV):")
    println("  Eb1: mean = $(round(mean(df.Eb1_keV), digits=1)), std = $(round(std(df.Eb1_keV), digits=1))")
    println("  Eb2: mean = $(round(mean(df.Eb2_keV), digits=1)), std = $(round(std(df.Eb2_keV), digits=1))")
    println("───────────────────────────────────────────────────────────────")
    println("Asymmetry:")
    println("  mean = $(round(mean(df.asymmetry), digits=4))")
    println("  std  = $(round(std(df.asymmetry), digits=4))")
    println("───────────────────────────────────────────────────────────────")
    println("Extreme Distances (mm):")
    println("  d1: mean = $(round(mean(df.d1_mm), digits=2)), std = $(round(std(df.d1_mm), digits=2))")
    println("  d2: mean = $(round(mean(df.d2_mm), digits=2)), std = $(round(std(df.d2_mm), digits=2))")
    println("───────────────────────────────────────────────────────────────")
    println("N peaks distribution:")
    for n in sort(unique(df.n_peaks))
        count = sum(df.n_peaks .== n)
        pct = round(100 * count / nrow(df), digits=1)
        println("  $n peaks: $count events ($pct%)")
    end
    println("═══════════════════════════════════════════════════════════════")
end
