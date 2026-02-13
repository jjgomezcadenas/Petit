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


function apply_blob_cut(data; eb2cut=400.0)
	bb, xe =select_blob(data.bb, data.xe, eb2cut=eb2cut)
										
    eff_bb = size(bb)[1]/size(data.bb)[1]
	eff_xe = size(xe)[1]/size(data.xe)[1]
	(bb=bb, xe=xe, effbb = eff_bb, effxe = eff_xe)
end
	

function select_blob(bb, xe;  eb2cut=400.0)
		bbsel = filter(row -> row.Eb2_keV > eb2cut, bb )
		xesel = filter(row -> row.Eb2_keV > eb2cut, xe )
		return bbsel, xesel
end


####
## plots
###


function plot_eb2_cut_eff_and_fom(bbdf, xedf;
                                  Eb2="Eb2_keV",
                                  length=100,
                                  rinf = 100,
                                  rsup = 800,
                                  fom_bias=0.1,
                                  suptitle="Blob2 Energy Cut Analysis")
	eb2_sig = bbdf[:,Eb2]
	eb2_bkg = xedf[:,Eb2]
	p1 = Petit.plot_efficiency_vs_cut(eb2_sig,
                                eb2_bkg;
                                cuts = range(rinf, rsup, length=length),
                                xlabel = "Eblob2 (keV)",
                                title = "Efficiency vs Cut",
                                signal_label = "bb0nu",
                                background_label = "Xe137")
	p2 = Petit.plot_fom_vs_cut(eb2_sig,
                         eb2_bkg;
                         cuts = range(rinf, rsup, length=length),
                         xlabel = "Eblob2 (keV)",
                         title= "FOM vs Cut",
                         mark_optimal = true,
                         bias=fom_bias)
	plot(p1, p2;
         layout=(1,2),
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

