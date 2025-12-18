function copy_eff!(effbkg, effiso)
	effbkg.efffid =effiso.efffid
	effbkg.effroi =effiso.effroi
	effbkg.effblb =effiso.effblb
end

total_eff(e) = reduce(*, (e.eff1tr, e.efffid, e.effroi, e.effblb))

function histo_b1_b2(xeb1, xeb2)
	
	h_b1, p_b1 = jn.Petit.step_hist(xeb1;
	         nbins = 20,
	         xlabel = "Eb1",
	        ylabel = "Frequency",
	         title=" Eb1 ")
	h_b2, p_b2 = jn.Petit.step_hist(xeb2;
	         nbins = 20,
	         xlabel = "Eb2",
	        ylabel = "Frequency",
	         title=" Eb2 ")
	
	plot(p_b1, p_b2)
end

function energy_trk_length(xe, xtl)
	p1 = etrk_plot(xe, nbins=50, title="Energy (keV)")
	p2 = trk_plot(xtl; nbins=30, title="Track Length (mm)")
    plot(p1, p2, 
        layout = (1, 2),
        size = (800, 500),
        plot_title = "Energy and Track Length",
        plot_titlefontsize = 16,
        margin = 5Plots.mm
    )
end

function eb1_vs_eb2(xeb1, xeb2, effdf; eblob_cut=600.0)
	p1 = eb1_eb2_plot(xeb1, xeb2, title="Bi-214", eblob_cut=eblob_cut)
	p2 = fom_plot(effdf, "Eff: Bi-214")
    plot(p1, p2,
        layout = (1, 2),
        size = (900, 500),
        #plot_title = "energy blob1 vs energy blob2",
        plot_titlefontsize = 16,
        margin = 5Plots.mm
    )
end

function eb1_eb2_plot(eblob1, eblob2; title="", eblob_cut=600.0)

    p3 = scatter(eblob1, eblob2,
        xlabel = "Blob 1 Energy (keV)",
        ylabel = "Blob 2 Energy (keV)",
        title = "Blob Energy Correlation",
        label = nothing,
        markersize = 4,
        markercolor = :purple,
        markeralpha = 0.6,
        markerstrokewidth = 0.5,
        markerstrokecolor = :black,
        xlims = (0, maximum([maximum(eblob1), maximum(eblob2)]) * 1.05),
        ylims = (0, maximum([maximum(eblob1), maximum(eblob2)]) * 1.05),
        aspect_ratio = :equal,
        grid = true,
        gridstyle = :dot,
        gridalpha = 0.3
    )

    # Add diagonal reference line
    max_energy = maximum([maximum(eblob1), maximum(eblob2)])
    plot!(p3, [0, max_energy], [0, max_energy],
        line = :dash,
        linecolor = :red,
        linewidth = 1,
        label = "y = x",
        alpha = 0.5
    )

    # Add horizontal line at eblob_cut level
    plot!(p3, [0, max_energy], [eblob_cut, eblob_cut],
        line = :dash,
        linecolor = :blue,
        linewidth = 2,
        label = "Cut: $(eblob_cut) keV",
        alpha = 0.7
    )
	p3


end

function fom_plot(results_df::DataFrame, title::String="")
    # Create the plot with error bars
    p = scatter(results_df.eblob2, results_df.eff,
                yerror = results_df.err,
                xlabel = "Blob 2 Energy Cut (keV)",
                ylabel = "Efficiency",
                title = isempty(title) ? "Efficiency vs Energy Cut" : title,
                label = "Data",
                markersize = 5,
                markercolor = :blue,
                markerstrokewidth = 1,
                markerstrokecolor = :black,
                linecolor = :blue,
                linewidth = 1,
                ylims = (0, maximum(results_df.eff) * 1.1),
                grid = true,
                gridstyle = :dot,
                gridalpha = 0.3,
                legend = :topright)

    # Add a smooth line through the points
    plot!(p, results_df.eblob2, results_df.eff,
          line = :solid,
          linewidth = 2,
          label = nothing,
          color = :blue,
          alpha = 0.5)

    # Add percentage labels on secondary y-axis
    plot!(p, yticks = (0:0.1:1.0, ["$(Int(y*100))%" for y in 0:0.1:1.0]))

    return p
end

function fom_plots(results::Vector{DataFrame}, 
				   labels::Vector{String},title::String="")
	p = scatter(
                xlabel = "Blob 2 Energy Cut (keV)",
                ylabel = "Efficiency",
                title = isempty(title) ? "Efficiency vs Energy Cut" : title,
                label = "Data",
                markersize = 5,
                markercolor = :blue,
                markerstrokewidth = 1,
                markerstrokecolor = :black,
                linecolor = :blue,
                linewidth = 1,
                grid = true,
                gridstyle = :dot,
                gridalpha = 0.3,
                legend = :topright)
	
    for (i, result_df) in enumerate(results)
    	p = scatter!(p, result_df.eblob2, result_df.eff,	 
                	yerror = result_df.err,
					ylims = (0, 0.2),
                	label = labels[i])

    	# Add a smooth line through the points
   		p = plot!(p, result_df.eblob2, result_df.eff,
          	line = :solid,
          	linewidth = 2,
          	label = nothing,
          	color = :blue,
          	alpha = 0.5)
	end

    # Add percentage labels on secondary y-axis
    #plot!(p, yticks = (0:0.1:1.0, ["$(Int(y*100))%" for y in 0:0.1:1.0]))

    return p
end