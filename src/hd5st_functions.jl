function copy_eff!(effbkg, effiso)
	effbkg.efffid =effiso.efffid
	effbkg.effroi =effiso.effroi
	effbkg.effblb =effiso.effblb
end

total_eff(e) = reduce(*, (e.eff1tr, e.efffid, e.effroi, e.effblb))

function histo_b1_b2(xeb1, xeb2)
	h_b1, p_b1 = step_hist(xeb1;
	         nbins = 20,
	         xlabel = "Eb1",
	        ylabel = "Frequency",
	         title=" Eb1 ")
	h_b2, p_b2 = step_hist(xeb2;
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


function gaussian_efficiency_analysis(μ1::Float64, μ2::Float64, 
									  σ1::Float64, σ2::Float64; 
									  xmin::Float64=2420.0, xmax::Float64=2500.0, 
									  nsteps::Int=20, ff::Float64=1.0)

	  println("Signal: μ1 = $(μ1), σ1  =$(σ1)")
	  println("Bkgnd:  μ2 = $(μ2), σ2  =$(σ2)")

      # Create normalized Normal distributions
      ng1 = Normal(μ1, σ1)
      ng2 = Normal(μ2, σ2)


      # Combined range for plotting both distributions
      xmin_plot = xmin
      xmax_plot = xmax
      x_plot = range(xmin_plot, xmax_plot, length=200)

      # Plot 1: Both normalized distributions
      p1 = plot(x_plot, pdf.(ng1, x_plot),
                label="Gaussian 1: μ=$(round(μ1, digits=2)), σ=$(round(σ1, digits=2))",
                lw=2, legend=:best, xlabel="x", ylabel="Probability Density",
                title="Normalized Gaussian Distributions")
      plot!(p1, x_plot, pdf.(ng2, x_plot),
            label="Gaussian 2: μ=$(round(μ2, digits=2)), σ=$(round(σ2, digits=2))",
            lw=2)

      # Compute efficiency curves (10 steps)
     
      # For Gaussian 1
      step_size = (xmax - xmin) / nsteps
      eff1 = Float64[]
      steps = Float64[]

      for i in 1:nsteps
          x_upper = xmin + i * step_size
          # Integral from x_upper to xmax divided by total integral
          integral_partial = cdf(ng1, xmax) - cdf(ng1, x_upper)
          integral_total = cdf(ng1, xmax) - cdf(ng1, xmin)
          efficiency = integral_partial / integral_total
          push!(eff1, efficiency)
          push!(steps, x_upper)
      end

      eff2 = Float64[]
      
      for i in 1:nsteps
		  x_upper = xmin + i * step_size
          # Integral from x_upper to xmax divided by total integral
          integral_partial = cdf(ng2, xmax) - cdf(ng2, x_upper)
          integral_total = cdf(ng2, xmax) - cdf(ng2, xmin)
          efficiency = integral_partial / integral_total
          push!(eff2, ff*efficiency)
          #push!(steps, x_upper)
      end

      # Plot 2: Efficiency curves
      p2 = plot(steps, eff1,
                marker=:circle, markersize=6, lw=2,
                label="Gaussian 1 Efficiency",
                xlabel="Step", ylabel="Efficiency",
                title="Efficiency Curves",
                legend=:topright, ylim=(0, 1.05))
      plot!(p2, steps, eff2,
            marker=:square, markersize=6, lw=2,
            label="Gaussian 2 Efficiency")

      # Add horizontal line at efficiency = 1
      hline!(p2, [1.0], linestyle=:dash, color=:gray, label="", alpha=0.5)

      return ng1, ng2, p1, eff1, eff2, steps, p2
  end


  function gaussian_efficiency_analysis(fg1, fg2; xmin=2420.0, xmax=2500.0, 
									  nsteps=20, ff=1.0)
      # Extract parameters
      μ1, σ1 = fg1.mu[1], fg1.std[1]
      μ2, σ2 = fg2.mu[1], fg2.std[1]
	  gaussian_efficiency_analysis(μ1, μ2, σ1, σ2; 
								   xmin=xmin, xmax=xmax, 
								   nsteps=nsteps, ff=ff)

	  
  end


function histo_xyz(X,Y,Z,sample)
	_, px = step_hist(X; 
                   nbins = 50,
				   logy=true, 
                   #xlim = (emin, emax),
                   xlabel = "X (mm)",
                   ylabel = "Frequency",
                   title ="$(sample) : X  (mm)")
	_, py = step_hist(Y;
                   nbins = 50,
				   logy=true, 
                   #xlim = (emin, emax),
                   xlabel = "Y (mm)",
                   ylabel = "Frequency",
                   title ="$(sample) :Y (mm)")
	_, pz = step_hist(Z;
                   nbins = 50,
				   logy=true, 
                   #xlim = (emin, emax),
                   xlabel = "Z (mm)",
                   ylabel = "Frequency",
                   title ="$(sample) :Z (mm)")
	
	pxy  = scatter(X, Y, xlabel="X (mm)", ylabel="Y (mm)",
                     title="$(sample) : X vs Y",
                     markersize=2,
                     alpha=0.5,
                     legend=false)
	
	 plot(px, py, pz, pxy, layout=(2,2), size=(1400, 800))
end


function etrk_plot(etrk; nbins=30, title="")

    # Create subplots
    p1 = histogram(etrk,
        bins = nbins,
        xlabel = "Track Energy (keV)",
        ylabel = "Counts",
        title = title,
        label = nothing,
        fillcolor = :steelblue,
        linecolor = :black,
        alpha = 0.7,
        xlims = (2400, 2550)
    )

    # Add mean and std annotation
    mean_etrk = mean(etrk)
    std_etrk = std(etrk)
    annotate!(p1,
        :topright,
        text(@sprintf("μ = %.1f keV\nσ = %.1f keV", mean_etrk, std_etrk), 10, :left)
    )

    p1
end


function trk_plot(trkl; nbins=30, title="")
    p1 = histogram(trkl,
        bins = nbins,
        xlabel = "Track Length (mm)",
        ylabel = "Counts",
        title = title,
        label = nothing,
        fillcolor = :darkorange,
        linecolor = :black,
        alpha = 0.7,
        xlims = (minimum(trkl) * 0.95, maximum(trkl) * 1.05)
    )

    # Add mean and std annotation
    mean_trkl = mean(trkl)
    std_trkl = std(trkl)
    annotate!(p1,
        :topright,
        text(@sprintf("μ = %.1f mm\nσ = %.1f mm", mean_trkl, std_trkl), 10, :left)
    )

    p1
end


function fom_ecut(etrkt; min_cut=2400.0, max_cut=2500.0,step=20.0)

	function simple_binomial_error(k::Int, n::Int)
	    if n == 0
	        return (0.0, 0.0)
	    end
	
	    p = k / n
	
	    # Standard error for binomial proportion
	    # σ = sqrt(p(1-p)/n)
	    error = sqrt(p * (1 - p) / n)
	
	    return (p, error)
	end

	# only consider events in range 
	etrk = in_range(etrkt, min_cut, max_cut)
	ff = sum(etrk)/sum(etrkt)
    # Total number of events
    n_total = length(etrk)


    # Initialize results arrays
    cuts = Float64[]
    efficiencies = Float64[]
    errors = Float64[]

    # Iterate over cut values
    cut_value = min_cut
    while cut_value <= max_cut
        # Count events passing the cut
        n_pass = sum(etrk .>= cut_value)

        # Calculate efficiency and error
        
        eff, err = simple_binomial_error(n_pass, n_total)
        
        push!(cuts, cut_value)
        push!(efficiencies, eff*ff)
        push!(errors, err)

        cut_value += step
    end

    # Create output dataframe
    results_df = DataFrame(
        eblob2 = cuts,
        eff = round.(efficiencies, digits=4),
        err = round.(errors, digits=4)
    )

    return results_df, ff
end