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


"""
    read_cnn_efficiencies(csvdir::String)

Read CNN efficiency CSV files and return interpolation functions for signal and background.

# Arguments
- `csvdir::String`: Directory containing the CNN efficiency CSV files

# Returns
NamedTuple with:
- `interp_1mm`: NamedTuple with (signal, background) interpolation functions for 1mm
- `interp_3p5mm`: NamedTuple with (signal, background) interpolation functions for 3.5mm
- `interp_10mm`: NamedTuple with (signal, background) interpolation functions for 10mm
- `data_1mm`, `data_3p5mm`, `data_10mm`: Raw DataFrames

Each interpolation function takes a threshold value and returns the efficiency.
"""
function read_cnn_efficiencies(csvdir::String)
    # File names
    file_1mm = joinpath(csvdir, "efficiency_data_MC_truth_1mm_15bar_214Bi.csv")
    file_3p5mm = joinpath(csvdir, "efficiency_data_MC_truth_3.5mm_15bar_214Bi.csv")
    file_10mm = joinpath(csvdir, "efficiency_data_MC_truth_10mm_15bar_214Bi.csv")

    # Read CSV files
    df_1mm = CSV.read(file_1mm, DataFrame)
    df_3p5mm = CSV.read(file_3p5mm, DataFrame)
    df_10mm = CSV.read(file_10mm, DataFrame)

    # Create interpolation functions
    # Sort by threshold (ascending) for proper interpolation
    function make_interpolators(df)
        sorted_df = sort(df, :threshold)
        thresh = Float64.(sorted_df.threshold)
        sig_eff = Float64.(sorted_df.signal_efficiency)
        bkg_eff = Float64.(sorted_df.background_efficiency)

        # Linear interpolation
        sig_interp = linear_interpolation(thresh, sig_eff, extrapolation_bc=Flat())
        bkg_interp = linear_interpolation(thresh, bkg_eff, extrapolation_bc=Flat())

        return (signal=sig_interp, background=bkg_interp)
    end

    interp_1mm = make_interpolators(df_1mm)
    interp_3p5mm = make_interpolators(df_3p5mm)
    interp_10mm = make_interpolators(df_10mm)

    return (
        interp_1mm = interp_1mm,
        interp_3p5mm = interp_3p5mm,
        interp_10mm = interp_10mm,
        data_1mm = df_1mm,
        data_3p5mm = df_3p5mm,
        data_10mm = df_10mm
    )
end


"""
    compute_roc(interp; thresholds=range(0.0, 1.0, length=100))

Compute ROC curve from signal and background interpolation functions.

# Arguments
- `interp`: NamedTuple with (signal, background) interpolation functions
- `thresholds`: Range of threshold values to evaluate (default: 0 to 1)

# Returns
NamedTuple with:
- `tpr`: True positive rate (signal efficiency)
- `fpr`: False positive rate (background efficiency)
- `thresholds`: The threshold values used
- `auc`: Area under the ROC curve
"""
function compute_roc(interp; thresholds=range(0.0, 1.0, length=100))
    thresh = collect(thresholds)
    tpr = [interp.signal(t) for t in thresh]      # True positive rate
    fpr = [interp.background(t) for t in thresh]  # False positive rate

    # Compute AUC using trapezoidal rule
    # Sort by FPR for proper integration
    sorted_idx = sortperm(fpr)
    fpr_sorted = fpr[sorted_idx]
    tpr_sorted = tpr[sorted_idx]

    auc = 0.0
    for i in 2:length(fpr_sorted)
        auc += (fpr_sorted[i] - fpr_sorted[i-1]) * (tpr_sorted[i] + tpr_sorted[i-1]) / 2
    end

    return (tpr=tpr, fpr=fpr, thresholds=thresh, auc=auc)
end


"""
    compute_all_rocs(cnn; thresholds=range(0.0, 1.0, length=100))

Compute ROC curves for all three detector resolutions.

# Arguments
- `cnn`: Result from read_cnn_efficiencies
- `thresholds`: Range of threshold values

# Returns
NamedTuple with roc_1mm, roc_3p5mm, roc_10mm (each containing tpr, fpr, thresholds, auc)
"""
function compute_all_rocs(cnn; thresholds=range(0.0, 1.0, length=100))
    roc_1mm = compute_roc(cnn.interp_1mm; thresholds=thresholds)
    roc_3p5mm = compute_roc(cnn.interp_3p5mm; thresholds=thresholds)
    roc_10mm = compute_roc(cnn.interp_10mm; thresholds=thresholds)

    return (roc_1mm=roc_1mm, roc_3p5mm=roc_3p5mm, roc_10mm=roc_10mm)
end


"""
    plot_roc(rocs; labels=["1mm", "3.5mm", "10mm"])

Plot ROC curves for multiple detector resolutions.

# Arguments
- `rocs`: Result from compute_all_rocs
- `labels`: Labels for each curve

# Returns
- Plot object
"""
function plot_roc(rocs; labels=["1mm (σ)", "3.5mm (σ)", "10mm (σ)"])
    p = plot(xlabel="True Positive Rate (Signal Eff.)",
             ylabel="Background Rejection (1 - Bkg Eff.)",
             title="ROC Curves",
             legend=:bottomleft,
             xlims=(0, 1), ylims=(0, 1),
             aspect_ratio=:equal,
             grid=true, gridstyle=:dot, gridalpha=0.3)

    # Plot diagonal (random classifier)
    plot!(p, [0, 1], [1, 0], linestyle=:dash, color=:gray,
          label="Random", alpha=0.5)

    # Plot each ROC curve (TPR on X, 1-FPR on Y)
    colors = [:blue, :green, :red]
    for (i, (roc, label, color)) in enumerate(zip(
            [rocs.roc_1mm, rocs.roc_3p5mm, rocs.roc_10mm], labels, colors))
        plot!(p, roc.tpr, 1.0 .- roc.fpr,
              label="$label (AUC=$(round(roc.auc, digits=3)))",
              linewidth=2, color=color)
    end

    return p
end


"""
    plot_efficiency_vs_threshold(cnn; thresholds=range(0.0, 1.0, length=100))

Plot signal and background efficiency vs threshold for all detector resolutions.

# Returns
- (plot_signal, plot_background): Tuple of two plots
"""
function plot_efficiency_vs_threshold(cnn; thresholds=range(0.0, 1.0, length=100))
    thresh = collect(thresholds)

    # Signal efficiency plot
    p_sig = plot(xlabel="Threshold", ylabel="Signal Efficiency",
                 title="Signal Efficiency vs Threshold",
                 legend=:topright, grid=true, ylims=(0, 1))

    plot!(p_sig, thresh, cnn.interp_1mm.signal.(thresh),
          label="1mm", linewidth=2, color=:blue)
    plot!(p_sig, thresh, cnn.interp_3p5mm.signal.(thresh),
          label="3.5mm", linewidth=2, color=:green)
    plot!(p_sig, thresh, cnn.interp_10mm.signal.(thresh),
          label="10mm", linewidth=2, color=:red)

    # Background efficiency plot
    p_bkg = plot(xlabel="Threshold", ylabel="Background Efficiency",
                 title="Background Efficiency vs Threshold",
                 legend=:topright, grid=true, ylims=(0, 1))

    plot!(p_bkg, thresh, cnn.interp_1mm.background.(thresh),
          label="1mm", linewidth=2, color=:blue)
    plot!(p_bkg, thresh, cnn.interp_3p5mm.background.(thresh),
          label="3.5mm", linewidth=2, color=:green)
    plot!(p_bkg, thresh, cnn.interp_10mm.background.(thresh),
          label="10mm", linewidth=2, color=:red)

    return (p_sig, p_bkg)
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