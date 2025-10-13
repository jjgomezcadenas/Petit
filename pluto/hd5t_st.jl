### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 349825ff-7ffe-4fa1-ba26-a772041f0323
begin
	using Revise
	using PlutoUI
	using CSV
	using DataFrames
	using Plots
	using Printf
	using HDF5
	using Markdown
	using InteractiveUtils
	using Statistics
	using StatsBase
	using Distributions
	using StatsPlots
	#using DataFramesMeta
	
	import Glob
	#using Interpolations
	#using QuadGK
	#using LsqFit
	#using Chain
	#using Unitful 
end

# ╔═╡ 04b446d6-f34f-11ed-2565-0b15d65b6781
PlutoUI.TableOfContents(title="HD5t analysis", indent=true)


# ╔═╡ 871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
begin
	cmdir=joinpath(ENV["DATA"], "HD5t/precdr")
	pdir =joinpath(ENV["PROJECTS"], "Petit")
end

# ╔═╡ 947c237c-9852-40e9-a83f-c23666db90aa
begin
      using Pkg
      Pkg.activate(pdir)
      Pkg.instantiate()
  end

# ╔═╡ 7504d7aa-a780-4956-99a5-08a7f9a462b2
function ingredients(path::String)
    # this is from the Julia source code (evalfile in base/loading.jl)
    # but with the modification that it returns the module instead of the last object
    name = Symbol(basename(path))
    m = Module(name)
    Core.eval(m,
        Expr(:toplevel,
                :(eval(x) = $(Expr(:core, :eval))($name, x)),
                :(include(x) = $(Expr(:top, :include))($name, x)),
                :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
                :(include($path))))
    m
end

# ╔═╡ c9fc0547-0e73-4629-9909-e59c3d75169d
begin
	jn = ingredients(string(pdir,"/src/Petit.jl"))
end

# ╔═╡ 8290afe2-e5cb-4793-b108-a313f2677119
md"""
## Functions
"""

# ╔═╡ 982a3890-7318-4034-b7fd-decafc288362
function voxel_size(dx, dy, p, l)
	sqrt(l)*(dx + dy) /(2.0*sqrt(p))
end

# ╔═╡ b581c086-865c-4761-864d-ca8e549f8a75
function copy_eff!(effbkg, effiso)
	effbkg.efffid =effiso.efffid
	effbkg.effroi =effiso.effroi
	effbkg.effblb =effiso.effblb
end


# ╔═╡ 0d6dfea6-48ea-4645-b00c-568ea4e204af
total_eff(e) = reduce(*, (e.eff1tr, e.efffid, e.effroi, e.effblb))

# ╔═╡ 46e73d58-dc53-42bb-add7-088ac144ca24
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

# ╔═╡ 65ec723a-6a38-4a2c-88e5-d6e1c59375bf
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

# ╔═╡ 1052ab37-288a-4c84-9bbe-3407a6a0dfaf
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


# ╔═╡ d2b5f0e9-4a79-4c33-afb2-7e970fb5565d
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

# ╔═╡ 641e9a21-5999-4eaa-8904-8e03fa952b49
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
					#ylims = (0, maximum(result_df.eff) * 1.1),
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
    plot!(p, yticks = (0:0.1:1.0, ["$(Int(y*100))%" for y in 0:0.1:1.0]))

    return p
end


# ╔═╡ d61ec9b7-4922-45de-9b6f-e0a9f4386740
md"""
## Analysis
"""

# ╔═╡ f4e98f28-5f3a-49f5-a053-8ed532e92a60
begin
	vhe = voxel_size(0.75, 1.6, 15.0, 100.0)
	md"""
	- voxel size for 10 % He = $@sprintf("%.1e", vhe) mm
	"""
end

# ╔═╡ 8fc7e228-fb0d-4c6c-ace9-8f6cb057a989
begin
	erex = 12.5 # keV
	rblob = 10.0 # mm
	i = 1
	md"""
	- Energy resolution (keV) = $(erex)
	- rblob (mm) = $(rblob)
	- example track to inspect number = $(i)
	"""
end

# ╔═╡ 22c204c8-9e36-41e0-a0d4-0b86c689739c
begin
	confCut= 0.40
	trklCut = 10.0
	eblb2Max = 1000.0
	md"""
	- Confidence cut = $(confCut)
	- Track Length  cut = $(trklCut)
	- Max energy of Blob 2 = $(eblb2Max)
	"""
end

# ╔═╡ de068020-a64c-43b4-9e78-8ade8a4f6274
begin
	fnBi = "blobsBi.csv"
	fnTl = "blobsTl.csv"
	fnBB = "blobsBB.csv"
	fnXe = "blobsXe.csv"
	md"""
	- file blobs Bi214 = $(fnBi)
	- file blobs Tl208 = $(fnTl)
	- file blobs bb0nu = $(fnBB)
	- file blobs Xe-137 = $(fnXe)
	"""
end

# ╔═╡ b5f850b0-19ec-4e96-96cc-4c43a91e05f3
begin
	trksBiCopper = jn.Petit.get_bkgnd_tracks(cmdir, isotope="bi214", 
						             bkgnd="copperbkg", 
						             tag="*st3mm*")
	effBiCopper = jn.Petit.Eff(trksBiCopper.n1trk/trksBiCopper.ntot, 1.0, 1.0, 1.0)

	trksBiInner = jn.Petit.get_bkgnd_tracks(cmdir, isotope="bi214", 
						             bkgnd="innerbkg", 
						             tag="*st3mm*")
	effBiInner = jn.Petit.Eff(trksBiInner.n1trk/trksBiInner.ntot, 1.0, 1.0, 1.0)
	
	md"""
	### Bi 214
	#### Copper
	- Total number of events generated = $@sprintf("%.2e", trksBiCopper.ntot)  
	- Total number of events 1 trk  = $(trksBiCopper.n1trk)
	- Selection efficiency: $@sprintf("%.2e", effBiCopper.eff1tr) 

	#### Inner
	- Total number of events generated = $@sprintf("%.2e", trksBiInner.ntot)  
	- Total number of events 1 trk  = $(trksBiInner.n1trk)
	- Selection efficiency: $@sprintf("%.2e", effBiInner.eff1tr) 
	"""
end

# ╔═╡ a60ea1b2-1c4c-4565-ac34-d2580dd016e3
md"""
### Build Tracks: example
"""

# ╔═╡ 99258d62-017c-45e2-b9e5-63da00815166
md"""
#### Compute Eblobs for all tracks
"""

# ╔═╡ 4f254ee3-c000-4c49-9f6b-0d7ef947c67d
md"Run blob analysis for Bi214? $(@bind babi CheckBox(default=false))"

# ╔═╡ bc7634bf-2d76-4966-95c1-70d3ea97e0ff
blobsBi = CSV.read(fnBi, DataFrame)

# ╔═╡ 067d7785-ce69-4cb6-b783-9fbeab4b9f32
begin
	hEnergyBiTrue, pEnergyBiTrue = jn.Petit.step_hist(collect(blobsBi.energyKeV);
	                                                 nbins = 40,
										             xlim   = (2300.0, 2700.0),
	                                                 xlabel = " E (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" E true (keV)")
	
	_, pConfidence = jn.Petit.step_hist(collect(blobsBi.confidence);
	                                                 nbins = 20,
										             xlim   = (0.0, 1.0),
	                                                 xlabel = " Confidence ",
	                                                 ylabel = "Frequency",
	                                                 title=" Fit confidence")
	_, pTrkL = jn.Petit.step_hist(collect(blobsBi.trackLength);
	                                                 nbins = 40,
										             #xlim   = (2300.0, 2700.0),
	                                                 xlabel = " Track Length (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" E true (keV)")
	
	eBi = jn.Petit.smear_histogram(hEnergyBiTrue, erex)
	
	hEnergyBi, pEnergyBi = jn.Petit.step_hist(eBi;
                                              nbins = 40,
                                              xlim   = (2300.0, 2700.0),
                                              xlabel = " E (keV)",
                                              ylabel = "Frequency",
                                              title=" Bi 214 E ")
	plot(pEnergyBiTrue, pEnergyBi, pConfidence, pTrkL )
end

# ╔═╡ c5374971-03ed-4bec-88b6-e636ac400348
md"""
#### Clean up 
- Reject tracks with confidence > $(confCut) 
- Reject tracks with track length < $(trklCut) 
- Reject tracks with track E blob2 > $(eblb2Max) 
"""

# ╔═╡ 73cb4e71-d232-4c97-ba8f-27d170823073
begin
	fgbi, pfgbi =jn.Petit.plot_fit_gauss(eBi, "Track Energy", "Counts",
                        30, 2400.0, 2500.0;
                        xgmin=2400.0, xgmax=2500.0, gbins=30)
	plot(pfgbi)
end

# ╔═╡ 046ec784-bc63-4dde-920d-36fc85f1ee8e


# ╔═╡ 48b32bbf-abbb-480e-9ba0-0df5f0777336
begin
	trksTlCopper = jn.Petit.get_bkgnd_tracks(cmdir, isotope="tl208", 
						             bkgnd="copperbkg", 
						             tag="*st3mm*")
	effTlCopper = jn.Petit.Eff(trksTlCopper.n1trk/trksTlCopper.ntot, 1.0, 1.0, 1.0)

	trksTlInner = jn.Petit.get_bkgnd_tracks(cmdir, isotope="tl208", 
						             bkgnd="innerbkg", 
						             tag="*st3mm*")
	effTlInner = jn.Petit.Eff(trksTlInner.n1trk/trksTlInner.ntot, 1.0, 1.0, 1.0)
	
	md"""
	### Tl 208
	#### Copper
	- Total number of events generated = $@sprintf("%.2e", trksTlCopper.ntot)  
	- Total number of events 1 trk  = $(trksTlCopper.n1trk)
	- Selection efficiency: $@sprintf("%.2e", effTlCopper.eff1tr) 

	#### Inner
	- Total number of events generated = $@sprintf("%.2e", trksTlInner.ntot)  
	- Total number of events 1 trk  = $(trksTlInner.n1trk)
	- Selection efficiency: $@sprintf("%.2e", effTlInner.eff1tr) 
	"""
end

# ╔═╡ 9cd7ccf0-a4e7-4445-a6b1-78ffa573eb3b
md"Run blob analysis for Bi214? $(@bind batl CheckBox(default=false))"


# ╔═╡ 37a7d81a-ca8e-4ac5-b2ca-04b261beaddb
blobsTl = CSV.read(fnTl, DataFrame)

# ╔═╡ 01297be8-6234-4bae-91d5-963891b755ef
begin
	hEnergyTlTrue, pEnergyTlTrue = jn.Petit.step_hist(collect(blobsTl.energyKeV);
	                                                 nbins = 40,
										             xlim   = (2300.0, 2700.0),
	                                                 xlabel = " E (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" E true (keV)")
	
	_, pConfidence2 = jn.Petit.step_hist(collect(blobsTl.confidence);
	                                                 nbins = 20,
										             xlim   = (0.0, 1.0),
	                                                 xlabel = " Confidence ",
	                                                 ylabel = "Frequency",
	                                                 title=" Fit confidence")
	_, pTrkL2 = jn.Petit.step_hist(collect(blobsTl.trackLength);
	                                                 nbins = 40,
										             #xlim   = (2300.0, 2700.0),
	                                                 xlabel = " Track Length (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" E true (keV)")
	
	eTl = jn.Petit.smear_histogram(hEnergyTlTrue, erex)
	
	hEnergyTl, pEnergyTl = jn.Petit.step_hist(eTl;
                                              nbins = 40,
                                              xlim   = (2300.0, 2700.0),
                                              xlabel = " E (keV)",
                                              ylabel = "Frequency",
                                              title=" Tl 208 E ")
	plot(pEnergyTlTrue, pEnergyTl, pConfidence2, pTrkL2 )
end

# ╔═╡ 2e1c70e9-6ecf-4e97-8705-6342028285aa
md"""
## bb0nu
"""

# ╔═╡ a6e42cf4-1389-497f-90da-1df8297cec74
begin
	trksbb = jn.Petit.get_bb0nu_tracks(cmdir)
	#trksbb.ntot = trksbb.metas[1]["events_processed"]

	effBB = jn.Petit.Eff(trksbb.n1trk/trksbb.ntot, 1.0, 1.0, 1.0)
	
	md"""
	- Total number of bb0nu events generated = $(trksbb.ntot)
	- Total number of bb0nu events 1 trk  = $(trksbb.n1trk)
	- Selection efficiency: $@sprintf("%.2e", effBB.eff1tr) 
	"""
end

# ╔═╡ f7b87c23-dee5-462d-b4b9-b964947b9ca7
md"Run blob analysis for bb0nu? $(@bind babb CheckBox(default=false))"


# ╔═╡ 18d8eaa0-5467-4123-bfc8-f52f377e7c72
blobsBB = CSV.read(fnBB, DataFrame)

# ╔═╡ ac45ca81-79bc-4631-bf37-c08bbbf6159e
begin
	hEnergyBBTrue, pEnergyBBTrue = jn.Petit.step_hist(collect(blobsBB.energyKeV);
	                                                 nbins = 40,
										             xlim   = (2300.0, 2700.0),
	                                                 xlabel = " E (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" E true (keV)")
	
	_, pConfidence3 = jn.Petit.step_hist(collect(blobsBB.confidence);
	                                                 nbins = 20,
										             xlim   = (0.0, 1.0),
	                                                 xlabel = " Confidence ",
	                                                 ylabel = "Frequency",
	                                                 title=" Fit confidence")
	_, pTrkL3 = jn.Petit.step_hist(collect(blobsBB.trackLength);
	                                                 nbins = 40,
										             #xlim   = (2300.0, 2700.0),
	                                                 xlabel = " Track Length (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" E true (keV)")
	
	eBB = jn.Petit.smear_histogram(hEnergyBBTrue, erex)
	
	hEnergyBB, pEnergyBB = jn.Petit.step_hist(eBB;
                                              nbins = 40,
                                              xlim   = (2300.0, 2700.0),
                                              xlabel = " E (keV)",
                                              ylabel = "Frequency",
                                              title=" Bi 214 E ")
	plot(pEnergyBBTrue, pEnergyBB, pConfidence3, pTrkL3 )
end

# ╔═╡ 6ea54a0d-8b35-47d9-ac29-43b8f1d103db
begin
	fgBB, pfgBB =jn.Petit.plot_fit_gauss(eBB, "Track Energy", "Counts",
                        30, 2400.0, 2500.0;
                        xgmin=2400.0, xgmax=2500.0, gbins=30)
	plot(pfgBB)
end


# ╔═╡ ffb9d0d1-bca6-4b9a-87a7-e8c6280488fa
md"""
## Signal and background efficiency
"""

# ╔═╡ 0397d23d-6c74-484b-9c8a-3989bab21239
md"""
## Xe-137
"""

# ╔═╡ d7729d46-79ac-48c2-948a-41c91f78e2e3
begin
	trksxe =  jn.Petit.get_xe137_tracks(cmdir)
	effXe = jn.Petit.Eff(trksxe.n1trk/trksxe.ntot, 1.0, 1.0, 1.0)
	
	md"""
	- Total number of Xe-137 events generated = $(trksxe.ntot)
	- Total number of Xe-137 events 1 trk  = $(trksxe.n1trk)
	- Selection efficiency: $@sprintf("%.2e", effXe.eff1tr) 
	"""
end

# ╔═╡ 2f87d310-8d44-4649-9dfb-440e160a24f4
md"Run blob analysis for xe137? $(@bind baxe CheckBox(default=false))"


# ╔═╡ 10a4a567-6337-418b-a70e-97647af0177b
blobsXe = CSV.read(fnXe, DataFrame)

# ╔═╡ f333b7c7-bebb-47da-bbfc-b350f0f34a0d
begin
	hEnergyXeTrue, pEnergyXeTrue = jn.Petit.step_hist(collect(blobsXe.energyKeV);
	                                                 nbins = 40,
										             xlim   = (2300.0, 2700.0),
	                                                 xlabel = " E (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" E true (keV)")
	
	_, pConfidence4 = jn.Petit.step_hist(collect(blobsXe.confidence);
	                                                 nbins = 20,
										             xlim   = (0.0, 1.0),
	                                                 xlabel = " Confidence ",
	                                                 ylabel = "Frequency",
	                                                 title=" Fit confidence")
	_, pTrkL4 = jn.Petit.step_hist(collect(blobsXe.trackLength);
	                                                 nbins = 40,
										             #xlim   = (2300.0, 2700.0),
	                                                 xlabel = " Track Length (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" E true (keV)")
	
	eXe = jn.Petit.smear_histogram(hEnergyXeTrue, erex)
	
	hEnergyXe, pEnergyXe = jn.Petit.step_hist(eXe;
                                              nbins = 40,
                                              xlim   = (2300.0, 2700.0),
                                              xlabel = " E (keV)",
                                              ylabel = "Frequency",
                                              title=" Xe-136E ")
	plot(pEnergyXeTrue, pEnergyXe, pConfidence4, pTrkL4 )
end

# ╔═╡ f15ad4b6-873a-4ce3-8321-5a38e9704df2
md"""
### Blob efficiencies
"""

# ╔═╡ d0f9e26b-2a3c-42a4-a826-7c8694f5d470
md"""
# Functions
"""

# ╔═╡ 8756b4ca-011b-4a6a-aff6-22488c8ac854
function eff_ecut(effbb, effbkg, ffbkg)
	effbb = effe0b.eff
	effbi = effbkg.eff
	energies = effe0b.eblob2
	
	d = Poisson.(ffbkg*effbi)
	fom  = effbb./std.(d) 

	pfom = plot(energies, fom,
                marker=:circle, markersize=5, lw=2,
                label="FOM ",
                xlabel="Step", ylabel="fom",
                title="FOM Curve",
                legend=:topright, ylim=(0, 1.5))

	return fom, pfom

end

# ╔═╡ f0ace932-f1d1-4e2b-9a6d-88b2cea4d849
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

# ╔═╡ 590f5634-843c-4589-ae4d-298293581732
function gaussian_efficiency_analysis(fg1, fg2; xmin=2420.0, xmax=2500.0, 
									  nsteps=20, ff=1.0)
      # Extract parameters
      μ1, σ1 = fg1.mu[1], fg1.std[1]
      μ2, σ2 = fg2.mu[1], fg2.std[1]
	  gaussian_efficiency_analysis(μ1, μ2, σ1, σ2; 
								   xmin=xmin, xmax=xmax, 
								   nsteps=nsteps, ff=ff)

	  
  end

# ╔═╡ 9bae1645-1485-411d-a6d3-61906c7c4194
function find_track_extremes(trk; i=1)
	xresult = jn.Petit.walk_track_from_extremes(trk[i])
	xstart_voxel, xend_voxel = xresult.extremes
  	xtrack_length = xresult.total_length
	energy_kev = 1e+3 * sum(trk[i].voxels.energy)
	return xresult, xstart_voxel, xend_voxel, xtrack_length, energy_kev
end

# ╔═╡ ffce2732-1271-4e7e-b3aa-4000f0355dff
function get_bi214_tracks()
	files = Glob.glob("*st3mm.h5", "/Users/jjgomezcadenas/Data/HD5t/precdr/copperbkg/bi214/")
	tracks, metadatas = jn.Petit.chain_track_files(files)
	nxbi214 = 0
	for i in 1:length(metadatas)
		nxbi214= nxbi214 + metadatas[i]["nevents_from_config"]
	end
	return tracks, metadatas, nxbi214, length(tracks)
end

# ╔═╡ 7dad709c-17e9-4738-b6c0-6d0c9bd70f09
function get_tl208_tracks(;bkgnd::String ="copperbkg", tag::String ="*st3mm.b*")

	path = joinpath(cmdir, bkgnd, "tl208") 
	files = Glob.glob(tag, path)
	tracks, metadatas = jn.Petit.chain_track_files(files)
	
	ntot = 0
	for i in 1:length(metadatas)
		ntot+= metadatas[i]["nevents_from_config"]
	end
	
	return tracks, metadatas, ntot, length(tracks)
end

# ╔═╡ 4f78d979-f707-44f9-892d-e278937c7c7e
cmdir

# ╔═╡ 883c8a64-9cc2-4a89-95d5-4e906a903954
function get_xe137_tracks()
	files = Glob.glob("*st3mm.b*", "/Users/jjgomezcadenas/Data/HD5t/precdr/xe137/")
	tracks, metadatas = jn.Petit.chain_track_files(files)
	
	ntl208= metadatas[1]["nevents_from_config"]
	## NB ntl208 takes only one file because we are splitting a very
	## large file in small segments, unlike the case for Bi214
	return tracks, metadatas, ntl208, length(tracks)
end

# ╔═╡ 75f229f8-1a11-4b9a-9ef6-331e98177caa
function get_tracks(path; voxel="3mm")
	str = "*st$(voxel).h5"
	files = Glob.glob(str, path)
	tracks, metadatas = jn.Petit.chain_track_files(files)
	nxbi214 = 0
	for i in 1:length(metadatas)
		nxbi214= nxbi214 + metadatas[i]["nevents_from_config"]
	end
	return tracks, metadatas, nxbi214, length(tracks)
end

# ╔═╡ 32f432a9-68c1-4a65-9b98-8f0174bf08ae


# ╔═╡ b59257e9-9a5e-467f-bbc3-f04cdf0e9a75
function histo_energy_trkl_conf(xcon, xtl, xe)
	h_c, p_c = jn.Petit.step_hist(xcon *1.0;
	         nbins = 20,
	        xlim   = (0.0, 1.0),
	         xlabel = "confidence",
	        ylabel = "Frequency",
	         title=" confidence ")
	h_tl, p_tl = jn.Petit.step_hist(xtl;
	         nbins = 20,
	        #xlim   = (0.0, 1.0),
	         xlabel = " Track Length",
	        ylabel = "Frequency",
	         title=" Track Length ")
	h_ek, p_ek = jn.Petit.step_hist(xe;
	         nbins = 20,
	        #xlim   = (0.0, 1.0),
	         xlabel = " E (keV)",
	        ylabel = "Frequency",
	         title=" E (keV)")
	return h_ek, p_ek, p_c, p_tl
end

# ╔═╡ a8d1a7d5-bfaa-4430-abde-01b167249f29


# ╔═╡ 98033f21-5720-4d22-94a8-a72e650a860b

	


# ╔═╡ 94e6ae79-2d97-461a-8a3a-bf100bb1771e
function histo_xyz(X,Y,Z,sample)
	_, px = jn.Petit.step_hist(X; 
                   nbins = 50,
				   logy=true, 
                   #xlim = (emin, emax),
                   xlabel = "X (mm)",
                   ylabel = "Frequency",
                   title ="$(sample) : X  (mm)")
	_, py = jn.Petit.step_hist(Y;
                   nbins = 50,
				   logy=true, 
                   #xlim = (emin, emax),
                   xlabel = "Y (mm)",
                   ylabel = "Frequency",
                   title ="$(sample) :Y (mm)")
	_, pz = jn.Petit.step_hist(Z;
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

# ╔═╡ 2de993fe-47dc-459d-8000-708d85ed425f


# ╔═╡ a1bac832-78d3-4ff5-8e05-bc8fc0eeef91
function histo_trak_energy(etrk, title; emin=2300.0, emax=2700.0)
	_, pe = jn.Petit.step_hist(etrk;
                   nbins = 30,
				   logy=true, 
                   xlim = (emin, emax),
                   xlabel = "energy (keV)",
                   ylabel = "Frequency",
                   title =title)
	plot(pe)
end

# ╔═╡ 0d4ea311-1602-43eb-8f57-f836434f4081
function track_positions(tracks)
	X = Float64[]
	Y = Float64[]
	Z = Float64[]
	for i in 1:length(tracks)
		append!(X, tracks[i].voxels.x)
		append!(Y, tracks[i].voxels.y)
		append!(Z, tracks[i].voxels.z)
	end
	X,Y,Z
end

# ╔═╡ 9bd589eb-0bbd-4964-b4ba-8a6f8cb7b38b
begin
	trkb214 = vcat(trksBiCopper.tracks, trksBiInner.tracks)
	ntrBi = jn.Petit.NTRKS(length(trkb214), 1, 1, 1) 
	effBi = jn.Petit.Eff(1.0, 1.0, 1.0, 1.0)
	Xbi, Ybi, Zbi =track_positions(trkb214)
	histo_xyz(Xbi,Ybi,Zbi,"Bi214")
	#ebi = track_energies_keV(trkb214)
	#histo_trak_energy(ebi, "Bi214-Energy")
end

# ╔═╡ ce0c0477-0b99-400a-a85e-f8e6cad5e094
let
	xresult, xstart_voxel, xend_voxel, xtrack_length, energy_kev=find_track_extremes(trkb214, i=i)
	blobs = jn.Petit.energy_in_spheres_around_extremes(trkb214[i], xresult, rblob)
	eb1 = blobs.blob1_energy * 1e+3
	eb2 = blobs.blob2_energy * 1e+3
	nb1 = blobs.blob1_voxel_count
	nb2 = blobs.blob2_voxel_count
	md"""
	#### Find blobs: Example 
	- confidence = $(xresult.confidence)
	- start voxel: x = $(xstart_voxel.x), y = $(xstart_voxel.y), z = $(xstart_voxel.z)
	- end voxel: x = $(xend_voxel.x), y = $(xend_voxel.y), z = $(xend_voxel.z)
	- track length L =$(xtrack_length)
	- energy = $(energy_kev) KeV
	- blob 1 energy = $(round(eb1, digits=1)) keV
	- blob 2 energy = $(round(eb2, digits=1)) keV
	- blob 1 # of voxels = $(nb1)
	- blob 2 # of voxels = $(nb2)
	"""
end

# ╔═╡ 1c202304-d6ca-4dde-a7ac-f5cb2b961b93
jn.Petit.plot_track_blobs(trkb214[i], rblob;
                         markersize_voxels=3.0,
                         show_connections=true,
                         alpha_connections=0.2,
                         alpha_spheres=0.3,
                         sphere_resolution=20)

# ╔═╡ 403b3c7f-d7c3-4762-8605-5d24fb087415
if babi
	let
		blobsBi = jn.Petit.blob_analysis(trkb214, rblob)
		jn.Petit.write_blobs_csv(blobsBi, fnBi)
	end
end

# ╔═╡ 898274d5-9add-4b0b-bf04-724bc6e06074
begin
  blobsBiFid = blobsBi[(blobsBi.confidence .> confCut) .& (blobsBi.trackLength .> trklCut) .& (blobsBi.eB2 .< eblb2Max), :]
	
	ntrBi.nfid = size(blobsBiFid)[1]
	effBi.efffid = ntrBi.nfid/ntrBi.ntot
end



# ╔═╡ d4bba489-68d1-4a21-994b-aebbb50fe601
begin
	histo_b1_b2(blobsBiFid.eB1, blobsBiFid.eB2)
end

# ╔═╡ 10ed8e9a-e5f4-4652-9f1d-8e59eb605c95
begin
	trktl208 = vcat(trksTlCopper.tracks, trksTlInner.tracks)
	ntrTl = jn.Petit.NTRKS(length(trktl208), 1, 1, 1) 
	effTl = jn.Petit.Eff(1.0, 1.0, 1.0, 1.0)
	XTl, YTl, ZTl =track_positions(trktl208)
	histo_xyz(XTl,YTl,ZTl,"Tl208")
end

# ╔═╡ ba010d8d-f836-4517-b9e3-3f154ce88497
let
	xresult, xstart_voxel, xend_voxel, xtrack_length, energy_kev=find_track_extremes(trktl208, i=i)
	blobs = jn.Petit.energy_in_spheres_around_extremes(trktl208[i], xresult, rblob)
	eb1 = blobs.blob1_energy * 1e+3
	eb2 = blobs.blob2_energy * 1e+3
	nb1 = blobs.blob1_voxel_count
	nb2 = blobs.blob2_voxel_count
	md"""
	#### Find blobs: Example 
	- confidence = $(xresult.confidence)
	- start voxel: x = $(xstart_voxel.x), y = $(xstart_voxel.y), z = $(xstart_voxel.z)
	- end voxel: x = $(xend_voxel.x), y = $(xend_voxel.y), z = $(xend_voxel.z)
	- track length L =$(xtrack_length)
	- energy = $(energy_kev) KeV
	- blob 1 energy = $(round(eb1, digits=1)) keV
	- blob 2 energy = $(round(eb2, digits=1)) keV
	- blob 1 # of voxels = $(nb1)
	- blob 2 # of voxels = $(nb2)
	"""
end

# ╔═╡ e32cb510-5628-440b-a3ed-a02c0f27005e
jn.Petit.plot_track_blobs(trktl208[i], rblob;
                         markersize_voxels=3.0,
                         show_connections=true,
                         alpha_connections=0.2,
                         alpha_spheres=0.3,
                         sphere_resolution=20)

# ╔═╡ 4bcf1a8b-7101-41e9-bb30-9874275c20a6
if batl
	let
		blobsTl = jn.Petit.blob_analysis(trktl208, rblob)
		jn.Petit.write_blobs_csv(blobsTl, fnTl)
	end
end

# ╔═╡ 9d7572dd-1ebe-426d-a8f9-3dd540775733
begin
  blobsTlFid = blobsTl[(blobsTl.confidence .> confCut) .& (blobsTl.trackLength .> trklCut) .& (blobsTl.eB2 .< eblb2Max), :]
	
	ntrTl.nfid = size(blobsTlFid)[1]
	effTl.efffid = ntrTl.nfid/ntrTl.ntot
end

# ╔═╡ aa3c48c6-3c83-47c3-ac59-ce12de38c371
begin
	histo_b1_b2(blobsTlFid.eB1, blobsTlFid.eB2)
end

# ╔═╡ 4eeabfad-fb46-4d90-a586-aaf712c9eb47
let
	copy_eff!(effBiCopper , effBi)
	copy_eff!(effBiInner , effBi)
	copy_eff!(effTlCopper , effTl)
	copy_eff!(effTlInner , effTl)
	#println("effBiCopper ->",effBiCopper)
	#println("effBiInner ->",effBiInner)
	#println("effTlCopper ->",effTlCopper)
	#println("effTlInner ->",effTlInner)

	eff_bi_copper = total_eff(effBiCopper)
	eff_bi_inner = total_eff(effBiInner)
	eff_tl_copper = total_eff(effTlCopper)
	eff_tl_inner = total_eff(effTlInner)
	eff_bb0 = total_eff(effBB)
	md"""
	#### Efficiencies for signal and background:
	- signal = $@sprintf("%.2e", eff_bb0)
	- bi from copper = $@sprintf("%.2e", eff_bi_copper)
	- bi from inner = $@sprintf("%.2e", eff_bi_inner)
	- tl from copper = $@sprintf("%.2e", eff_tl_copper)
	- tl from inner = $@sprintf("%.2e", eff_tl_inner)
	"""
end

# ╔═╡ d7648894-3e02-47ca-af42-a288190f699d
let
	copy_eff!(effBiCopper , effBi)
	copy_eff!(effBiInner , effBi)
	copy_eff!(effTlCopper , effTl)
	copy_eff!(effTlInner , effTl)
	

	eff_bi_copper = total_eff(effBiCopper)
	eff_bi_inner = total_eff(effBiInner)
	eff_tl_copper = total_eff(effTlCopper)
	eff_tl_inner = total_eff(effTlInner)
	eff_xe = total_eff(effXe)
	eff_bb0 = total_eff(effBB)
	md"""
	#### Efficiencies for signal and background:
	- signal = $@sprintf("%.2e", eff_bb0)
	- Xe137 = $@sprintf("%.2e", eff_xe)
	- bi from copper = $@sprintf("%.2e", eff_bi_copper)
	- bi from inner = $@sprintf("%.2e", eff_bi_inner)
	- tl from copper = $@sprintf("%.2e", eff_tl_copper)
	- tl from inner = $@sprintf("%.2e", eff_tl_inner)
	"""
end

# ╔═╡ 759c3202-b6b7-45e7-ab11-75ed778273f7
begin
	trksbb0 = trksbb.tracks
	ntrBB = jn.Petit.NTRKS(length(trksbb0), 1, 1, 1) 
	Xbb, Ybb, Zbb =track_positions(trksbb0)
	histo_xyz(Xbb,Ybb,Zbb,"bb0nu")
end

# ╔═╡ f398f371-5271-4af7-965b-73a7027eff0b
let
	xresult, xstart_voxel, xend_voxel, xtrack_length, energy_kev=find_track_extremes(trksbb0, i=i)
	blobs = jn.Petit.energy_in_spheres_around_extremes(trksbb0[i], xresult, rblob)
	eb1 = blobs.blob1_energy * 1e+3
	eb2 = blobs.blob2_energy * 1e+3
	nb1 = blobs.blob1_voxel_count
	nb2 = blobs.blob2_voxel_count
	md"""
	#### Find blobs: Example 
	- confidence = $(xresult.confidence)
	- start voxel: x = $(xstart_voxel.x), y = $(xstart_voxel.y), z = $(xstart_voxel.z)
	- end voxel: x = $(xend_voxel.x), y = $(xend_voxel.y), z = $(xend_voxel.z)
	- track length L =$(xtrack_length)
	- energy = $(energy_kev) KeV
	- blob 1 energy = $(round(eb1, digits=1)) keV
	- blob 2 energy = $(round(eb2, digits=1)) keV
	- blob 1 # of voxels = $(nb1)
	- blob 2 # of voxels = $(nb2)
	"""
end

# ╔═╡ 00910312-2f5c-49bf-abda-11a6acf86249
jn.Petit.plot_track_blobs(trksbb0[i], rblob;
                         markersize_voxels=3.0,
                         show_connections=true,
                         alpha_connections=0.2,
                         alpha_spheres=0.3,
                         sphere_resolution=20)

# ╔═╡ 131e6669-5a26-4d9c-9ef4-911aa779a20e
if babb
	let
		blobsBB = jn.Petit.blob_analysis(trksbb0, rblob)
		jn.Petit.write_blobs_csv(blobsBB, fnBB)
	end
end


# ╔═╡ cf3042b0-d5f3-4716-ab25-6ca736b491db
begin
  blobsBBFid = blobsBB[(blobsBB.confidence .> confCut) .& (blobsBB.trackLength .> trklCut) .& (blobsBB.eB2 .< eblb2Max), :]
	
	ntrBB.nfid = size(blobsBBFid)[1]
	effBB.efffid = ntrBB.nfid/ntrBB.ntot
end

# ╔═╡ 794ac8f8-3eee-4a35-96b0-b5dca3615a07
begin
	histo_b1_b2(blobsBBFid.eB1, blobsBBFid.eB2)
end


# ╔═╡ ac25b86e-97b3-4d54-b71a-a68a79a90995
begin
	trksxe1e = trksxe.tracks
	ntrXe = jn.Petit.NTRKS(length(trksxe1e), 1, 1, 1) 
	Xxe, Yxe, Zxe =track_positions(trksxe1e)
	histo_xyz(Xbb,Ybb,Zbb,"Xe-137")
end

# ╔═╡ 6ec7bab6-27c1-450d-a346-9abe523f1ca7
let
	xresult, xstart_voxel, xend_voxel, xtrack_length, energy_kev=find_track_extremes(trksxe1e, i=i)
	blobs = jn.Petit.energy_in_spheres_around_extremes(trksxe1e[i], xresult, rblob)
	eb1 = blobs.blob1_energy * 1e+3
	eb2 = blobs.blob2_energy * 1e+3
	nb1 = blobs.blob1_voxel_count
	nb2 = blobs.blob2_voxel_count
	md"""
	#### Find blobs: Example 
	- confidence = $(xresult.confidence)
	- start voxel: x = $(xstart_voxel.x), y = $(xstart_voxel.y), z = $(xstart_voxel.z)
	- end voxel: x = $(xend_voxel.x), y = $(xend_voxel.y), z = $(xend_voxel.z)
	- track length L =$(xtrack_length)
	- energy = $(energy_kev) KeV
	- blob 1 energy = $(round(eb1, digits=1)) keV
	- blob 2 energy = $(round(eb2, digits=1)) keV
	- blob 1 # of voxels = $(nb1)
	- blob 2 # of voxels = $(nb2)
	"""
end

# ╔═╡ 75547508-ba93-4df5-af81-6026f88499a5
jn.Petit.plot_track_blobs(trksxe1e[i], rblob;
                         markersize_voxels=3.0,
                         show_connections=true,
                         alpha_connections=0.2,
                         alpha_spheres=0.3,
                         sphere_resolution=20)

# ╔═╡ 93f37121-c331-4e84-80ac-66c4c6e2aef6
if baxe
	let
		blobsXe = jn.Petit.blob_analysis(trksxe1e, rblob)
		jn.Petit.write_blobs_csv(blobsXe, fnXe)
	end
end


# ╔═╡ d995912c-5e3d-46da-a19b-730498af7ff5
begin
  blobsXeFid = blobsXe[(blobsXe.confidence .> confCut) .& (blobsXe.trackLength .> trklCut) .& (blobsXe.eB2 .< eblb2Max), :]
	
	ntrXe.nfid = size(blobsXeFid)[1]
	effXe.efffid = ntrXe.nfid/ntrXe.ntot
end

# ╔═╡ f35f3201-6180-418d-aa4a-b495f4a9de09
begin
	histo_b1_b2(blobsXeFid.eB1, blobsXeFid.eB2)
end


# ╔═╡ 87c72baa-ba4c-4164-ba72-e68b6a7bb2d6
function track_energies_keV(tracks)
	E = Float64[]
	for i in 1:length(tracks)
		energy_kev = 1e+3 * sum(tracks[i].voxels.energy)
		push!(E, energy_kev)
	end
	E
end

# ╔═╡ 8f7a358e-b98a-41cc-9d30-2c7cce24aebb
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

# ╔═╡ 131fb43e-982e-43aa-83e1-028ae0e793a6
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

# ╔═╡ e086d2fc-e305-4b5b-a76b-f0f20808035e
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

# ╔═╡ 8e1206b0-024c-4866-a76c-1f963e2cb303
begin
	energy_trk_length(eBi, blobsBiFid.trackLength)
end

# ╔═╡ 63496eae-e07f-4259-86a7-459de5025613
begin
	energy_trk_length(eTl, blobsTlFid.trackLength)
end

# ╔═╡ e571208c-50dd-487d-acb5-e17d5785f92b
begin
	energy_trk_length(eBB, blobsBBFid.trackLength)
end


# ╔═╡ 64fa7f0f-5515-4990-839b-e92d5edf4f44
begin
	energy_trk_length(eXe, blobsXeFid.trackLength)
end

# ╔═╡ 198102c2-a488-4b26-a222-1cd2a6294528


# ╔═╡ 03dd8040-7c41-465b-b55c-00c7be539cb7


# ╔═╡ 65857e3a-32ac-4496-98a1-9d39a9694c37
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
	etrk = jn.Petit.in_range(etrkt, min_cut, max_cut)
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

# ╔═╡ 4657a7ef-3241-4ea9-81cd-4fb554a37c5b
begin
	effBlbBi = jn.Petit.fom_blobs(blobsBiFid.eB1, 
								  blobsBiFid.eB2,
								  min_cut=200.0,
						          max_cut=600.0,step=25.0)
	
	effEneBi, ffBi = fom_ecut(eBi; min_cut=2420.0, max_cut=2500.0,step=5.0)
	eb1_vs_eb2(blobsBiFid.eB1, blobsBiFid.eB2, effBlbBi,eblob_cut=425.0)
end


# ╔═╡ 960a2984-043e-4fab-a981-fe0f8e165e02
effBlbBi

# ╔═╡ c546afb2-f5cd-4290-ae6d-2805daa2a3fa
begin
  zng1, zng2, zp_dists, zeff1, zeff2, zsteps, zp_eff = gaussian_efficiency_analysis(2458.0,
							 2448.0,12.5, 12.5, nsteps=20, ff=ffBi)
  zp_dists  # Show distributions plot
 end

# ╔═╡ 4c76b220-b46a-4ba2-88f0-3bb787b802c2
zp_eff

# ╔═╡ 6a750201-47f6-4fcc-a551-7884c0280bc5
begin
	zeffbb = zeff1[1:end-1]
	zeffbi = zeff2[1:end-1]
	zenergies = zsteps[1:end-1]
	
	zd = Poisson.(ffBi*zeffbi)
	zfom  = zeffbb./std.(zd) 
	plot(zenergies, zfom,
                marker=:circle, markersize=5, lw=2,
                label="FOM ",
                xlabel="Step", ylabel="fom",
                title="FOM Curve",
                legend=:topright, ylim=(0, 1.5))
end

# ╔═╡ 56917df0-2e21-44c4-a7b7-3619a8661ab5
let
	ice = argmax(zfom)
	ecutx = zenergies[ice]
	efbb = zeffbb[ice]
	efbi =zeffbi[ice]
	md"""
	- Gaussian FOM maximizes at $(ecutx)
	- eff bb0nu = $(efbb)
	- eff bi214 = $(efbi)
	"""
end

# ╔═╡ 8c10a242-7411-4973-9fdc-a5e8ce461ee4
begin
	effBlbTl = jn.Petit.fom_blobs(blobsTlFid.eB1, 
										   blobsTlFid.eB2,
										   min_cut=200.0,
						                   max_cut=600.0,step=25.0)
	effEneTl, ffTl = fom_ecut(eTl; min_cut=2420.0, max_cut=2500.0,step=5.0)
	eb1_vs_eb2(blobsTlFid.eB1, blobsTlFid.eB2, effBlbTl,eblob_cut=425.0)
end

# ╔═╡ 178e3c90-ed84-43e9-b7d8-2282c3758961
begin
	pFomTlEne = fom_plot(effEneTl, " Efficiency Tl-208: Ecut")
end

# ╔═╡ e7ab75e6-33b7-4429-b9a2-b8600856a6fc
begin
	effBlbBB = jn.Petit.fom_blobs(blobsBBFid.eB1, 
								  blobsBBFid.eB2,
								  min_cut=200.0,
						          max_cut=600.0,step=25.0)
	
	effEneBB, ffBB = fom_ecut(eBB; min_cut=2420.0, max_cut=2500.0,step=5.0)
	eb1_vs_eb2(blobsBBFid.eB1, blobsBBFid.eB2, effBlbBB,eblob_cut=425.0)
end


# ╔═╡ 04bc7077-564a-470f-8138-c0ecf85d689f
begin
	pFomBBEne = fom_plot(effEneBB, " Efficiency bb0nu: Ecut")
end

# ╔═╡ 3a233b40-4e22-4a1b-9be5-c0990fc924bd
begin
	fomBiBBx = effEneBB.eff ./(sqrt.(effEneBi.eff))
	fomBiBB = replace(fomBiBBx, NaN => 0.0)
	plot(effEneBB.eblob2, fomBiBB,
          line = :solid,
          linewidth = 2,
          label = nothing,
          color = :blue,
          alpha = 0.5)
end

# ╔═╡ dafbd900-6a41-49bc-8afe-e9e31022ab49
let
	ice = argmax(fomBiBB)
	ecutx = effEneBB.eblob2[ice]
	efbb = effEneBB.eff[ice]
	efbi =effEneBi.eff[ice]
	eftl =effEneTl.eff[ice]
	effBi.effroi = efbi
	effTl.effroi = eftl
	effBB.effroi = efbb
	md"""
	- FOM maximizes at $(ecutx)
	- eff bb0nu = $(efbb)
	- eff bi214 = $(efbi)
	- eff tl208 = $(eftl)
	"""
end

# ╔═╡ 54497681-00c7-4940-9047-23cf1d86b236
begin
	fomBlbBiBBx = effBlbBB.eff ./(sqrt.(effBlbBi.eff))
	fomBlBBiBB = replace(fomBlbBiBBx, NaN => 0.0)
	plot(effBlbBB.eblob2, fomBlBBiBB,
          line = :solid,
          linewidth = 2,
          label = nothing,
          color = :blue,
          alpha = 0.5)
end

# ╔═╡ 50076b02-bd2c-4478-961b-9abee43d82ca
let
	ice = argmax(fomBlBBiBB)
	ecutx = effBlbBB.eblob2[ice]
	efbb = effBlbBB.eff[ice]
	efbi =effBlbBi.eff[ice]
	eftl =effBlbTl.eff[ice]
	effBi.effblb = efbi
	effTl.effblb = eftl
	effBB.effblb = efbb
	md"""
	### Blobs cut
	- FOM maximizes at $(ecutx)
	- eff bb0nu = $(efbb)
	- eff bi214 = $(efbi)
	- eff tl208 = $(eftl)
	"""
end

# ╔═╡ 9240273d-a92e-403e-a582-a7cb4789be60
begin
	effBlbXe = jn.Petit.fom_blobs(blobsXeFid.eB1, 
								  blobsXeFid.eB2,
								  min_cut=200.0,
						          max_cut=600.0,step=25.0)
	
	effEneXe, ffXe = fom_ecut(eXe; min_cut=2420.0, max_cut=2500.0,step=5.0)
	eb1_vs_eb2(blobsXeFid.eB1, blobsXeFid.eB2, effBlbXe,eblob_cut=425.0)
end

# ╔═╡ 48577121-4f6b-494f-8ae4-d8db338089fa
begin
	fomBlbXeBBx = effBlbBB.eff ./(sqrt.(effBlbXe.eff))
	fomBlBXeBB = replace(fomBlbXeBBx, NaN => 0.0)
	plot(effBlbBB.eblob2, fomBlBXeBB,
          line = :solid,
          linewidth = 2,
          label = nothing,
          color = :blue,
          alpha = 0.5)
end

# ╔═╡ 2837e158-e36c-49cd-9eaf-7ffd09f064b8
let
	ice = argmax(fomBlBXeBB)
	ecutx = effBlbBB.eblob2[ice]
	efbb = effBlbBB.eff[ice]
	efbi =effBlbBi.eff[ice]
	eftl =effBlbTl.eff[ice]
	efxe =effBlbXe.eff[ice]
	effBi.effblb = efbi
	effTl.effblb = eftl
	effBB.effblb = efbb
	effXe.effblb = efxe
	md"""
	### Blobs cut
	- FOM maximizes at $(ecutx)
	- eff bb0nu = $(efbb)
	- eff bi214 = $(efbi)
	- eff tl208 = $(eftl)
	- eff xe-137 = $(eftl)
	"""
end

# ╔═╡ ffc1b36a-7df4-48f6-8fff-5095d50d6039
let
	results = [effBlbBi, effBlbTl, effBlbXe]
	labels = ["Bi-214", "Tl-208", "1e"]
	fom_plots(results, labels)
end

# ╔═╡ e7060fc7-f3e6-477b-9963-072de8730de8


# ╔═╡ 85d7771b-fc92-4844-9086-0d8c3b41e210
function signal_eff(ehrx, rlow, rup; step=10.0) 
	countx = []
	norm = sum(ehrx.weights)
	for rx in rlow:step:rup
		push!(countx, jn.Petit.counts_in_range(ehrx, rx, rup))
	end
	countx = countx /norm
end

# ╔═╡ e18e2427-5818-4d3a-9dd8-e715aad2958f
begin
	function histogram_energy_primary(partdf)
	energies = 1e+3*jn.Petit.energy_primary(partdf)
	he, pe = jn.Petit.step_hist(energies;
                   nbins = 50,
                   xlim = (2300.0, 2700.0),
                   xlabel = "energy (keV)",
                   ylabel = "Frequency",
                   title ="Energy Primary")
	plot(pe)
end
function histogram_energy_deposited(hitsdf; emin=2300.0, emax=2700.0)
	energies = 1e+3*jn.Petit.energy_deposited(hitsdf)
	_, pe = jn.Petit.step_hist(energies;
                   nbins = 50,
				   logy=true, 
                   xlim = (emin, emax),
                   xlabel = "energy (keV)",
                   ylabel = "Frequency",
                   title ="Energy Deposited")
	plot(pe)
end
end

# ╔═╡ beb867b5-711f-47ba-b6d1-c3a6e63dd6f1
function histogram_alphas(df)
	
	h_x, p_x = jn.Petit.step_hist(df.x;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "X (mm)",
         ylabel = "Frequency",
         title=" X")

	h_y, p_y = jn.Petit.step_hist(df.y;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "Y (mm)",
         ylabel = "Frequency",
         title=" Y")

	h_z, p_z = jn.Petit.step_hist(df.z;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "Z (mm)",
         ylabel = "Frequency",
         title=" Z ")

	h_e, p_e = jn.Petit.step_hist(1e+3*df.energy;
         nbins = 50,
         xlim   = (0.0, 8000.0),
         xlabel = "energy (V)",
         ylabel = "Frequency",
         title=" E (keV) ")

	plot(p_x, p_y, p_z, p_e)
end

# ╔═╡ 33efdcc5-cc1d-48f4-81af-b5f65bb348f7
function histogram_stats(hitsdf)
	h_id, p_id = jn.Petit.step_hist(hitsdf.event_id;
         nbins = 50,
         xlim   = (0.0, 10000.0),
         xlabel = "event number",
         ylabel = "Frequency",
         title=" Events")
	
	nhits = jn.Petit.hits_per_all_events(hitsdf)
	h_nhits, p_nhits = jn.Petit.step_hist(nhits;
         nbins = 50,
         xlim   = (0.0, 500.0),
         xlabel = "number of hits per track",
         ylabel = "Frequency",
         title=" Hits per track")

	h_x, p_x = jn.Petit.step_hist(hitsdf.x;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "X (mm)",
         ylabel = "Frequency",
         title=" X")

	h_y, p_y = jn.Petit.step_hist(hitsdf.y;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "Y (mm)",
         ylabel = "Frequency",
         title=" Y")

	h_z, p_z = jn.Petit.step_hist(hitsdf.z;
         nbins = 50,
         xlim   = (-2000.0, 2000.0),
         xlabel = "Z (mm)",
         ylabel = "Frequency",
         title=" Z ")

	plot(p_nhits, p_x, p_y, p_z)
end

# ╔═╡ 6fc13d23-cb44-40a6-a52d-32b81855d6a0
function histogram_voxel_energy(df; emx=200.0, title="Voxel energy")
	h_vxe, p_vxe = jn.Petit.step_hist(1e+3*df.energy;
         nbins = 50,
         xlim   = (0.0, emx),
         xlabel = "Voxel energy (keV)",
         ylabel = "Frequency",
         title=title)
	plot(p_vxe)
end

# ╔═╡ 2d4732eb-47b0-498b-b0ff-c9dcf3951ecd
function histogram_distances(df; dmx=100.0, dcmx=20.0)
	vd = jn.Petit.voxel_distances(df)
	vcd = jn.Petit.voxel_closest_distance(df)
	_, p_vd = jn.Petit.step_hist(vd;
         nbins = 50,
         xlim   = (0.0, dmx),
         xlabel = "voxel distance (mm)",
         ylabel = "Frequency",
         title=" Voxel Distance ")
	_, p_vcd = jn.Petit.step_hist(vcd;
         nbins = 50,
         xlim   = (0.0, dcmx),
         xlabel = "voxel closest distance (mm)",
         ylabel = "Frequency",
         title=" Voxel Closest Distance ")
	
	plot(p_vd, p_vcd)
end

# ╔═╡ 769f772d-2695-468f-93a0-7655b5d08f10
function histogram_energies_trks(results)
	
	_, p_t1 = jn.Petit.step_hist(results.single_track.energy;
         nbins = 50,
         xlim   = (0.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E one track ")
	_, p_t2 = jn.Petit.step_hist(results.two_track_primary.energy;
         nbins = 50,
         xlim   = (0.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E two tracks ")
	_, p_t2s = jn.Petit.step_hist(results.two_track_secondary.energy;
         nbins = 50,
         xlim   = (0.0, 250.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E secondaries 2t ")
	_, p_t3s = jn.Petit.step_hist(results.three_track_secondary.energy;
         nbins = 50,
         xlim   = (0.0, 250.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" E secondaries 3 track ")
	
	plot(p_t1, p_t2, p_t2s, p_t3s)
end

# ╔═╡ 4bb62775-b172-4537-9087-c37cf9ea97f4
function plot_true_hits2(ghdf::GroupedDataFrame, index::Int; nbins=100)
    df = get_event(ghdf, index)
	plot_hits_evt(df; nbins)
end

# ╔═╡ 46d3b837-c9df-4e04-98ae-611054c9ad6d
function get_group_contents(filename::String; path::String = "/")
    h5open(filename, "r") do file
        if !haskey(file, path)
            error("Path $path not found in HDF5 file.")
        end
        group = file[path]
        return Dict(name => typeof(group[name]) for name in keys(group))
    end
end

# ╔═╡ a4d68f10-b2f3-4b48-b3d0-2baac6f93ca1
function get_subgroups(filename::String, path::String = "/")
    h5open(filename, "r") do file
        if !haskey(file, path)
            error("Path $path not found in HDF5 file.")
        end
        group = file[path]
        return [name for name in keys(group) if group[name] isa HDF5.Group]
    end
end

# ╔═╡ bbc7767d-2ab4-40ff-ad2a-e1f81b6ad368
function inspect_mc(filename::String)
    h5open(filename, "r") do fid
        dset = fid["MC"]
        println("Type: ", typeof(dset))
        data = read(dset)
        println("Read type: ", typeof(data))
        println("First element: ", data[1])
        return data
    end
end

# ╔═╡ 59c673be-1ead-488c-84bd-21bf14980313


# ╔═╡ 70914efe-c290-4f43-8c2a-f521281c977e


# ╔═╡ Cell order:
# ╠═947c237c-9852-40e9-a83f-c23666db90aa
# ╠═04b446d6-f34f-11ed-2565-0b15d65b6781
# ╠═871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
# ╠═349825ff-7ffe-4fa1-ba26-a772041f0323
# ╠═7504d7aa-a780-4956-99a5-08a7f9a462b2
# ╠═c9fc0547-0e73-4629-9909-e59c3d75169d
# ╠═8290afe2-e5cb-4793-b108-a313f2677119
# ╠═982a3890-7318-4034-b7fd-decafc288362
# ╠═b581c086-865c-4761-864d-ca8e549f8a75
# ╠═0d6dfea6-48ea-4645-b00c-568ea4e204af
# ╠═46e73d58-dc53-42bb-add7-088ac144ca24
# ╠═e086d2fc-e305-4b5b-a76b-f0f20808035e
# ╠═d2b5f0e9-4a79-4c33-afb2-7e970fb5565d
# ╠═65ec723a-6a38-4a2c-88e5-d6e1c59375bf
# ╠═1052ab37-288a-4c84-9bbe-3407a6a0dfaf
# ╠═641e9a21-5999-4eaa-8904-8e03fa952b49
# ╠═d61ec9b7-4922-45de-9b6f-e0a9f4386740
# ╠═f4e98f28-5f3a-49f5-a053-8ed532e92a60
# ╠═8fc7e228-fb0d-4c6c-ace9-8f6cb057a989
# ╠═22c204c8-9e36-41e0-a0d4-0b86c689739c
# ╠═de068020-a64c-43b4-9e78-8ade8a4f6274
# ╠═b5f850b0-19ec-4e96-96cc-4c43a91e05f3
# ╠═9bd589eb-0bbd-4964-b4ba-8a6f8cb7b38b
# ╠═a60ea1b2-1c4c-4565-ac34-d2580dd016e3
# ╠═ce0c0477-0b99-400a-a85e-f8e6cad5e094
# ╠═1c202304-d6ca-4dde-a7ac-f5cb2b961b93
# ╠═99258d62-017c-45e2-b9e5-63da00815166
# ╠═4f254ee3-c000-4c49-9f6b-0d7ef947c67d
# ╠═403b3c7f-d7c3-4762-8605-5d24fb087415
# ╠═bc7634bf-2d76-4966-95c1-70d3ea97e0ff
# ╠═067d7785-ce69-4cb6-b783-9fbeab4b9f32
# ╠═c5374971-03ed-4bec-88b6-e636ac400348
# ╠═898274d5-9add-4b0b-bf04-724bc6e06074
# ╠═d4bba489-68d1-4a21-994b-aebbb50fe601
# ╠═8e1206b0-024c-4866-a76c-1f963e2cb303
# ╠═73cb4e71-d232-4c97-ba8f-27d170823073
# ╠═4657a7ef-3241-4ea9-81cd-4fb554a37c5b
# ╠═960a2984-043e-4fab-a981-fe0f8e165e02
# ╠═046ec784-bc63-4dde-920d-36fc85f1ee8e
# ╠═48b32bbf-abbb-480e-9ba0-0df5f0777336
# ╠═10ed8e9a-e5f4-4652-9f1d-8e59eb605c95
# ╠═ba010d8d-f836-4517-b9e3-3f154ce88497
# ╠═e32cb510-5628-440b-a3ed-a02c0f27005e
# ╠═9cd7ccf0-a4e7-4445-a6b1-78ffa573eb3b
# ╠═4bcf1a8b-7101-41e9-bb30-9874275c20a6
# ╠═37a7d81a-ca8e-4ac5-b2ca-04b261beaddb
# ╠═01297be8-6234-4bae-91d5-963891b755ef
# ╠═9d7572dd-1ebe-426d-a8f9-3dd540775733
# ╠═aa3c48c6-3c83-47c3-ac59-ce12de38c371
# ╠═63496eae-e07f-4259-86a7-459de5025613
# ╠═8c10a242-7411-4973-9fdc-a5e8ce461ee4
# ╠═178e3c90-ed84-43e9-b7d8-2282c3758961
# ╠═2e1c70e9-6ecf-4e97-8705-6342028285aa
# ╠═a6e42cf4-1389-497f-90da-1df8297cec74
# ╠═759c3202-b6b7-45e7-ab11-75ed778273f7
# ╠═f398f371-5271-4af7-965b-73a7027eff0b
# ╠═00910312-2f5c-49bf-abda-11a6acf86249
# ╠═f7b87c23-dee5-462d-b4b9-b964947b9ca7
# ╠═131e6669-5a26-4d9c-9ef4-911aa779a20e
# ╠═18d8eaa0-5467-4123-bfc8-f52f377e7c72
# ╠═ac45ca81-79bc-4631-bf37-c08bbbf6159e
# ╠═cf3042b0-d5f3-4716-ab25-6ca736b491db
# ╠═794ac8f8-3eee-4a35-96b0-b5dca3615a07
# ╠═e571208c-50dd-487d-acb5-e17d5785f92b
# ╠═6ea54a0d-8b35-47d9-ac29-43b8f1d103db
# ╠═e7ab75e6-33b7-4429-b9a2-b8600856a6fc
# ╠═04bc7077-564a-470f-8138-c0ecf85d689f
# ╠═3a233b40-4e22-4a1b-9be5-c0990fc924bd
# ╠═dafbd900-6a41-49bc-8afe-e9e31022ab49
# ╠═54497681-00c7-4940-9047-23cf1d86b236
# ╠═50076b02-bd2c-4478-961b-9abee43d82ca
# ╠═c546afb2-f5cd-4290-ae6d-2805daa2a3fa
# ╠═4c76b220-b46a-4ba2-88f0-3bb787b802c2
# ╠═6a750201-47f6-4fcc-a551-7884c0280bc5
# ╠═56917df0-2e21-44c4-a7b7-3619a8661ab5
# ╠═ffb9d0d1-bca6-4b9a-87a7-e8c6280488fa
# ╠═4eeabfad-fb46-4d90-a586-aaf712c9eb47
# ╠═0397d23d-6c74-484b-9c8a-3989bab21239
# ╠═d7729d46-79ac-48c2-948a-41c91f78e2e3
# ╠═ac25b86e-97b3-4d54-b71a-a68a79a90995
# ╠═6ec7bab6-27c1-450d-a346-9abe523f1ca7
# ╠═75547508-ba93-4df5-af81-6026f88499a5
# ╠═2f87d310-8d44-4649-9dfb-440e160a24f4
# ╠═93f37121-c331-4e84-80ac-66c4c6e2aef6
# ╠═10a4a567-6337-418b-a70e-97647af0177b
# ╠═f333b7c7-bebb-47da-bbfc-b350f0f34a0d
# ╠═d995912c-5e3d-46da-a19b-730498af7ff5
# ╠═f35f3201-6180-418d-aa4a-b495f4a9de09
# ╠═64fa7f0f-5515-4990-839b-e92d5edf4f44
# ╠═9240273d-a92e-403e-a582-a7cb4789be60
# ╠═48577121-4f6b-494f-8ae4-d8db338089fa
# ╠═2837e158-e36c-49cd-9eaf-7ffd09f064b8
# ╠═d7648894-3e02-47ca-af42-a288190f699d
# ╠═f15ad4b6-873a-4ce3-8321-5a38e9704df2
# ╠═ffc1b36a-7df4-48f6-8fff-5095d50d6039
# ╠═d0f9e26b-2a3c-42a4-a826-7c8694f5d470
# ╠═8756b4ca-011b-4a6a-aff6-22488c8ac854
# ╠═f0ace932-f1d1-4e2b-9a6d-88b2cea4d849
# ╠═590f5634-843c-4589-ae4d-298293581732
# ╠═9bae1645-1485-411d-a6d3-61906c7c4194
# ╠═ffce2732-1271-4e7e-b3aa-4000f0355dff
# ╠═7dad709c-17e9-4738-b6c0-6d0c9bd70f09
# ╠═4f78d979-f707-44f9-892d-e278937c7c7e
# ╠═883c8a64-9cc2-4a89-95d5-4e906a903954
# ╠═75f229f8-1a11-4b9a-9ef6-331e98177caa
# ╠═32f432a9-68c1-4a65-9b98-8f0174bf08ae
# ╠═b59257e9-9a5e-467f-bbc3-f04cdf0e9a75
# ╠═a8d1a7d5-bfaa-4430-abde-01b167249f29
# ╠═98033f21-5720-4d22-94a8-a72e650a860b
# ╠═94e6ae79-2d97-461a-8a3a-bf100bb1771e
# ╠═2de993fe-47dc-459d-8000-708d85ed425f
# ╠═a1bac832-78d3-4ff5-8e05-bc8fc0eeef91
# ╠═0d4ea311-1602-43eb-8f57-f836434f4081
# ╠═87c72baa-ba4c-4164-ba72-e68b6a7bb2d6
# ╠═8f7a358e-b98a-41cc-9d30-2c7cce24aebb
# ╠═131fb43e-982e-43aa-83e1-028ae0e793a6
# ╠═198102c2-a488-4b26-a222-1cd2a6294528
# ╠═03dd8040-7c41-465b-b55c-00c7be539cb7
# ╠═65857e3a-32ac-4496-98a1-9d39a9694c37
# ╠═e7060fc7-f3e6-477b-9963-072de8730de8
# ╠═85d7771b-fc92-4844-9086-0d8c3b41e210
# ╠═e18e2427-5818-4d3a-9dd8-e715aad2958f
# ╠═beb867b5-711f-47ba-b6d1-c3a6e63dd6f1
# ╠═33efdcc5-cc1d-48f4-81af-b5f65bb348f7
# ╠═6fc13d23-cb44-40a6-a52d-32b81855d6a0
# ╠═2d4732eb-47b0-498b-b0ff-c9dcf3951ecd
# ╠═769f772d-2695-468f-93a0-7655b5d08f10
# ╠═4bb62775-b172-4537-9087-c37cf9ea97f4
# ╠═46d3b837-c9df-4e04-98ae-611054c9ad6d
# ╠═a4d68f10-b2f3-4b48-b3d0-2baac6f93ca1
# ╠═bbc7767d-2ab4-40ff-ad2a-e1f81b6ad368
# ╠═59c673be-1ead-488c-84bd-21bf14980313
# ╠═70914efe-c290-4f43-8c2a-f521281c977e
