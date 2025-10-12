### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

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

# ╔═╡ fda162d5-3552-4716-968f-d2c72a093c39
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

# ╔═╡ 8290afe2-e5cb-4793-b108-a313f2677119
md"""
## Bi-214
"""

# ╔═╡ 5bcffb43-a8f7-448f-86c2-f3261d489bbf
md""" 
- E, X, Y, Z
"""

# ╔═╡ 84e6f57f-50a7-43fc-82a4-4ef238449fc7
md"""
- voxel size for 10 % He
"""

# ╔═╡ 982a3890-7318-4034-b7fd-decafc288362
function voxel_size(dx, dy, p, l)
	sqrt(l)*(dx + dy) /(2.0*sqrt(p))
end

# ╔═╡ b7e69022-a4d9-4d4b-841f-9257bffe0b07
vhe = voxel_size(0.75, 1.6, 15.0, 100.0)

# ╔═╡ a60ea1b2-1c4c-4565-ac34-d2580dd016e3
md"""
### Build Tracks: example
"""

# ╔═╡ 03cd5a27-1ace-43b3-8a31-4e34629b082d
md"""
### Find blobs: example
"""

# ╔═╡ 99258d62-017c-45e2-b9e5-63da00815166
md"""
#### Compute E, TRKL, EBLOB1, EBLOB2 for all tracks
"""

# ╔═╡ 2e1c70e9-6ecf-4e97-8705-6342028285aa
md"""
## bb0nu
"""

# ╔═╡ ffb9d0d1-bca6-4b9a-87a7-e8c6280488fa
md"""
## Signal and background efficiency
"""

# ╔═╡ d3fd5874-8de1-451c-89b1-3dc05a65cef4
md"""
## Bi-214 and bbonu efficiencies
"""

# ╔═╡ e0a4575c-afd3-4619-8fad-e321d573d1ec
md"""
## Tl208
"""

# ╔═╡ bd5f963d-e833-41a0-a46f-f014f47dd1ca
md"""
#### Energy cut
"""

# ╔═╡ a77652cb-6ab9-4f85-b8a2-5a21607444cd
md"""
#### Blob cuts
"""

# ╔═╡ 0397d23d-6c74-484b-9c8a-3989bab21239
md"""
## Xe-137
"""

# ╔═╡ ad4d3db9-43f6-48ee-84ba-a5a360cb258b
md"""
#### Blob cuts
"""

# ╔═╡ d0f9e26b-2a3c-42a4-a826-7c8694f5d470
md"""
# Functions
"""

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

# ╔═╡ b5f850b0-19ec-4e96-96cc-4c43a91e05f3
begin
	trkb214, metabi214,  nbi214tot, nbi2141t  =  get_bi214_tracks()
	
	md"""
	- Total number of Bi-214 events generated = $(nbi214tot)
	- Total number of Bi-214 events 1 trk  = $(nbi2141t)
	"""
end

# ╔═╡ ce0c0477-0b99-400a-a85e-f8e6cad5e094
begin
	xresult, xstart_voxel, xend_voxel, xtrack_length, energy_kev=find_track_extremes(trkb214, i=i)
	md"""
	#### Find track Extremes
	- confidence = $(xresult.confidence)
	- start voxel: x = $(xstart_voxel.x), y = $(xstart_voxel.y), z = $(xstart_voxel.z)
	- end voxel: x = $(xend_voxel.x), y = $(xend_voxel.y), z = $(xend_voxel.z)
	- track length L =$(xtrack_length)
	- energy = $(energy_kev) KeV
	"""
end

# ╔═╡ bf77b535-485a-4608-aa4b-1c15950f3f9b
begin
	
	blobs = jn.Petit.energy_in_spheres_around_extremes(trkb214[i], xresult, rblob)
	eb1 = blobs.blob1_energy * 1e+3
	eb2 = blobs.blob2_energy * 1e+3
	nb1 = blobs.blob1_voxel_count
	nb2 = blobs.blob2_voxel_count
	md"""
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

# ╔═╡ 4853838b-7292-49ad-b3d8-f116a69da7d6
begin
	trktl1, metatl1, ntltot1, ntl1t1 =  get_tl208_tracks(bkgnd="copperbkg", tag="*st3mm.b*")

	trktl2, metatl2, ntltot2, ntl1t2 =  get_tl208_tracks(bkgnd="innerbkg", tag="*st3mm.b*")

	trktl = vcat(trktl1, trktl2)
	ntltot = ntltot1 + ntltot2 
	ntl1t = ntl1t1 + ntl1t2
	
	md"""
	- Total number of tl208 events generated in copper = $@sprintf("%.3e", ntl1t1) 
	- Total number of tl208 events 1 trk from copper  = $(ntltot1)
	- Total number of tl208 events generated in inner = $(ntl1t2)
	- Total number of tl208 events (inner) 1 trk  = $(ntltot2)
	- Total number of tl208 1 trk in detector = $(length(trktl))
	"""
end

# ╔═╡ 25650a53-6284-4c57-a590-71283d15a3fd
begin
	xresulttl, xstart_voxeltl, xend_voxeltl, xtrack_lengthtl, energy_kevtl=find_track_extremes(trktl, i=i)
	md"""
	#### Find track Extremes
	- confidence = $(xresulttl.confidence)
	- start voxel: x = $(xstart_voxeltl.x), y = $(xstart_voxeltl.y), z = $(xstart_voxeltl.z)
	- end voxel: x = $(xend_voxeltl.x), y = $(xend_voxeltl.y), z = $(xend_voxeltl.z)
	- track length L =$(xtrack_lengthtl)
	- energy = $(energy_kevtl) keV
	"""
end

# ╔═╡ ede62a96-39b6-4d2b-ab8a-528048c56556
begin
	
	blobstl = jn.Petit.energy_in_spheres_around_extremes(trktl[i], xresulttl, rblob)
	etlb1 = blobstl.blob1_energy * 1e+3
	etlb2 = blobstl.blob2_energy * 1e+3
	ntlb1 = blobstl.blob1_voxel_count
	ntlb2 = blobstl.blob2_voxel_count
	md"""
	- blob 1 energy = $(round(etlb1, digits=1)) keV
	- blob 2 energy = $(round(etlb2, digits=1)) keV
	- blob 1 # of voxels = $(ntlb1)
	- blob 2 # of voxels = $(ntlb2)
	"""
end

# ╔═╡ e7e8e6a4-f101-4f55-b8bf-fb9ffd5fa023
jn.Petit.plot_track_blobs(trktl[i], rblob;
                         markersize_voxels=3.0,
                         show_connections=true,
                         alpha_connections=0.2,
                         alpha_spheres=0.3,
                         sphere_resolution=20)

# ╔═╡ 883c8a64-9cc2-4a89-95d5-4e906a903954
function get_xe137_tracks()
	files = Glob.glob("*st3mm.b*", "/Users/jjgomezcadenas/Data/HD5t/precdr/xe137/")
	tracks, metadatas = jn.Petit.chain_track_files(files)
	
	ntl208= metadatas[1]["nevents_from_config"]
	## NB ntl208 takes only one file because we are splitting a very
	## large file in small segments, unlike the case for Bi214
	return tracks, metadatas, ntl208, length(tracks)
end

# ╔═╡ d7729d46-79ac-48c2-948a-41c91f78e2e3
begin
	trk1e, meta1e, n1etot, n1e1t =  get_xe137_tracks()
	md"""
	- Total number of tl208 events generated = $(n1e1t)
	- Total number of tl208 events 1 trk  = $(n1etot)
	"""
end

# ╔═╡ 4478f154-a096-4bdb-a758-c1e300a1de45
begin
	xresult1e, xstart_voxel1e, xend_voxel1e, xtrack_length1e, energy_kev1e=find_track_extremes(trk1e, i=i)
	md"""
	#### Find track Extremes
	- confidence = $(xresult1e.confidence)
	- start voxel: x = $(xstart_voxel1e.x), y = $(xstart_voxel1e.y), z = $(xstart_voxel1e.z)
	- end voxel: x = $(xend_voxel1e.x), y = $(xend_voxel1e.y), z = $(xend_voxel1e.z)
	- track length L =$(xtrack_length1e)
	- energy = $(energy_kev1e) keV
	"""
end

# ╔═╡ a7d23b8d-0854-4751-a216-d519e85e2d8b
begin
	
	blobs1e = jn.Petit.energy_in_spheres_around_extremes(trk1e[i], xresult1e, rblob)
	e1eb1 = blobs1e.blob1_energy * 1e+3
	e1eb2 = blobs1e.blob2_energy * 1e+3
	n1eb1 = blobs1e.blob1_voxel_count
	n1eb2 = blobs1e.blob2_voxel_count
	md"""
	- blob 1 energy = $(round(e1eb1, digits=1)) keV
	- blob 2 energy = $(round(e1eb2, digits=1)) keV
	- blob 1 # of voxels = $(n1eb1)
	- blob 2 # of voxels = $(n1eb2)
	"""
end

# ╔═╡ 63fcc371-d276-4e15-b812-da1974b22df3
jn.Petit.plot_track_blobs(trk1e[i], rblob;
                         markersize_voxels=3.0,
                         show_connections=true,
                         alpha_connections=0.2,
                         alpha_spheres=0.3,
                         sphere_resolution=20)

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

# ╔═╡ a6e42cf4-1389-497f-90da-1df8297cec74
begin
	trkbb, metabb, nbbtot, nbbb1t =  get_tracks("/Users/jjgomezcadenas/Data/HD5t/precdr/bb0nu")

	nbbtot = 1000 # only read 1000 events from bb0nu file
	md"""
	- Total number of bb0nu events generated = $(nbbb1t)
	- Total number of bb0nu events 1 trk  = $(nbbtot)
	"""
end

# ╔═╡ c8bf9933-3676-4b6f-b1c4-9fe73036b8f2
begin
	
	xresultb0, xstart_voxelb0, xend_voxelb0, xtrack_lengthb0, energy_kevb0=find_track_extremes(trkbb, i=i)
	md"""
	#### Find track Extremes
	- confidence = $(xresultb0.confidence)
	- start voxel: x = $(xstart_voxelb0.x), y = $(xstart_voxelb0.y), z = $(xstart_voxelb0.z)
	- end voxel: x = $(xend_voxelb0.x), y = $(xend_voxelb0.y), z = $(xend_voxelb0.z)
	- track length L =$(xtrack_lengthb0)
	- energy = $(energy_kevb0) keV
	"""
end

# ╔═╡ db13c2ed-535b-4beb-9c71-bdb7f9fbd2c1
begin
	
	blobsb0 = jn.Petit.energy_in_spheres_around_extremes(trkbb[i], xresultb0, rblob)
	ebb1 = blobsb0.blob1_energy * 1e+3
	ebb2 = blobsb0.blob2_energy * 1e+3
	nbb1 = blobsb0.blob1_voxel_count
	nbb2 = blobsb0.blob2_voxel_count
	md"""
	- blob 1 energy = $(round(ebb1, digits=1)) keV
	- blob 2 energy = $(round(ebb2, digits=1)) keV
	- blob 1 # of voxels = $(nbb1)
	- blob 2 # of voxels = $(nbb2)
	"""
end

# ╔═╡ b6353c16-9a6c-474c-84c2-dc8baf0c2dfa
jn.Petit.plot_track_blobs(trkbb[i], rblob;
                         markersize_voxels=3.0,
                         show_connections=true,
                         alpha_connections=0.2,
                         alpha_spheres=0.3,
                         sphere_resolution=20)

# ╔═╡ 32f432a9-68c1-4a65-9b98-8f0174bf08ae
function histo_b1_b2(xeb1, xeb2, xcon; xcut=0.7)
	yeb1 = xeb1[xcon .> xcut]
	yeb2 = xeb2[xcon .> xcut]
	h_b1, p_b1 = jn.Petit.step_hist(yeb1;
	         nbins = 20,
	        #xlim   = (0.0, 1.0),
	         xlabel = "Eb1",
	        ylabel = "Frequency",
	         title=" Eb1 ")
	h_b2, p_b2 = jn.Petit.step_hist(yeb2;
	         nbins = 20,
	        #xlim   = (0.0, 1.0),
	         xlabel = "Eb2",
	        ylabel = "Frequency",
	         title=" Eb2 ")
	
	plot(p_b1, p_b2)
end

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

# ╔═╡ 24518656-35e4-48ff-9b1e-7da70c3a5994
begin
	Xbi, Ybi, Zbi =track_positions(trkb214)
	histo_xyz(Xbi,Ybi,Zbi,"Bi214")
end

# ╔═╡ 00c71b09-25bc-41a3-b621-929caaddd626
begin
	Xbb, Ybb, Zbb =track_positions(trkbb)
	histo_xyz(Xbb,Ybb,Zbb,"bb0nu")
end

# ╔═╡ 887e37cd-8e50-44c3-b239-1936a9b0baef
begin
	Xtl, Ytl, Ztl =track_positions(trktl)
	histo_xyz(Xtl,Ytl,Ztl,"Tl208")
end

# ╔═╡ b59703c8-2851-494b-b22d-d626ebf9a2e3
let
	Xtl, Ytl, Ztl =track_positions(trk1e)
	histo_xyz(Xtl,Ytl,Ztl,"1e")
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

# ╔═╡ 55b40e99-5ded-4757-ba20-28e80c144777
begin
	ebi = track_energies_keV(trkb214)
	histo_trak_energy(ebi, "Bi214-Energy")
end

# ╔═╡ b6b7d4ab-0d75-49f1-bfd2-766fe60bf585
begin
	ebb = track_energies_keV(trkbb)
	histo_trak_energy(ebb, "bb0nu")
end

# ╔═╡ 78d16892-dab9-4b2a-97b9-06e49592f0c6
begin
	etl = track_energies_keV(trktl)
	histo_trak_energy(etl, "Tl208-Energy")
end

# ╔═╡ ec09be7e-90d7-48ee-b7d2-9125c8ca7329
let
	etl = track_energies_keV(trk1e)
	histo_trak_energy(etl, "1e-Energy")
end

# ╔═╡ 98033f21-5720-4d22-94a8-a72e650a860b
function blob_analysis(strks, r)
	eB1=Float64[]
	eB2=Float64[]
	CON = Float64[]
	E = track_energies_keV(strks)
	TL = Float64[]
	for strk in strks
		xresult = jn.Petit.walk_track_from_extremes(strk)
	  	xtrack_length = xresult.total_length
		confidence = xresult.confidence
		
		blobs = jn.Petit.energy_in_spheres_around_extremes(strk, xresult, r)
		eb1 = blobs.blob1_energy * 1e+3
		eb2 = blobs.blob2_energy * 1e+3
		push!(eB1, eb1)
		push!(eB2, eb2)
		push!(CON, confidence)
		push!(TL, xtrack_length)
	end
	return CON, TL, E, eB1, eB2
		
	
end

# ╔═╡ 403b3c7f-d7c3-4762-8605-5d24fb087415
begin
	xcon, xtl, xe, xeb1, xeb2 = blob_analysis(trkb214, rblob)
	hek, p_ek, p_c, p_tl = histo_energy_trkl_conf(xcon, xtl, xe)
	
	xebis = jn.Petit.smear_histogram(hek, erex)
	ehbi, phbi = jn.Petit.step_hist(xebis;
         nbins = 40,
         xlim   = (2300.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" Bi 214 E ")
	plot(p_c, p_tl, p_ek, phbi)
end

# ╔═╡ d4bba489-68d1-4a21-994b-aebbb50fe601
histo_b1_b2(xeb1, xeb2, xcon)

# ╔═╡ 73cb4e71-d232-4c97-ba8f-27d170823073
begin
fgbi, pfgbi =jn.Petit.plot_fit_gauss(xebis, "Track Energy", "Counts",
                        30, 2400.0, 2500.0;
                        xgmin=2400.0, xgmax=2500.0, gbins=30)
	plot(pfgbi)
end

# ╔═╡ d1c7141a-1f18-4079-ada2-463593c5e20e
begin
	xconb0, xtlb0, xeb0, xebb1, xebb2 = blob_analysis(trkbb, rblob)
	hekb, p_ekb, p_cb, p_tlb = histo_energy_trkl_conf(xconb0, xtlb0, xeb0)
	xebbs = jn.Petit.smear_histogram(hekb, erex)
	ehbb, phbb = jn.Petit.step_hist(xebbs;
         nbins = 40,
         xlim   = (2300.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" Bi 214 E ")
	plot(p_cb, p_tlb, p_ekb, phbb)
end

# ╔═╡ f0f3e0a6-6a6a-4e95-b889-6a8995faf266
histo_b1_b2(xebb1, xebb2,xconb0)

# ╔═╡ f5900835-ec28-427b-837f-09620735976b
begin
fgbb, pfgbb =jn.Petit.plot_fit_gauss(xebbs, "Track Energy", "Counts",
                        30, 2410.0, 2500.0;
                        xgmin=2400.0, xgmax=2500.0, gbins=30)
plot(pfgbi, pfgbb,layout=(1,2), size=(1400, 800))
end

# ╔═╡ f1f036c6-d455-499e-b7c6-2a6716cd3485
begin
	xtlcon, xtltl, xtle, xtleb1, xtleb2 = blob_analysis(trktl, rblob)
	htlek, p_tlek, p_tlc, p_tltl = histo_energy_trkl_conf(xtlcon, xtltl, xtle)
	
	xetls = jn.Petit.smear_histogram(htlek, erex)
	ehtl, phtl = jn.Petit.step_hist(xetls;
         nbins = 40,
         xlim   = (2300.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" Tl-208 E ")
	plot(p_tlc, p_tltl, p_tlek, phtl)
end

# ╔═╡ 6de3499c-3cf8-4985-9a2b-e69f173044a8
length(xtlcon)

# ╔═╡ f32501b3-a66c-4c3d-aff0-ac28f05abc8d
histo_b1_b2(xtleb1, xtleb2, xtlcon)

# ╔═╡ dd12ed13-1cd4-48d8-bf28-25baa2a090cd
begin
	x1econ, x1etl, x1e, x1eeb1, x1eeb2 = blob_analysis(trk1e, rblob)
	h1eek, p_1eek, p_1ec, p_1etl = histo_energy_trkl_conf(x1econ, x1etl, x1e)
	
	xe1es = jn.Petit.smear_histogram(h1eek, erex)
	eh1e, ph1e = jn.Petit.step_hist(xe1es;
         nbins = 40,
         xlim   = (2300.0, 2700.0),
         xlabel = " E (keV)",
         ylabel = "Frequency",
         title=" TXe-137 E ")
	plot(p_1ec, p_1etl, p_1eek, ph1e)
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

# ╔═╡ a8d1a7d5-bfaa-4430-abde-01b167249f29
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
energy_trk_length(xebis, xtl)

# ╔═╡ 9353df52-9492-44b6-bb80-bdaec11ffb27
energy_trk_length(xebbs, xtlb0)

# ╔═╡ d208210b-ffd9-4bd1-9cd3-c828bcdef732
energy_trk_length(xetls, xtltl)

# ╔═╡ 198102c2-a488-4b26-a222-1cd2a6294528
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
        xlims = (0, maximum([maximum(eblob1), maximum(eblob2)]) * 1.1),
        ylims = (0, maximum([maximum(eblob1), maximum(eblob2)]) * 1.1),
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

# ╔═╡ 03dd8040-7c41-465b-b55c-00c7be539cb7
function fom_blobs(eblob1, eblob2; min_cut=100.0, max_cut=1000.0,step=50.0)

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
    # Total number of events
    n_total = length(eblob1)


    # Initialize results arrays
    cuts = Float64[]
    efficiencies = Float64[]
    errors = Float64[]

    # Iterate over cut values
    cut_value = min_cut
    while cut_value <= max_cut
        # Count events passing the cut
        n_pass = sum(eblob2 .>= cut_value)

        # Calculate efficiency and error
        
        eff, err = simple_binomial_error(n_pass, n_total)
        

        push!(cuts, cut_value)
        push!(efficiencies, eff)
        push!(errors, err)

        cut_value += step
    end

    # Create output dataframe
    results_df = DataFrame(
        eblob2 = cuts,
        eff = round.(efficiencies, digits=4),
        err = round.(errors, digits=4)
    )

    return results_df
end

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

# ╔═╡ e7060fc7-f3e6-477b-9963-072de8730de8
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


# ╔═╡ 2de993fe-47dc-459d-8000-708d85ed425f
function eb1_vs_eb2(xeb1, xeb2, effdf; eblob_cut=600.0)
	p1 = eb1_eb2_plot(xeb1, xeb2, title="Bi-214", eblob_cut=eblob_cut)
	p2 = fom_plot(effdf, "FOM: Bi-214")
    plot(p1, p2,
        layout = (1, 2),
        size = (900, 500),
        plot_title = "energy blob1 vs energy blob2",
        plot_titlefontsize = 16,
        margin = 5Plots.mm
    )
end

# ╔═╡ 4657a7ef-3241-4ea9-81cd-4fb554a37c5b
begin
	effbBi = fom_blobs(xeb1, xeb2; min_cut=200.0, max_cut=500.0,step=25.0)
	effeBi, ffBi = fom_ecut(xebis; min_cut=2420.0, max_cut=2500.0,step=5.0)
	eb1_vs_eb2(xeb1, xeb2, effbBi,eblob_cut=425.0)
end


# ╔═╡ 21d074a4-0ce2-481b-9169-f3361eac108d
pfombie = fom_plot(effeBi, " Efficiency Bi-214: Ecut")

# ╔═╡ 3a790c47-681a-4880-adbe-7172e4f8fd13
effeBi

# ╔═╡ 2643efff-0077-4e37-a283-146194a56a75
effeBi.eff

# ╔═╡ e24e87ee-f6dc-4dac-ab76-c24d2deb4e3f
 begin
  ng1, ng2, p_dists, eff1, eff2, steps, p_eff = gaussian_efficiency_analysis(fgbb, fgbi,xmin=2420.0, xmax=2500.0, nsteps=20, ff=ffBi)
  p_dists  # Show distributions plot
 end

# ╔═╡ 5c33c5db-0583-4c66-92bd-c3cb985d30f7
p_eff    # Show efficiency plot

# ╔═╡ 270136ac-f0d2-447b-b4a5-da44709ee79a
begin
	effbb = eff1[1:end-1]
	effbi = eff2[1:end-1]
	energies = steps[1:end-1]
	
	d = Poisson.(ffBi*effbi)
	fom  = effbb./std.(d) 
	plot(energies, fom,
                marker=:circle, markersize=5, lw=2,
                label="FOM ",
                xlabel="Step", ylabel="fom",
                title="FOM Curve",
                legend=:topright, ylim=(0, 1.5))
end

# ╔═╡ ea7e8b85-7a18-4f9f-b1e0-60643ea44a29
begin
	ie = argmax(fom)
	ecut = energies[ie]
	cefbb = effbb[ie]
	cefbi = effbi[ie]
	md"""
	- FOM maximizes at $(ecut)
	- eff bb0nu = $(cefbb)
	- eff bi214 = $(cefbi)
	"""
end

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

# ╔═╡ 988b5474-d301-4527-b553-cebe8decbec5
begin
	zie = argmax(zfom)
	zecut = energies[zie]
	zcefbb = effbb[zie]
	zcefbi = effbi[zie]
	md"""
	- FOM maximizes at $(zecut)
	- eff bb0nu = $(zcefbb)
	- eff bi214 = $(zcefbi)
	"""
end

# ╔═╡ 785df551-c770-421b-a2c1-91d064b47de5
begin
	effb0b = fom_blobs(xebb1, xebb2; min_cut=200.0, max_cut=500.0,step=25.0)
	effe0b, ffbb = fom_ecut(xebbs; min_cut=2420.0, max_cut=2500.0,step=5.0)
	eb1_vs_eb2(xebb1, xebb2, effb0b, eblob_cut=425.0)
end

# ╔═╡ 0d365351-e12c-43e2-87a7-6f46c391aaa2
effe0b

# ╔═╡ 73ff1c94-7317-4a56-a96c-3dd9c6cf69a8
effe0b.eff

# ╔═╡ 995dd661-1dac-4f62-b91c-18f2310c797d
let
	ie = 8
	ecut = effe0b.eblob2[ie]
	ceffbi = effeBi.eff[ie]
	cefbb = effe0b.eff[ie]
	md"""
	- FOM maximizes bkgnd suppresion at $(ecut)
	- eff bb0nu = $(cefbb)
	- eff bi214 = $(ceffbi)
	"""
end

# ╔═╡ 811eeb28-c626-414d-b6fb-41dce9de1d05
begin
	pfombbe = fom_plot(effe0b, "Efficiency: bb0nu: Ecut")
	plot(pfombie, pfombbe)
end

# ╔═╡ 35c2ad9f-a206-45e3-8258-5e04c645dd13
begin
	yeffbb = effb0b.eff
	yeffbi = effbBi.eff
	eblb = effb0b.eblob2
	
	yfom  = yeffbb./sqrt.(yeffbi) 
	plot(eblb, yfom,
                marker=:circle, markersize=5, lw=2,
                label="FOM ",
                xlabel="Step", ylabel="fom",
                title="FOM Curve blobs",
                legend=:topright, ylim=(0, 5.5))
end

# ╔═╡ 77e38166-f9a8-4a00-b034-ddc69a52a345
begin
	yie = argmax(yfom)
	yecut = eblb[yie]
	ycefbb = yeffbb[yie]
	ycefbi = yeffbi[yie]
	md"""
	- FOM maximizes at $(yecut)
	- eff bb0nu = $(ycefbb)
	- eff bi214 = $(ycefbi)
	"""
end

# ╔═╡ 73eed3c2-8dbe-45bf-85af-5cc48ed6667d
begin
	eff1tbi = nbi2141t/nbi214tot
	effroibi = eff1tbi * zcefbi
	effblbbi = effroibi * ycefbi
	md"""
	### Efficiency for Bi-214 (maximize Bi-214 suppression)
	- 1 trk no blobs = $(eff1tbi)
	- in ROI   = $(effroibi)
	- blob cut = $(effblbbi)
	"""
end

# ╔═╡ 0d3b2f02-e66b-469e-96e9-1bb702466f87
begin
	eff1tbb = nbbb1t/nbbtot 
	effroibb = eff1tbb * zcefbb
	effblbbb = effroibb * ycefbb
	md"""
	### Efficiency for bb0nu (maximize Bi-214 suppression)
	- 1 trk no blobs = $(eff1tbb)
	- in ROI   = $(effroibb)
	- blob cut = $(effblbbb)
	"""
end

# ╔═╡ 60e2d0eb-07dc-43ee-9c7b-67c46765a31a
begin
	imx2 = 5
	y2ecut = eblb[imx2]
	y2cefbb = yeffbb[imx2]
	y2cefbi = yeffbi[imx2]
	md"""
	- FOM alternative at $(y2ecut)
	- eff bb0nu = $(y2cefbb)
	- eff bi214 = $(y2cefbi)
	"""
end

# ╔═╡ 37ca7c05-f216-4302-99dd-c41b295e100c
let
	eff1tbi = nbi2141t/nbi214tot
	effroibi = eff1tbi * zcefbi
	effblbbi = effroibi * y2cefbi
	md"""
	### Efficiency for Bi-214 (maximize bbonu efficiency)
	- 1 trk no blobs = $(eff1tbi)
	- in ROI   = $(effroibi)
	- blob cut = $(effblbbi)
	"""
end

# ╔═╡ c958c249-5c8b-4e59-befe-e69b98d84cf5
let
	eff1tbb = nbbb1t/nbbtot 
	effroibb = eff1tbb * zcefbb
	effblbbb = effroibb * y2cefbb
	md"""
	### Efficiency for bb0nu (maximize bb0nu efficiency)
	- 1 trk no blobs = $(eff1tbb)
	- in ROI   = $(effroibb)
	- blob cut = $(effblbbb)
	"""
end

# ╔═╡ fd487373-3d9c-4f79-9a30-3765e7ee242f
eblb

# ╔═╡ 6a4318cf-6897-4ff0-9d82-cde4ab6c9e40
let
	ie = 8
	effetl, fftl = fom_ecut(xetls; min_cut=2420.0, max_cut=2500.0,step=5.0)
	ecut = effe0b.eblob2[ie]
	cefftl = effetl.eff[ie]
	cefbb = effe0b.eff[ie]
	md"""
	- FOM maximizes bkgnd suppresion at $(ecut)
	- eff bb0nu = $(cefbb)
	- eff bi214 = $(cefftl)
	"""
end

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

# ╔═╡ 899c2928-2457-4748-9fe6-3f60c85a8889
begin
	efom, pefom = eff_ecut(effe0b, effeBi, ffBi)
	plot(pefom)
end

# ╔═╡ 7987b9c8-275b-4c74-baeb-1a5c5ca42fde
begin
	effbtl = fom_blobs(xtleb1, xtleb2; min_cut=200.0, max_cut=500.0,step=25.0)
	eb1_vs_eb2(xtleb1, xtleb2, effbtl,eblob_cut=425.0)
end

# ╔═╡ 020c589a-ac92-4031-bf29-f482c525caf6
begin
	yefftl = effbtl.eff
	
	yfomtl  = yeffbb./sqrt.(yefftl) 
	plot(eblb, yfomtl,
                marker=:circle, markersize=5, lw=2,
                label="FOM ",
                xlabel="Step", ylabel="fom",
                title="FOM Tl-208 blobs",
                legend=:topright, ylim=(0, 5.5))
end

# ╔═╡ 6589dac6-db2d-4e2a-82c8-11524df319db
yefftl

# ╔═╡ 83cd563c-5743-4a89-8fea-1d512cd2939a
let
	imx2 = 10
	y2ecut = eblb[imx2]
	y2cefbb = yeffbb[imx2]
	y2cefbi = yefftl[imx2]
	md"""
	- FOM alternative at $(y2ecut)
	- eff bb0nu = $(y2cefbb)
	- eff bi214 = $(y2cefbi)
	"""
end

# ╔═╡ e5c0d31f-7a7f-4093-bcde-1e66fd50b5c7
begin
	effb1e = fom_blobs(x1eeb1, x1eeb2; min_cut=200.0, max_cut=500.0,step=25.0)
	eb1_vs_eb2(x1eeb1, x1eeb2, effb1e,eblob_cut=425.0)
end

# ╔═╡ dfe62938-b4e9-460a-accb-9f852509f038
begin
	yeff1e = effb1e.eff
	
	yfom1e  = yeffbb./sqrt.(yeff1e) 
	plot(eblb, yfom1e,
                marker=:circle, markersize=5, lw=2,
                label="FOM ",
                xlabel="Step", ylabel="fom",
                title="FOMXe-137 blobs",
                legend=:topright, ylim=(0, 5.5))
end

# ╔═╡ 98207c81-9f05-476f-846e-e819398ca81a
let
	imx2 = 10
	y2ecut = eblb[imx2]
	y2cefbb = yeffbb[imx2]
	y2cefbi = yeff1e[imx2]
	md"""
	- FOM alternative at $(y2ecut)
	- eff bb0nu = $(y2cefbb)
	- eff bi214 = $(y2cefbi)
	"""
end

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
# ╠═fda162d5-3552-4716-968f-d2c72a093c39
# ╠═8290afe2-e5cb-4793-b108-a313f2677119
# ╠═b5f850b0-19ec-4e96-96cc-4c43a91e05f3
# ╠═5bcffb43-a8f7-448f-86c2-f3261d489bbf
# ╠═55b40e99-5ded-4757-ba20-28e80c144777
# ╠═24518656-35e4-48ff-9b1e-7da70c3a5994
# ╠═84e6f57f-50a7-43fc-82a4-4ef238449fc7
# ╠═982a3890-7318-4034-b7fd-decafc288362
# ╠═b7e69022-a4d9-4d4b-841f-9257bffe0b07
# ╠═a60ea1b2-1c4c-4565-ac34-d2580dd016e3
# ╠═ce0c0477-0b99-400a-a85e-f8e6cad5e094
# ╠═03cd5a27-1ace-43b3-8a31-4e34629b082d
# ╠═bf77b535-485a-4608-aa4b-1c15950f3f9b
# ╠═1c202304-d6ca-4dde-a7ac-f5cb2b961b93
# ╠═99258d62-017c-45e2-b9e5-63da00815166
# ╠═403b3c7f-d7c3-4762-8605-5d24fb087415
# ╠═d4bba489-68d1-4a21-994b-aebbb50fe601
# ╠═8e1206b0-024c-4866-a76c-1f963e2cb303
# ╠═73cb4e71-d232-4c97-ba8f-27d170823073
# ╠═4657a7ef-3241-4ea9-81cd-4fb554a37c5b
# ╠═21d074a4-0ce2-481b-9169-f3361eac108d
# ╠═3a790c47-681a-4880-adbe-7172e4f8fd13
# ╠═2e1c70e9-6ecf-4e97-8705-6342028285aa
# ╠═a6e42cf4-1389-497f-90da-1df8297cec74
# ╠═b6b7d4ab-0d75-49f1-bfd2-766fe60bf585
# ╠═00c71b09-25bc-41a3-b621-929caaddd626
# ╠═c8bf9933-3676-4b6f-b1c4-9fe73036b8f2
# ╠═db13c2ed-535b-4beb-9c71-bdb7f9fbd2c1
# ╠═b6353c16-9a6c-474c-84c2-dc8baf0c2dfa
# ╠═d1c7141a-1f18-4079-ada2-463593c5e20e
# ╠═f0f3e0a6-6a6a-4e95-b889-6a8995faf266
# ╠═9353df52-9492-44b6-bb80-bdaec11ffb27
# ╠═f5900835-ec28-427b-837f-09620735976b
# ╠═785df551-c770-421b-a2c1-91d064b47de5
# ╠═0d365351-e12c-43e2-87a7-6f46c391aaa2
# ╠═2643efff-0077-4e37-a283-146194a56a75
# ╠═73ff1c94-7317-4a56-a96c-3dd9c6cf69a8
# ╠═899c2928-2457-4748-9fe6-3f60c85a8889
# ╠═995dd661-1dac-4f62-b91c-18f2310c797d
# ╠═811eeb28-c626-414d-b6fb-41dce9de1d05
# ╠═e24e87ee-f6dc-4dac-ab76-c24d2deb4e3f
# ╠═5c33c5db-0583-4c66-92bd-c3cb985d30f7
# ╠═270136ac-f0d2-447b-b4a5-da44709ee79a
# ╠═ea7e8b85-7a18-4f9f-b1e0-60643ea44a29
# ╠═c546afb2-f5cd-4290-ae6d-2805daa2a3fa
# ╠═4c76b220-b46a-4ba2-88f0-3bb787b802c2
# ╠═6a750201-47f6-4fcc-a551-7884c0280bc5
# ╠═988b5474-d301-4527-b553-cebe8decbec5
# ╠═35c2ad9f-a206-45e3-8258-5e04c645dd13
# ╠═77e38166-f9a8-4a00-b034-ddc69a52a345
# ╠═60e2d0eb-07dc-43ee-9c7b-67c46765a31a
# ╠═ffb9d0d1-bca6-4b9a-87a7-e8c6280488fa
# ╠═d3fd5874-8de1-451c-89b1-3dc05a65cef4
# ╠═73eed3c2-8dbe-45bf-85af-5cc48ed6667d
# ╠═0d3b2f02-e66b-469e-96e9-1bb702466f87
# ╠═37ca7c05-f216-4302-99dd-c41b295e100c
# ╠═c958c249-5c8b-4e59-befe-e69b98d84cf5
# ╠═e0a4575c-afd3-4619-8fad-e321d573d1ec
# ╠═4853838b-7292-49ad-b3d8-f116a69da7d6
# ╠═78d16892-dab9-4b2a-97b9-06e49592f0c6
# ╠═887e37cd-8e50-44c3-b239-1936a9b0baef
# ╠═25650a53-6284-4c57-a590-71283d15a3fd
# ╠═ede62a96-39b6-4d2b-ab8a-528048c56556
# ╠═e7e8e6a4-f101-4f55-b8bf-fb9ffd5fa023
# ╠═f1f036c6-d455-499e-b7c6-2a6716cd3485
# ╠═6de3499c-3cf8-4985-9a2b-e69f173044a8
# ╠═f32501b3-a66c-4c3d-aff0-ac28f05abc8d
# ╠═d208210b-ffd9-4bd1-9cd3-c828bcdef732
# ╠═bd5f963d-e833-41a0-a46f-f014f47dd1ca
# ╠═6a4318cf-6897-4ff0-9d82-cde4ab6c9e40
# ╠═a77652cb-6ab9-4f85-b8a2-5a21607444cd
# ╠═7987b9c8-275b-4c74-baeb-1a5c5ca42fde
# ╠═020c589a-ac92-4031-bf29-f482c525caf6
# ╠═6589dac6-db2d-4e2a-82c8-11524df319db
# ╠═fd487373-3d9c-4f79-9a30-3765e7ee242f
# ╠═83cd563c-5743-4a89-8fea-1d512cd2939a
# ╠═0397d23d-6c74-484b-9c8a-3989bab21239
# ╠═d7729d46-79ac-48c2-948a-41c91f78e2e3
# ╠═ec09be7e-90d7-48ee-b7d2-9125c8ca7329
# ╠═b59703c8-2851-494b-b22d-d626ebf9a2e3
# ╠═4478f154-a096-4bdb-a758-c1e300a1de45
# ╠═a7d23b8d-0854-4751-a216-d519e85e2d8b
# ╠═63fcc371-d276-4e15-b812-da1974b22df3
# ╠═dd12ed13-1cd4-48d8-bf28-25baa2a090cd
# ╠═ad4d3db9-43f6-48ee-84ba-a5a360cb258b
# ╠═e5c0d31f-7a7f-4093-bcde-1e66fd50b5c7
# ╠═dfe62938-b4e9-460a-accb-9f852509f038
# ╠═98207c81-9f05-476f-846e-e819398ca81a
# ╠═d0f9e26b-2a3c-42a4-a826-7c8694f5d470
# ╠═8756b4ca-011b-4a6a-aff6-22488c8ac854
# ╠═f0ace932-f1d1-4e2b-9a6d-88b2cea4d849
# ╠═590f5634-843c-4589-ae4d-298293581732
# ╠═9bae1645-1485-411d-a6d3-61906c7c4194
# ╠═ffce2732-1271-4e7e-b3aa-4000f0355dff
# ╠═7dad709c-17e9-4738-b6c0-6d0c9bd70f09
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
