### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ 212d7244-cad6-11f0-8084-f5e1b60d0b39
begin
	using Markdown
	using InteractiveUtils
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
	using NearestNeighbors
	import Glob
end

# ╔═╡ e7f260ce-6c29-45c6-a66d-a6d7e2d306ca
begin
	cmdir=joinpath(ENV["DATA"], "HD5t/itaca")
	pdir =joinpath(ENV["PROJECTS"], "Petit")
end

# ╔═╡ 10f1d6b3-7a5b-46ef-8129-f29715cd5dc3
begin
      using Pkg
      Pkg.activate(pdir)
      Pkg.instantiate()
  end

# ╔═╡ 6cef6b7d-0f42-4e0b-b949-be4b79426a81
PlutoUI.TableOfContents(title="Itaca analysis", indent=true)

# ╔═╡ 7a4f0352-ed82-4674-9530-64151d245c44
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

# ╔═╡ c2cad4cd-3800-4d81-98ae-59f4effe5ac6
begin
	jn = ingredients(string(pdir,"/src/Petit.jl"))
end

# ╔═╡ b53a30f0-fda8-44b0-bed2-c663af6f01b1
md"""
## Functions
"""

# ╔═╡ 9c0e82c3-1d71-4164-bf94-a356538eee8d
function voxel_size(dx, dy, p, l)
	sqrt(l)*(dx + dy) /(2.0*sqrt(p))
end

# ╔═╡ 65f44517-5794-4182-aae0-4e2385a1b74b
function consecutive_voxel_distances(df::DataFrame)
      n = nrow(df)
      n < 2 && return Float64[]

      [sqrt((df.x[i+1] - df.x[i])^2 +
            (df.y[i+1] - df.y[i])^2 +
            (df.z[i+1] - df.z[i])^2) for i in 1:(n-1)]
  end

# ╔═╡ dc25b509-24f9-4e92-aa86-42e59585b8c0
  function closest_voxel_distances_fast(df::DataFrame)
      n = nrow(df)
      n < 2 && return Float64[]

      points = reduce(hcat, [[df.x[i], df.y[i], df.z[i]] for i in 1:n])
      tree = KDTree(points)
      _, dists = knn(tree, points, 2, true)

      return [d[2] for d in dists]
  end

# ╔═╡ 94473bfe-81b7-4f39-b5e4-60eb4a840a9a
function select_single_tracks(TRKS)
	ES = Float64[]
	NT=Integer[]
	TRKS1 = []
	for (iv,trks) in enumerate(TRKS)
		#println("event =$(iv), event id = $(trks[1].voxels.event_id[1]) number of tracks = $(length(trks))")
		push!(NT,length(trks) )
		if length(trks) ==1
			push!(TRKS1, trks)
		else
			for (i,trk) in enumerate(trks)
				#println("track =$i")
				#println("track length =$(length(trk.voxels.x))")
				if length(trk.voxels.x) == 1
					es = trk.voxels.energy[1]*1e+3
					push!(ES, es)
				end
			end
		end
	end
	return TRKS1, NT, ES
end

# ╔═╡ e91c9f90-63a6-48cf-ab8a-3d8d4ac320c4
function load_data(input_file)
	input_path = joinpath(cmdir, input_file)
	dfs = jn.Petit.get_dataset_dfs(input_path)
	hitsdf = dfs["hits"]

    # Count events from loaded data
    println("Counting events from loaded data...")
    ntot = length(unique(hitsdf.event_id))
    println("Number of events with hits: $ntot")
	return hitsdf
end

# ╔═╡ a1694272-6f8d-43e1-a3ed-c575a8e93f7f
function select_events(hitsdf; iev, lev, voxel_size_mm, energy_threshold_kev,
					  emin, emax, nprint=10)

	max_distance_mm = 3.0 * voxel_size_mm

	println("++voxel_size_mm =$voxel_size_mm")
    println("++max_distance_mm =$max_distance_mm")
    println("++energy_threshold_kev =$energy_threshold_kev")
    println("++emin =$emin")
    println("++emax =$emax")
	
	TRKS = []
	ENE = Float64[]
	for ievt in iev:lev
		bbt1mm = jn.Petit.select_events(hitsdf, ievt;
                                 voxel_size_mm=voxel_size_mm,
                                 max_distance_mm=max_distance_mm,
                                 energy_threshold_kev=energy_threshold_kev,
                                 emin=emin,
                                 emax=emax)
		if ievt % nprint == 0
         	 println("++nevent =$ievt")
			println("number of tracks found = $(length(bbt1mm))")
      	end

		ene = 0.0
		for tt in bbt1mm
			ene+=sum(tt.voxels.energy)
		end
		push!(ENE, ene*1e+3)
		push!(TRKS, bbt1mm)
	end
	TRKS, ENE
end
		

# ╔═╡ 84eb8703-228b-46f2-8868-5086f576572e
function print_reco(trk, walk, rblob)
	extremes, _, _, xtrack_length, confidence = walk
	xstart_voxel, xend_voxel = extremes
	blobs = jn.Petit.energy_in_spheres_around_extremes(trk, walk, rblob)
	eb1 = blobs.blob1_energy * 1e+3
	nb1 = blobs.blob1_voxel_count
	eb2 = blobs.blob2_energy * 1e+3
	nb2 = blobs.blob2_voxel_count
	md"""
	#### Find blobs: xe137 with radius $(rblob) 
	- confidence = $(confidence)
	- start voxel: x = $(xstart_voxel.x), y = $(xstart_voxel.y), z = $(xstart_voxel.z)
	- end voxel: x = $(xend_voxel.x), y = $(xend_voxel.y), z = $(xend_voxel.z)
	- track length L =$(xtrack_length)
	- blob 1 energy = $(round(eb1, digits=1)) keV
	- blob 2 energy = $(round(eb2, digits=1)) keV
	- blob 1 # of voxels = $(nb1)
	- blob 2 # of voxels = $(nb2)
	"""
end

# ╔═╡ fc9cc4a7-56b5-4852-bc19-483fa7379ab2
md"""
## Analysis 15 bar
"""

# ╔═╡ 132042e7-f7ad-46ea-88d6-396bd072c5a8
md"""
### Size of voxels
"""

# ╔═╡ 4fd05d01-16c4-4849-8486-ccfad7e3290d
begin
	vxe = voxel_size(0.9, 3.5, 15.0, 200.0)
	md"""
	- voxel size for pure xenon = $@sprintf("%.1e", vxe) mm
	"""
end

# ╔═╡ 7dea81cf-f51d-4411-ae68-bb86bb04dcf3
begin
	vhe = voxel_size(0.75, 1.6, 15.0, 100.0)
	md"""
	- voxel size for 10 % He = $@sprintf("%.1e", vhe) mm
	"""
end

# ╔═╡ 6a0c75ee-563a-4034-aa07-ad1e11ef1d48
md"""
### blob radius
"""

# ╔═╡ 880b9287-6ba9-44cb-8036-4b46a8780177
begin
	rblob5=5.0
	rblob10 = 10.0
end

# ╔═╡ 7b3a8a4d-f649-4175-b881-453b7e44e598
md"""
### Load data
"""

# ╔═╡ d0fecc3e-5e0e-4522-9efa-b4b60fd151a5
bidf = load_data("gammas_2458keV_15bar_100mum.next.h5")

# ╔═╡ b32e2b11-788e-4017-b5cf-4e4d3968398d
bbdf = load_data("0nubb_15bar_100mum.next.h5")

# ╔═╡ 6950e911-3003-4778-98bd-2b1a0a819844
xedf = load_data("electrons_2400_2500_15bar_100mum.next.h5")

# ╔═╡ 41653083-87fb-46f7-8d1f-5cf9d34c2ea7
begin
	# Event number to look at
	nevent = 0
end

# ╔═╡ 47dc42ac-9dd9-4318-9a29-7655541d3af8
md"""
#### bb0nu example, pure MC data
"""

# ╔═╡ 83a718af-c2d9-4743-827d-8682a7379c65
begin
	bbevt = jn.Petit.get_event(bbdf, nevent)
	jn.Petit.plot_event(bbevt)
end

# ╔═╡ 7ba1eacc-88f9-4bc3-ba75-1efce41dadb7
let
	hd, pd = jn.Petit.step_hist(consecutive_voxel_distances(bbevt);
	                                                 nbins = 40,
										             xlim   = (0.0, 1.0),
	                                                 xlabel = " d (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" Distance voxels (mm)")
	he, pe = jn.Petit.step_hist(bbevt.energy*1e+3;
	                                                 nbins = 40,
										             #xlim   = (0.0, 50.0),
	                                                 xlabel = " E  (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" E voxel (keV)")
	plot(pd, pe)
end

# ╔═╡ 38474add-cbdf-4772-ba95-c7da00485b31
md"""
#### Xe137 example, pure MC data
"""

# ╔═╡ 86bf06e8-9caa-4c18-80ff-160cbd945986
begin
	bievt = jn.Petit.get_event(xedf, nevent)
	jn.Petit.plot_event(bievt)
end

# ╔═╡ 3fb45e85-18cb-455c-b034-b67f8513384a
let
	hd, pd = jn.Petit.step_hist(consecutive_voxel_distances(bidf);
	                                                 nbins = 40,
										             xlim   = (0.0, 1.0),
	                                                 xlabel = " d (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" Distance voxels (mm)")
	he, pe = jn.Petit.step_hist(bidf.energy*1e+3;
	                                                 nbins = 40,
										             #xlim   = (0.0, 50.0),
	                                                 xlabel = " E  (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" E voxel (keV)")
	plot(pd, pe)
end

# ╔═╡ 8c173831-8e37-47d4-a183-d95be8773781
md"""
#### bb0nu example, voxels = 1mm
"""

# ╔═╡ d92a0d53-d5f0-4e77-9911-a97e93d5b7b4
begin
	bbv1mm = jn.Petit.voxelize_event(bbdf, nevent, 1.0)
	jn.Petit.plot_event(bbv1mm)
end

# ╔═╡ 972b17ca-810f-4279-b945-7e7acc6c2f1f
let
	
	hd, pd = jn.Petit.step_hist(closest_voxel_distances_fast(bbv1mm);
	                                                 nbins = 40,
										             #xlim   = (0.0, 1.0),
	                                                 xlabel = " d (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" Distance voxels (mm)")
	he, pe = jn.Petit.step_hist(bbv1mm.energy*1e+3;
	                                                 nbins = 40,
										             #xlim   = (0.0, 50.0),
	                                                 xlabel = " E  (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" E voxel (keV)")
	plot(pd, pe)
	
end

# ╔═╡ 3bc7bdae-e234-4c92-858c-b00fbd8ce861
md"""
#### bb0nu example, voxels = 2mm
"""

# ╔═╡ fa3f37e2-d477-4c0b-942c-00c6d649601b
begin
	bbv2mm = jn.Petit.voxelize_event(bbdf, nevent, 2.0)
	jn.Petit.plot_event(bbv2mm)
end

# ╔═╡ 0c01706f-e6fb-46ca-842c-0150c01fec10
let
	
	hd, pd = jn.Petit.step_hist(closest_voxel_distances_fast(bbv2mm);
	                                                 nbins = 40,
										             xlim   = (0.0, 10.),
	                                                 xlabel = " d (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" Distance voxels (mm)")
	he, pe = jn.Petit.step_hist(bbv2mm.energy*1e+3;
	                                                 nbins = 40,
										             xlim   = (0.0, 50.0),
	                                                 xlabel = " E  (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" E voxel (keV)")
	plot(pd, pe)
end

# ╔═╡ 74e35417-b378-4009-8270-c7502fd98721
md"""
#### bb0nu example, voxels = 4mm
"""

# ╔═╡ 169dd459-9acc-4631-b8eb-9e75dd79810a
begin
	bbv4mm = jn.Petit.voxelize_event(bbdf, nevent, 4.0)
	jn.Petit.plot_event(bbv4mm)
end

# ╔═╡ 71102ff9-eac0-4833-b2ee-6f5c570a175d
let
	
	hd, pd = jn.Petit.step_hist(closest_voxel_distances_fast(bbv4mm);
	                                                 nbins = 40,
										             xlim   = (0.0, 5.0),
	                                                 xlabel = " d (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" Distance voxels (mm)")
	he, pe = jn.Petit.step_hist(bbv4mm.energy*1e+3;
	                                                 nbins = 40,
										             xlim   = (0.0, 50.0),
	                                                 xlabel = " E  (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" E voxel (keV)")
	plot(pd, pe)
end

# ╔═╡ a440b111-9440-4a0b-90f9-b53d0a13c27e
md"""
#### bb0nu example, voxels = 8mm
"""

# ╔═╡ e2686385-dd8b-4a43-afff-60bb53ef9595
begin
	bbv8mm = jn.Petit.voxelize_event(bbdf, nevent, 8.0)
	jn.Petit.plot_event(bbv8mm)
end

# ╔═╡ 9fccc980-fa99-4e6c-acf2-9ed71ffc8685
let
	
	hd, pd = jn.Petit.step_hist(closest_voxel_distances_fast(bbv8mm);
	                                                 nbins = 20,
										             xlim   = (0.0, 10.0),
	                                                 xlabel = " d (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" Distance voxels (mm)")
	he, pe = jn.Petit.step_hist(bbv8mm.energy*1e+3;
	                                                 nbins = 40,
										             xlim   = (0.0, 150.0),
	                                                 xlabel = " E  (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" E voxel (keV)")
	plot(pd, pe)
end

# ╔═╡ 105bd81d-963d-4ac1-9bfc-1f3179f90202
md"""
### Select bb0nu and Xe137 events, returning tracks
"""

# ╔═╡ 18308109-43e7-463a-b004-ec3da05303a8
bbtrks1mm, ene_bbtrks1mm = select_events(bbdf; iev=1, lev=100, 
						  voxel_size_mm=1.0, 
						  energy_threshold_kev=1.0,
					  	  emin=2400.0, emax=2500.0, 
						  nprint=20)

# ╔═╡ 26e87126-da00-4894-9aa9-eae9deef4666
begin
xetrks1mm, ene_xetrks1mm = select_events(xedf; iev=1, lev=100, 
						  voxel_size_mm=1.0, 
						  energy_threshold_kev=1.0,
					  	  emin=2400.0, emax=2600.0, 
						  nprint=20)
end

# ╔═╡ 46bc27bd-e06c-4ed0-97dd-3e31879b3535
begin
hene_xetrks1mm, pene_xetrks1mm = jn.Petit.step_hist(ene_xetrks1mm;
	                                            nbins = 40,
										    	xlim   = (2400.0, 2600.0),
	                                            xlabel = " energy of tracks",
	                                            ylabel = "Frequency",
	                                            title=" Xe-137 energy of tracks")
hene_bbtrks1mm, pene_bbtrks1mm = jn.Petit.step_hist(ene_bbtrks1mm;
	                                            nbins = 40,
										    	xlim   = (2400.0, 2600.0),
	                                            xlabel = " energy of tracks",
	                                            ylabel = "Frequency",
	                                            title=" Xe-137 energy of tracks")

	
	plot(pene_bbtrks1mm, pene_xetrks1mm)
end

# ╔═╡ 085ebb25-027a-4027-9f4c-a69013aee707
md"""
### Select single-tracks for bb0nu and Bi214
"""

# ╔═╡ e7e5b39a-c679-46a2-b78c-0b48db8a4027
begin
	bbstrk1mm, ntstrk1mm, egstrk1mm = select_single_tracks(bbtrks1mm)
	println(" bb0nu: # of events with a single track = $(count(==(1), ntstrk1mm))")
end

# ╔═╡ a793ebb6-0ee9-4117-a3e8-717063c2e8a0
let
	
	he, pe = jn.Petit.step_hist(egstrk1mm;
	                                            nbins = 40,
										        #xlim   = (0.0, 1.0),
	                                            xlabel = " E (keV)",
	                                            ylabel = "Frequency",
	                                            title=" bb0nu E satellite γ's (keV)")
	hnt, pnt = jn.Petit.step_hist(ntstrk1mm;
	                                            nbins = 40,
										    	#xlim   = (0.0, 1.0),
	                                            xlabel = " # of tracks",
	                                            ylabel = "Frequency",
	                                            title=" bb0nu: # of tracks")
	
	plot(pe, pnt)
end

# ╔═╡ f2f7a839-6c22-40fd-98ce-dbc3322f194a
begin
	xe_strk1mm, xe_ntstrk1mm, xe_egstrk1mm = select_single_tracks(xetrks1mm)
	println(" Bi214: # of events with a single track = $(count(==(1), xe_ntstrk1mm))")
end

# ╔═╡ 285e8019-e433-4837-ab5d-edea5ad5031c
let
	
	he, pe = jn.Petit.step_hist(xe_egstrk1mm;
	                                            nbins = 40,
										        #xlim   = (1.0, 10.0),
	                                            xlabel = " E (keV)",
	                                            ylabel = "Frequency",
	                                            title=" Xe-137 electrons (keV)")
	hnt, pnt = jn.Petit.step_hist(xe_ntstrk1mm;
	                                            nbins = 40,
										    	xlim   = (1.0, 10.0),
	                                            xlabel = " # of tracks",
	                                            ylabel = "Frequency",
	                                            title=" Xe-137 electrons : # of tracks")
	
	plot(pe, pnt)
end

# ╔═╡ e56d5fa6-02d9-4260-aca5-812164409a48
md"""
### Blobs
"""

# ╔═╡ dafac731-bbfa-4bec-b6cc-26e7bf180b0b
begin
	bbstrk1mm[1]
end

# ╔═╡ 56a621ff-cb18-4356-9dfb-c4dae25cca2e
bbstrk1mm[1][1]

# ╔═╡ 4067b243-9679-4402-8570-5f351b489805
ievent=2

# ╔═╡ 07c6b6cd-681b-498b-9448-4d49014820f9
typeof(bbstrk1mm[ievent][1])

# ╔═╡ 985ff432-1d22-4e22-86e1-6382eaee243b
begin
	bbtrk = bbstrk1mm[ievent][1]
	xout =jn.Petit.walk_track_from_extremes(bbtrk)
	print_reco(bbtrk, xout, rblob5)
end

# ╔═╡ d8e42ebf-e7ec-42fb-9867-3af552ffb842
jn.Petit.plot_track_with_extremes(bbtrk, xout;
                                  show_connections=false,
                                  alpha_connections=0.3)

# ╔═╡ b4ce8dcb-52fd-462a-b864-4e2b11a11990
begin
	p2d, p3d = jn.Petit.plot_track_blobs(bbtrk, rblob5;
                         markersize_voxels=3.0,
                         show_connections=true,
                         alpha_connections=0.2,
                         alpha_spheres=0.3,
                         sphere_resolution=20)
	plot(p2d)
end

# ╔═╡ 975b6fa8-fac9-48ea-9e31-cf7d397ebf1c
plot(p3d, size=(800, 600))

# ╔═╡ ae4ab185-1685-4d3f-bda2-c359cf9ca91a


# ╔═╡ dd011045-aca5-49a9-8227-aa66097f98ca
begin
	xetrk = xe_strk1mm[ievent][1]
	xoutxe =jn.Petit.walk_track_from_extremes(xetrk)
	print_reco(xetrk, xoutxe, rblob5)
end
	

# ╔═╡ fec4dcbb-adab-4435-ad6a-6d2b51b5d205


# ╔═╡ d656e40b-c9e7-4be2-91eb-8f56fbca9495
jn.Petit.plot_track_with_extremes(xetrk, xoutxe;
                                  show_connections=false,
                                  alpha_connections=0.3)

# ╔═╡ 7db0dde5-2760-47aa-a3c4-e5fb93e3f36b
begin
	x2d, x3d = jn.Petit.plot_track_blobs(xetrk, rblob5;
                         markersize_voxels=3.0,
                         show_connections=true,
                         alpha_connections=0.2,
                         alpha_spheres=0.3,
                         sphere_resolution=20)
	plot(x2d)
end

# ╔═╡ ef541827-3bcc-497d-8a1c-9a2d05d5e75a
plot(x3d)

# ╔═╡ Cell order:
# ╠═e7f260ce-6c29-45c6-a66d-a6d7e2d306ca
# ╠═10f1d6b3-7a5b-46ef-8129-f29715cd5dc3
# ╠═212d7244-cad6-11f0-8084-f5e1b60d0b39
# ╠═6cef6b7d-0f42-4e0b-b949-be4b79426a81
# ╠═7a4f0352-ed82-4674-9530-64151d245c44
# ╠═c2cad4cd-3800-4d81-98ae-59f4effe5ac6
# ╠═b53a30f0-fda8-44b0-bed2-c663af6f01b1
# ╠═9c0e82c3-1d71-4164-bf94-a356538eee8d
# ╠═65f44517-5794-4182-aae0-4e2385a1b74b
# ╠═dc25b509-24f9-4e92-aa86-42e59585b8c0
# ╠═94473bfe-81b7-4f39-b5e4-60eb4a840a9a
# ╠═e91c9f90-63a6-48cf-ab8a-3d8d4ac320c4
# ╠═a1694272-6f8d-43e1-a3ed-c575a8e93f7f
# ╠═84eb8703-228b-46f2-8868-5086f576572e
# ╠═fc9cc4a7-56b5-4852-bc19-483fa7379ab2
# ╠═132042e7-f7ad-46ea-88d6-396bd072c5a8
# ╠═4fd05d01-16c4-4849-8486-ccfad7e3290d
# ╠═7dea81cf-f51d-4411-ae68-bb86bb04dcf3
# ╠═6a0c75ee-563a-4034-aa07-ad1e11ef1d48
# ╠═880b9287-6ba9-44cb-8036-4b46a8780177
# ╠═7b3a8a4d-f649-4175-b881-453b7e44e598
# ╠═d0fecc3e-5e0e-4522-9efa-b4b60fd151a5
# ╠═b32e2b11-788e-4017-b5cf-4e4d3968398d
# ╠═6950e911-3003-4778-98bd-2b1a0a819844
# ╠═41653083-87fb-46f7-8d1f-5cf9d34c2ea7
# ╠═47dc42ac-9dd9-4318-9a29-7655541d3af8
# ╠═83a718af-c2d9-4743-827d-8682a7379c65
# ╠═7ba1eacc-88f9-4bc3-ba75-1efce41dadb7
# ╠═38474add-cbdf-4772-ba95-c7da00485b31
# ╠═86bf06e8-9caa-4c18-80ff-160cbd945986
# ╠═3fb45e85-18cb-455c-b034-b67f8513384a
# ╠═8c173831-8e37-47d4-a183-d95be8773781
# ╠═d92a0d53-d5f0-4e77-9911-a97e93d5b7b4
# ╠═972b17ca-810f-4279-b945-7e7acc6c2f1f
# ╠═3bc7bdae-e234-4c92-858c-b00fbd8ce861
# ╠═fa3f37e2-d477-4c0b-942c-00c6d649601b
# ╠═0c01706f-e6fb-46ca-842c-0150c01fec10
# ╠═74e35417-b378-4009-8270-c7502fd98721
# ╠═169dd459-9acc-4631-b8eb-9e75dd79810a
# ╠═71102ff9-eac0-4833-b2ee-6f5c570a175d
# ╠═a440b111-9440-4a0b-90f9-b53d0a13c27e
# ╠═e2686385-dd8b-4a43-afff-60bb53ef9595
# ╠═9fccc980-fa99-4e6c-acf2-9ed71ffc8685
# ╠═105bd81d-963d-4ac1-9bfc-1f3179f90202
# ╠═18308109-43e7-463a-b004-ec3da05303a8
# ╠═26e87126-da00-4894-9aa9-eae9deef4666
# ╠═46bc27bd-e06c-4ed0-97dd-3e31879b3535
# ╠═085ebb25-027a-4027-9f4c-a69013aee707
# ╠═e7e5b39a-c679-46a2-b78c-0b48db8a4027
# ╠═a793ebb6-0ee9-4117-a3e8-717063c2e8a0
# ╠═f2f7a839-6c22-40fd-98ce-dbc3322f194a
# ╠═285e8019-e433-4837-ab5d-edea5ad5031c
# ╠═e56d5fa6-02d9-4260-aca5-812164409a48
# ╠═dafac731-bbfa-4bec-b6cc-26e7bf180b0b
# ╠═07c6b6cd-681b-498b-9448-4d49014820f9
# ╠═56a621ff-cb18-4356-9dfb-c4dae25cca2e
# ╠═4067b243-9679-4402-8570-5f351b489805
# ╠═985ff432-1d22-4e22-86e1-6382eaee243b
# ╠═d8e42ebf-e7ec-42fb-9867-3af552ffb842
# ╠═b4ce8dcb-52fd-462a-b864-4e2b11a11990
# ╠═975b6fa8-fac9-48ea-9e31-cf7d397ebf1c
# ╠═ae4ab185-1685-4d3f-bda2-c359cf9ca91a
# ╠═dd011045-aca5-49a9-8227-aa66097f98ca
# ╠═fec4dcbb-adab-4435-ad6a-6d2b51b5d205
# ╠═d656e40b-c9e7-4be2-91eb-8f56fbca9495
# ╠═7db0dde5-2760-47aa-a3c4-e5fb93e3f36b
# ╠═ef541827-3bcc-497d-8a1c-9a2d05d5e75a
