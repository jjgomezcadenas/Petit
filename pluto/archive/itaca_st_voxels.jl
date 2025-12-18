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

# ╔═╡ c5f2e0aa-a7ee-406c-be13-e6219c409638


# ╔═╡ 84eb8703-228b-46f2-8868-5086f576572e
function print_reco(trk, walk, rblob, label)
	extremes, _, _, xtrack_length, confidence = walk
	xstart_voxel, xend_voxel = extremes
	blobs = jn.Petit.energy_in_spheres_around_extremes(trk, walk, rblob)
	eb1 = blobs.blob1_energy * 1e+3
	nb1 = blobs.blob1_voxel_count
	eb2 = blobs.blob2_energy * 1e+3
	nb2 = blobs.blob2_voxel_count
	db = abs((eb1 -eb2))/(eb1+eb2)
	md"""
	#### Find blobs: $(label) with radius $(rblob) 
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

# ╔═╡ 0289428e-6db3-4bfc-82fe-5b9e71dd3706
function print_reco_blobs(trk, walk, rblob, label)
	extremes, _, _, xtrack_length, confidence = walk
	blobs = jn.Petit.energy_in_spheres_around_extremes(trk, walk, rblob)
	eb1 = blobs.blob1_energy * 1e+3
	eb2 = blobs.blob2_energy * 1e+3
	db = abs((eb1 -eb2))/(eb1+eb2)
	md"""
	#### Find blobs: $(label) with radius $(rblob) 
	- confidence = $(confidence)
	- track length L =$(xtrack_length)
	- blob 1 energy = $(round(eb1, digits=1)) keV
	- blob 2 energy = $(round(eb2, digits=1)) keV
	- blob asymmetry = $(db)
	"""
end

# ╔═╡ 726badaf-c389-484e-a46d-cdda4d4d578d
function get_energy_blobs(trk, walk, rblob)
	extremes, _, _, xtrack_length, confidence = walk
	xstart_voxel, xend_voxel = extremes
	blobs = jn.Petit.energy_in_spheres_around_extremes(trk, walk, rblob)
	eb1 = blobs.blob1_energy * 1e+3
	nb1 = blobs.blob1_voxel_count
	eb2 = blobs.blob2_energy * 1e+3
	nb2 = blobs.blob2_voxel_count
	return eb1, eb2
end

# ╔═╡ 67ac7cde-e9a4-4013-8e30-61ac29d28175


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

# ╔═╡ 105bd81d-963d-4ac1-9bfc-1f3179f90202
md"""
### Select bb0nu and Xe137 events, returning tracks
"""

# ╔═╡ 7364c781-7714-46f8-b412-2dd04ee4b8a5
begin
	trksbb1mm = jn.Petit.get_bb0nu_tracks(cmdir; tag="*v1mm*")
	bbekev = [sum(tt.voxels.energy) * 1e+3 for tt in trksbb1mm.tracks]
	bbevox = [e * 1e+3 for tt in trksbb1mm.tracks for e in tt.voxels.energy]
	bbxvox = [x for tt in trksbb1mm.tracks for x in tt.voxels.x]
	bbyvox = [x for tt in trksbb1mm.tracks for x in tt.voxels.y]
	bbzvox = [x for tt in trksbb1mm.tracks for x in tt.voxels.z]
	md"""
	read $(trksbb1mm.n1trk) tracks
	"""
end

# ╔═╡ caa5c1d8-5812-42f4-892a-34c14599b3e0
begin
	trksxe1mm = jn.Petit.get_xe137_tracks(cmdir; xedir="xe137", tag="*v1mm*")
	xeekev = [sum(tt.voxels.energy) * 1e+3 for tt in trksxe1mm.tracks]
	xeevox = [e * 1e+3 for tt in trksxe1mm.tracks for e in tt.voxels.energy]
	xexvox = [x for tt in trksxe1mm.tracks for x in tt.voxels.x]
	xeyvox = [x for tt in trksxe1mm.tracks for x in tt.voxels.y]
	xezvox = [x for tt in trksxe1mm.tracks for x in tt.voxels.z]
	md"""
	read $(trksxe1mm.n1trk) tracks
	"""
end

# ╔═╡ d089e072-9ae9-4737-9e1e-44d84e354d28
let
	
	_, pekev = jn.Petit.step_hist(bbekev;
	                                                 nbins = 40,
										             xlim   = (2400.0, 2500.0),
	                                                 xlabel = " E (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" bb track energy")
	_, pevox = jn.Petit.step_hist(bbevox;
	                                                 nbins = 40,
										             xlim   = (0.0, 100.0),
	                                                 xlabel = " E  (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" bb energy voxel (keV)")
	_, pxvox = jn.Petit.step_hist(bbxvox;
	                                                 nbins = 40,
										             #xlim   = (0.0, 100.0),
	                                                 xlabel = " X  (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" bb X (mm)")
	_, pyvox = jn.Petit.step_hist(bbyvox;
	                                                 nbins = 40,
										             #xlim   = (0.0, 100.0),
	                                                 xlabel = " Y  (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" bb Y (mm)")
	_, pzvox = jn.Petit.step_hist(bbzvox;
	                                                 nbins = 20,
										             #xlim   = (0.0, 100.0),
	                                                 xlabel = " Z  (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" bb Z (mm)")
	pxyvox = scatter(bbxvox, bbyvox, markersize=2, markercolor=:black,
      xlabel="x (mm)",
      ylabel="y (mm)",
      title="bb Voxel positions",
      label=nothing)

	plot(pekev, pevox, pxyvox, pzvox)
	
end

# ╔═╡ cf81a428-8cb7-4c19-8e58-f87778e867fd
let
	
	_, pekev = jn.Petit.step_hist(xeekev;
	                                                 nbins = 40,
										             xlim   = (2400.0, 2500.0),
	                                                 xlabel = " E (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" xe track energy")
	_, pevox = jn.Petit.step_hist(xeevox;
	                                                 nbins = 40,
										             xlim   = (0.0, 100.0),
	                                                 xlabel = " E  (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" xe energy voxel (keV)")
	_, pxvox = jn.Petit.step_hist(xexvox;
	                                                 nbins = 40,
										             #xlim   = (0.0, 100.0),
	                                                 xlabel = " X  (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" xe X (mm)")
	_, pyvox = jn.Petit.step_hist(xeyvox;
	                                                 nbins = 40,
										             #xlim   = (0.0, 100.0),
	                                                 xlabel = " Y  (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" xe Y (mm)")
	_, pzvox = jn.Petit.step_hist(xezvox;
	                                                 nbins = 20,
										             #xlim   = (0.0, 100.0),
	                                                 xlabel = " Z  (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" xe Z (mm)")
	pxyvox = scatter(xexvox, xeyvox, markersize=2, markercolor=:black,
      xlabel="x (mm)",
      ylabel="y (mm)",
      title="bb Voxel positions",
      label=nothing)

	plot(pekev, pevox, pxyvox, pzvox)
	
end

# ╔═╡ 085ebb25-027a-4027-9f4c-a69013aee707
md"""
### Select single-tracks for bb0nu and Bi214
"""

# ╔═╡ e56d5fa6-02d9-4260-aca5-812164409a48
md"""
### Blobs
"""

# ╔═╡ e12c925b-101b-4829-bb58-e98debb454a9
trksbb1mm.n1trk

# ╔═╡ 4067b243-9679-4402-8570-5f351b489805
begin
	ievent=1
	rblob = 5.0
end

# ╔═╡ 985ff432-1d22-4e22-86e1-6382eaee243b
begin
	#bbtrk = bbstrk1mm[ievent][1]
	bbtrk = trksbb1mm.tracks[ievent]
	
	xout =jn.Petit.walk_track_from_extremes(bbtrk)
	print_reco_blobs(bbtrk, xout, rblob, "bb0nu")
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

# ╔═╡ dd011045-aca5-49a9-8227-aa66097f98ca
begin
	xetrk = trksxe1mm.tracks[ievent] 
	xoutxe =jn.Petit.walk_track_from_extremes(xetrk)
	print_reco_blobs(xetrk, xoutxe, rblob, "xe137")
end
	

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

# ╔═╡ d37787c8-8511-40b6-aee5-8064cf5ed268
function blob_asymmetry(trks; i0=1, il=10, r0=5, rl =15)
	ie = il-i0 + 1
	ir = rl - r0 + 1
	DB = Matrix{Float64}(undef, ie, ir)
	
	for ievent in i0:il
		bbtrk = trks[ievent]
		xout =jn.Petit.walk_track_from_extremes(bbtrk)
		
		for (ir,rr) in enumerate(r0:rl)
			eb1, eb2 = get_energy_blobs(bbtrk, xout, rr*1.0)
			db = abs((eb1 -eb2))/(eb1+eb2)
			DB[ievent, ir] = db
			#println("DB[$ievent, $ir]=$db")
		end
	end
	
	vec(mean(DB, dims=1)), vec(std(DB, dims=1))
end


# ╔═╡ 51575ff2-a24a-4177-90fe-c904a20c2d57
begin
	r0 = 5
	rl = 15
	bbAm,bbAs = blob_asymmetry(trksbb1mm.tracks; i0=1, il=100, r0=r0, rl =rl)
	xeAm,xeAs = blob_asymmetry(trksxe1mm.tracks; i0=1, il=100, r0=r0, rl =rl)
	
end

# ╔═╡ c79c2d21-a22d-4ed2-b7a7-3a5f7a6b9b4c
begin 
	xr = collect(r0:rl)
	scatter(xr, bbAm, yerror=bbAs, label="bb0nu blob asymmetry")
  	plot!(xr, bbAm, label=nothing, linewidth=1)
	scatter!(xr, xeAm, yerror=xeAs, label="xe-137 blob asymmetry")
  	plot!(xr, xeAm, label=nothing, linewidth=1)
end

# ╔═╡ e9b013a3-3353-4035-8830-5737d952bf52
function blob_asymmetry_rblob(trks; i0=1, il=10, r0)
	ie = il-i0 + 1
	EB1 = Vector{Float64}(undef, ie)
	EB2 = Vector{Float64}(undef, ie)
	DB = Vector{Float64}(undef, ie)
	ll = min(il, length(trks))
	for ievent in i0:ll
		bbtrk = trks[ievent]
		xout =jn.Petit.walk_track_from_extremes(bbtrk)
	
		eb1, eb2 = get_energy_blobs(bbtrk, xout, r0)
		db = 1.0
		if (eb1+eb2) >0
			db = abs((eb1 -eb2))/(eb1+eb2)
		end
		EB1[ievent] = eb1
		EB2[ievent] = eb2
		DB[ievent] = db
	end
	
	EB1, EB2, DB
end

# ╔═╡ 9a89e70a-f40b-4fda-a17d-8a875a8c327a
begin
	
	bbeb1,bbeb2, bbdb = blob_asymmetry_rblob(trksbb1mm.tracks; i0=1, il=500, r0=5.0)
	xeeb1,xeeb2, xedb = blob_asymmetry_rblob(trksxe1mm.tracks; i0=1, il=500, r0=5.0)
	bbeb1r10,bbeb2r10, _ = blob_asymmetry_rblob(trksbb1mm.tracks; i0=1, il=500, r0=10.0)
	xeeb1r10,xeeb2r10, _ = blob_asymmetry_rblob(trksxe1mm.tracks; i0=1, il=500, r0=10.0)
	
end

# ╔═╡ 79e0190f-054a-4455-9df8-7ba381eb6f89
bbeb1

# ╔═╡ 4995d151-11e6-46b0-a19d-c61ca0de8b8c
bbeb2

# ╔═╡ ccfa2e3c-92d8-4108-a21e-4f09bb1d4878
let
	pbbeb1eb2 = scatter(bbeb1,bbeb2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" bb0nu: Eb1 vs Eb2",
      label=nothing)
	pxeeb1eb2 = scatter(xeeb1,xeeb2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title="Xe137: Eb1 vs Eb2",
      label=nothing)
	pbbeb1eb2r10 = scatter(bbeb1r10,bbeb2r10, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" bb0nu: Eb1 vs Eb2",
      label=nothing)
	pxeeb1eb2r10 = scatter(xeeb1r10,xeeb2r10, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title="Xe137: Eb1 vs Eb2",
      label=nothing)
	_, pbbdb = jn.Petit.step_hist(bbdb;
	                                                 nbins = 40,
										             #xlim   = (0.0, 100.0),
	                                                 xlabel = " ΔEb bb0nu",
	                                                 ylabel = "Frequency",
	                                                 title=" ΔEb bb0nu")
	_, pxedb = jn.Petit.step_hist(xedb;
	                                                 nbins = 40,
										             #xlim   = (0.0, 100.0),
	                                                 xlabel = " ΔEb Xe137",
	                                                 ylabel = "Frequency",
	                                                 title=" ΔEb bb0nu")
	pdbbbxe = scatter(bbdb,xedb, markersize=2, markercolor=:black,
      xlabel="ΔEb bb0nu (keV)",
      ylabel="ΔEb Xe137 (keV)",
      title="ΔEb bb0nu vs Xe137",
      label=nothing)

	plot(pbbeb1eb2, pxeeb1eb2, pbbeb1eb2r10, pxeeb1eb2r10)
end

# ╔═╡ 66fc61a5-3fa0-43aa-a005-d6d0abf7e491
begin
	bbr3df = jn.Petit.merge_csv_files("bb0nu_blobs_1mm_rblob_3")
	bbr5df = jn.Petit.merge_csv_files("bb0nu_blobs_1mm_rblob_5")
	bbr8df = jn.Petit.merge_csv_files("bb0nu_blobs_1mm_rblob_8")
	bbr10df = jn.Petit.merge_csv_files("bb0nu_blobs_1mm_rblob_10")
	xer3df = jn.Petit.merge_csv_files("xe137_blobs_1mm_rblob_3") 
	xer5df = jn.Petit.merge_csv_files("xe137_blobs_1mm_rblob_5") 
	xer8df = jn.Petit.merge_csv_files("xe137_blobs_1mm_rblob_8") 
	xer10df = jn.Petit.merge_csv_files("xe137_blobs_1mm_rblob_10") 
end

# ╔═╡ 1e9da122-9d8a-4b9e-a5df-fac798a5f8a8
let
	pbbr3eb1b2 = scatter(bbr3df.eB1,bbr3df.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" bb0nu (r= 3mm) ",
      label=nothing)
	pxer3eb1b2 = scatter(xer3df.eB1,xer3df.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" Xe137 (r= 3mm) ",
      label=nothing)
	pbbr5eb1b2 = scatter(bbr5df.eB1,bbr5df.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" bb0nu (r= 5mm) ",
      label=nothing)
	pxer5eb1b2 = scatter(xer5df.eB1,xer5df.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" Xe137 (r= 5mm) ",
      label=nothing)
	

	plot(pbbr3eb1b2,pxer3eb1b2, pbbr5eb1b2,pxer5eb1b2 )
end

# ╔═╡ 84208e5e-6209-427b-a5d8-9c2b5612f0d9
let
	pbbr8eb1b2 = scatter(bbr8df.eB1,bbr8df.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" bb0nu (r= 8mm) ",
      label=nothing)
	pxer8eb1b2 = scatter(xer8df.eB1,xer8df.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" Xe137 (r= 8mm) ",
      label=nothing)
	pbbr10eb1b2 = scatter(bbr10df.eB1,bbr10df.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" bb0nu (r= 10 mm) ",
      label=nothing)
	pxer10eb1b2 = scatter(xer10df.eB1,xer10df.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" Xe137 (r= 10mm) ",
      label=nothing)
	

	plot(pbbr8eb1b2,pxer8eb1b2, pbbr10eb1b2,pxer10eb1b2 )
end

# ╔═╡ b8b4f20a-59a6-4346-b85b-33cfc02ac6fe
begin
	bbr10v4mmdf = jn.Petit.merge_csv_files("bb0nu_blobs_4mm_rblob_10")
	xer10v4mmdf = jn.Petit.merge_csv_files("xe137_blobs_4mm_rblob_10") 
end

# ╔═╡ d7409794-6178-4ba4-b364-55ca6261f503
let
	pbbr10v4eb1b2 = scatter(bbr10v4mmdf.eB1,bbr10v4mmdf.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" bb0nu (r= 10mm, v=4mm) ",
      label=nothing)
	pxer10v4eb1b2 = scatter(xer10v4mmdf.eB1,xer10v4mmdf.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" Xe137 (r= 10mm, v=4mm) ",
      label=nothing)
	
	

	plot(pbbr10v4eb1b2,pxer10v4eb1b2)
end

# ╔═╡ 98f2f9c4-29f1-4c7d-9f63-f98f2659c8b9
begin
	bbr10v8mmdf = jn.Petit.merge_csv_files("bb0nu_blobs_8mm_rblob_10")
	xer10v8mmdf = jn.Petit.merge_csv_files("xe137_blobs_8mm_rblob_10") 
end

# ╔═╡ 3c512f54-7896-4ca8-9c58-7002919fffd7
let
	pbbr10v8eb1b2 = scatter(bbr10v8mmdf.eB1,bbr10v8mmdf.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" bb0nu (r= 10mm, v=8mm) ",
      label=nothing)
	pxer10v8eb1b2 = scatter(xer10v8mmdf.eB1,xer10v8mmdf.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" Xe137 (r= 10mm, v=8mm) ",
      label=nothing)
	
	

	plot(pbbr10v8eb1b2,pxer10v8eb1b2)
end

# ╔═╡ 7ec35c71-d5e9-4089-ad6d-71e0237155c6
begin
	bbr15v8mmdf = jn.Petit.merge_csv_files("bb0nu_blobs_8mm_rblob_15")
	xer15v8mmdf = jn.Petit.merge_csv_files("xe137_blobs_8mm_rblob_15") 
	pbbr15v8eb1b2 = scatter(bbr15v8mmdf.eB1,bbr15v8mmdf.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" bb0nu (r= 15mm, v=8mm) ",
      label=nothing)
	pxer15v8eb1b2 = scatter(xer15v8mmdf.eB1,xer15v8mmdf.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" Xe137 (r= 15mm, v=8mm) ",
      label=nothing)
	
	

	plot(pbbr15v8eb1b2,pxer15v8eb1b2)
end

# ╔═╡ 7d33aed5-8731-45a1-910e-fcc7f372354a
begin
	bbr20v8mmdf = jn.Petit.merge_csv_files("bb0nu_blobs_8mm_rblob_20")
	xer20v8mmdf = jn.Petit.merge_csv_files("xe137_blobs_8mm_rblob_20") 
	pbbr20v8eb1b2 = scatter(bbr20v8mmdf.eB1,bbr20v8mmdf.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" bb0nu (r= 20mm, v=8mm) ",
      label=nothing)
	pxer20v8eb1b2 = scatter(xer20v8mmdf.eB1,xer20v8mmdf.eB2, markersize=2, markercolor=:black,
      xlabel="eblob1 (keV)",
      ylabel="eblob2 (keV)",
      title=" Xe137 (r= 20mm, v=8mm) ",
      label=nothing)
	
	
	plot(pbbr20v8eb1b2,pxer20v8eb1b2)
end

# ╔═╡ df0aab75-2c70-40e5-b4ce-5f0088c959a0
md"""
## Functions
"""

# ╔═╡ 4945ede2-bba5-4a6d-818f-482e66acd675
function find_at_working_point(results::Vector{DataFrame}, value::Float64, mode::String)

      output = DataFrame(matched_value=Float64[], result=Float64[], cut=Float64[])

      for df in results
          if mode == "efficiency"
              idx = argmin(abs.(df.eff_signal .- value))
              push!(output, (df.eff_signal[idx], df.bkg_rejection[idx], df.cut[idx]))
          elseif mode == "rejection"
              idx = argmin(abs.(df.bkg_rejection .- value))
              push!(output, (df.bkg_rejection[idx], df.eff_signal[idx], df.cut[idx]))
          else
              error("mode must be 'efficiency' or 'rejection'")
          end
      end

      return output
  end

# ╔═╡ b64344c7-5b7f-41a8-be36-bd24a3d0256f

  function compute_efficiency(df::DataFrame, cut_values::Vector{Float64};
							  column::Symbol=:eB2)
      total = nrow(df)
      total == 0 && return zeros(length(cut_values))

      [count(>(cut), df[!, column]) / total for cut in cut_values]
  end



# ╔═╡ 4dbfcaba-82d6-4c40-ad8b-1df81f6387e3
function plot_roc_curves(df_signals::Vector{DataFrame}, df_bkgs::Vector{DataFrame},
                           labels::Vector{String};
                           column::Symbol=:eB2,
                           cut_min::Float64=0.0,
                           cut_max::Float64=500.0,
                           ncuts::Int=50,
                           title::String="Signal Efficiency vs Background Rejection")

      @assert length(df_signals) == length(df_bkgs) == length(labels) "Vectors must have same length"

      cut_values = collect(range(cut_min, cut_max, length=ncuts))
      results = DataFrame[]

      # Initialize plot
      p = plot(
          xlabel = "Signal Efficiency",
          ylabel = "Background Rejection",
          title = title,
          xlims = (0, 1),
          ylims = (0, 1),
          aspect_ratio = :equal,
          grid = true,
          gridalpha = 0.3,
          legend = :bottomleft
      )

      # Plot each ROC curve
      for (i, (df_sig, df_bkg, label)) in enumerate(zip(df_signals, df_bkgs, labels))
          eff_signal = compute_efficiency(df_sig, cut_values; column=column)
          eff_bkg = compute_efficiency(df_bkg, cut_values; column=column)
          bkg_rejection = 1.0 .- eff_bkg

          # Store results
          push!(results, DataFrame(
              cut = cut_values,
              eff_signal = eff_signal,
              eff_bkg = eff_bkg,
              bkg_rejection = bkg_rejection
          ))

          # Add curve to plot
          plot!(p, eff_signal, bkg_rejection,
              linewidth = 2,
              label = label
          )
      end

      # Add diagonal reference line (random classifier)
      #plot!(p, [0, 1], [1, 0], linestyle=:dash, color=:gray, linewidth=1, label="random")

      return p, results
  end

# ╔═╡ 32ec055d-1963-4cc7-9d22-943e14274914
begin
	df_signals = [bbr3df, bbr5df, bbr8df, bbr10df]
	df_bkgs = [ xer3df,  xer5df,  xer8df,  xer10df]
  	labels = ["r=3mm", "r=5mm", "r=8mm", "r=10mm"]

  proc, rresults = plot_roc_curves(df_signals, df_bkgs, labels;
      column = :eB2,
      cut_min = 0.0,
      cut_max = 700.0,
      ncuts = 300
  )

  plot(proc)
end

# ╔═╡ 6cf21ba6-fcc2-4fa3-9ba6-6c7f6d548cac
begin
	_, rdf = plot_roc_curves(df_signals, df_bkgs, labels;
	      column = :eB2,
	      cut_min = 0.0,
	      cut_max = 700.0,
	      ncuts =30
	  )
	rocr3mm = rdf[1]
	rocr5mm = rdf[2]
	rocr8mm = rdf[3]
	rocr10mm = rdf[4]
	rocr3mm
end

# ╔═╡ 25d558c5-7c67-4290-bd91-c32e5def857d
begin
	df_signals2 = [bbr10v4mmdf]
	df_bkgs2 = [xer10v4mmdf]
  	labels2 = ["r=10mm: v = 4mm"]

  proc2, rresults2 = plot_roc_curves(df_signals2, df_bkgs2, labels2;
      column = :eB2,
      cut_min = 0.0,
      cut_max = 700.0,
      ncuts = 300
  )
	rocv4mmr10mm =rresults2[1]
  plot(proc2)
end

# ╔═╡ ba4535c5-0cac-4565-b3dd-9159d4d2632a
begin
	df_signals3 = [bbr10df, bbr10v4mmdf]
	df_bkgs3 = [xer10df, xer10v4mmdf]
  	labels3 = ["v 1mm", "v 4mm"]

  proc3, rresults3 = plot_roc_curves(df_signals3, df_bkgs3, labels3;
      column = :eB2,
      cut_min = 0.0,
      cut_max = 700.0,
      ncuts = 300
  )

  plot(proc3)
end

# ╔═╡ 90ebbadc-bccf-4d46-bd5f-c68475b6f2ff
begin
	df_signals4 = [bbr10df, bbr10v4mmdf, bbr15v8mmdf, bbr20v8mmdf]
	df_bkgs4 = [xer10df, xer10v4mmdf, xer15v8mmdf, xer20v8mmdf]
  	labels4 = ["V1mmR10", "V4mmR10", "V8mmR15", "V8mmR20"]

  proc4, rresults4 = plot_roc_curves(df_signals4, df_bkgs4, labels4;
      column = :eB2,
      cut_min = 0.0,
      cut_max = 1000.0,
      ncuts = 400
  )
	rocv8mmr20mm = rresults4[4]
  plot(proc4)
end

# ╔═╡ c4287e08-6881-4be6-af0e-9b88a46844ad
# Get background rejection at 80% signal efficiency for each dataset
  bkg_rejections = find_at_working_point([rocr10mm,rocv4mmr10mm,rocv8mmr20mm], 0.8, "efficiency")

 

# ╔═╡ 584072ea-623f-4400-be26-12429825a049
 # Get signal efficiency at 90% background rejection for each dataset
  efficiencies = find_at_working_point([rocr10mm,rocv4mmr10mm, rocv8mmr20mm], 0.9, "rejection")

# ╔═╡ f4ad81a7-b6c4-43a0-9ccf-88b6171593ef
 
function plot_signal_eff_vs_bkg_rejection(df_signal::DataFrame, df_bkg::DataFrame;
                                             column::Symbol=:eB2,
                                             cut_min::Float64=0.0,
                                             cut_max::Float64=500.0,
                                             ncuts::Int=50)

      cut_values = collect(range(cut_min, cut_max, length=ncuts))

      eff_signal = compute_efficiency(df_signal, cut_values; column=column)
      eff_bkg = compute_efficiency(df_bkg, cut_values; column=column)
      bkg_rejection = 1.0 .- eff_bkg

      # Create results DataFrame
      results = DataFrame(
          cut = cut_values,
          eff_signal = eff_signal,
          eff_bkg = eff_bkg,
          bkg_rejection = bkg_rejection
      )

      # Plot signal efficiency vs background rejection
      p = plot(eff_signal, bkg_rejection,
          xlabel = "Signal Efficiency",
          ylabel = "Background Rejection",
          title = "Signal Efficiency vs Background Rejection",
          label = nothing,
          linewidth = 2,
          marker = :circle,
          markersize = 3,
          xlims = (0, 1),
          ylims = (0, 1),
          aspect_ratio = :equal,
          grid = true,
          gridalpha = 0.3
      )

      # Add diagonal reference line (random classifier)
      #plot!(p, [0, 1], [1, 0], linestyle=:dash, color=:gray, label="random")

      return p, results
  end

 

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
# ╠═e91c9f90-63a6-48cf-ab8a-3d8d4ac320c4
# ╠═c5f2e0aa-a7ee-406c-be13-e6219c409638
# ╠═84eb8703-228b-46f2-8868-5086f576572e
# ╠═0289428e-6db3-4bfc-82fe-5b9e71dd3706
# ╠═726badaf-c389-484e-a46d-cdda4d4d578d
# ╠═67ac7cde-e9a4-4013-8e30-61ac29d28175
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
# ╠═105bd81d-963d-4ac1-9bfc-1f3179f90202
# ╠═7364c781-7714-46f8-b412-2dd04ee4b8a5
# ╠═caa5c1d8-5812-42f4-892a-34c14599b3e0
# ╠═d089e072-9ae9-4737-9e1e-44d84e354d28
# ╠═cf81a428-8cb7-4c19-8e58-f87778e867fd
# ╠═085ebb25-027a-4027-9f4c-a69013aee707
# ╠═e56d5fa6-02d9-4260-aca5-812164409a48
# ╠═e12c925b-101b-4829-bb58-e98debb454a9
# ╠═4067b243-9679-4402-8570-5f351b489805
# ╠═985ff432-1d22-4e22-86e1-6382eaee243b
# ╠═d8e42ebf-e7ec-42fb-9867-3af552ffb842
# ╠═b4ce8dcb-52fd-462a-b864-4e2b11a11990
# ╠═975b6fa8-fac9-48ea-9e31-cf7d397ebf1c
# ╠═dd011045-aca5-49a9-8227-aa66097f98ca
# ╠═d656e40b-c9e7-4be2-91eb-8f56fbca9495
# ╠═7db0dde5-2760-47aa-a3c4-e5fb93e3f36b
# ╠═ef541827-3bcc-497d-8a1c-9a2d05d5e75a
# ╠═d37787c8-8511-40b6-aee5-8064cf5ed268
# ╠═51575ff2-a24a-4177-90fe-c904a20c2d57
# ╠═c79c2d21-a22d-4ed2-b7a7-3a5f7a6b9b4c
# ╠═e9b013a3-3353-4035-8830-5737d952bf52
# ╠═79e0190f-054a-4455-9df8-7ba381eb6f89
# ╠═4995d151-11e6-46b0-a19d-c61ca0de8b8c
# ╠═9a89e70a-f40b-4fda-a17d-8a875a8c327a
# ╠═ccfa2e3c-92d8-4108-a21e-4f09bb1d4878
# ╠═66fc61a5-3fa0-43aa-a005-d6d0abf7e491
# ╠═1e9da122-9d8a-4b9e-a5df-fac798a5f8a8
# ╠═84208e5e-6209-427b-a5d8-9c2b5612f0d9
# ╠═32ec055d-1963-4cc7-9d22-943e14274914
# ╠═6cf21ba6-fcc2-4fa3-9ba6-6c7f6d548cac
# ╠═b8b4f20a-59a6-4346-b85b-33cfc02ac6fe
# ╠═d7409794-6178-4ba4-b364-55ca6261f503
# ╠═25d558c5-7c67-4290-bd91-c32e5def857d
# ╠═ba4535c5-0cac-4565-b3dd-9159d4d2632a
# ╠═98f2f9c4-29f1-4c7d-9f63-f98f2659c8b9
# ╠═3c512f54-7896-4ca8-9c58-7002919fffd7
# ╠═7ec35c71-d5e9-4089-ad6d-71e0237155c6
# ╠═7d33aed5-8731-45a1-910e-fcc7f372354a
# ╠═90ebbadc-bccf-4d46-bd5f-c68475b6f2ff
# ╠═c4287e08-6881-4be6-af0e-9b88a46844ad
# ╠═584072ea-623f-4400-be26-12429825a049
# ╠═df0aab75-2c70-40e5-b4ce-5f0088c959a0
# ╠═4945ede2-bba5-4a6d-818f-482e66acd675
# ╠═b64344c7-5b7f-41a8-be36-bd24a3d0256f
# ╠═4dbfcaba-82d6-4c40-ad8b-1df81f6387e3
# ╠═f4ad81a7-b6c4-43a0-9ccf-88b6171593ef
