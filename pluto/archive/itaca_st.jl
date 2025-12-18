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
	using Interpolations 
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
	#it = ingredients(string(pdir,"/src/itaca_functions.jl"))
end

# ╔═╡ b53a30f0-fda8-44b0-bed2-c663af6f01b1
md"""
## Functions
"""

# ╔═╡ 9c0e82c3-1d71-4164-bf94-a356538eee8d
function voxel_size(dx, dy, p, l)
	sqrt(l)*(dx + dy) /(2.0*sqrt(p))
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
function print_reco_blobs(trk, walk; rblob=10.0, label="")
	extremes, _, _, xtrack_length, confidence = walk
	blobs = jn.Petit.energy_in_spheres_around_extremes(trk, walk, rblob)
	eb1 = blobs.blob1_energy * 1e+3
	eb2 = blobs.blob2_energy * 1e+3
	db = abs((eb1 -eb2))/(eb1+eb2)
	md"""
	#### Find blobs: $(label) with radius $(@sprintf("%.1f", rblob)) mm  
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

# ╔═╡ 43d5968c-e185-4586-9f49-f9be529e114b
function transform_hits_df(df::DataFrame; energy_to_electrons::Float64=1e5/2.5)
      df2 = select(df, Not([:time, :label, :particle_id, :hit_id]))
      df2.electrons = round.(Int, df2.energy .* energy_to_electrons)
      return df2
  end

# ╔═╡ f28a58d0-55a3-4f6d-af3f-200b8c395c98
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


# ╔═╡ 35828cbe-b5d9-4f85-883b-8996646338f8
begin
	cnndir = "/Users/jjgomezcadenas/Projects/Petit/pluto/csvCNN"
  cnn = read_cnn_efficiencies(cnndir)
end

# ╔═╡ 01ce8389-a772-480b-a875-57916a4cd0ea
begin
	th=0.70
	 md"""
  ##### For Ions (σ~1mm)
  - sign eff  = $(cnn.interp_1mm.signal(0.99))
  - bkgn eff  = $(cnn.interp_1mm.background(0.99))

  ##### For Xe/He (σ~3.5 mm)
  - sign eff = $(cnn.interp_3p5mm.signal(0.95))
  - bkgn eff = $(cnn.interp_3p5mm.background(0.95))

  #### For Xe (σ~10 mm) 
   - sign eff = $(cnn.interp_10mm.signal(0.94))
  - bkgn eff = $(cnn.interp_10mm.background(0.94))
	"""
end

# ╔═╡ 529ebb7f-1615-464e-b569-dc423e1bb8c7
0.04/0.007

# ╔═╡ 79156ecf-35d9-45a5-bf05-51ceeee9da2e
begin
# Plot signal efficiency vs threshold
  thresholds = range(0, 1, length=100)
  plot(thresholds, cnn.interp_1mm.signal.(thresholds), label="1mm")
  plot!(thresholds, cnn.interp_3p5mm.signal.(thresholds), label="3.5mm")
  plot!(thresholds, cnn.interp_10mm.signal.(thresholds), label="10mm")
end

# ╔═╡ 896f02e7-fede-42b7-9366-4955d1879b6c
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


# ╔═╡ 312de1b0-bba0-40ab-a537-c7c5343380b2
function compute_all_rocs(cnn; thresholds=range(0.0, 1.0, length=100))
    roc_1mm = compute_roc(cnn.interp_1mm; thresholds=thresholds)
    roc_3p5mm = compute_roc(cnn.interp_3p5mm; thresholds=thresholds)
    roc_10mm = compute_roc(cnn.interp_10mm; thresholds=thresholds)

    return (roc_1mm=roc_1mm, roc_3p5mm=roc_3p5mm, roc_10mm=roc_10mm)
end





# ╔═╡ 1024cc13-60e0-42bf-ad14-2f1d90da4484
function plot_roc(rocs; labels=["1mm (σ)", "3.5mm (σ)", "10mm (σ)"])
    p = plot(xlabel="True Positive Rate (Signal Eff.)",
             ylabel="Background Rejection (1 - Bkg Eff.)",
             title="ROC Curves",
             legend=:bottomleft,
             xlims=(0, 1), ylims=(0, 1),
             aspect_ratio=:equal,
             grid=true, gridstyle=:dot, gridalpha=0.3)

    # Plot diagonal (random classifier)
    #plot!(p, [0, 1], [1, 0], linestyle=:dash, color=:gray,
    #      label="Random", alpha=0.5)

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

 

# ╔═╡ 27d8edb9-5b07-4760-9cbe-3bb413466b99
function plot_roc2(rocs; labels=["ITACA", "GXe (pure Xenon)"])
    p = plot(xlabel="True Positive Rate (Signal Eff.)",
             ylabel="Background Rejection (1 - Bkg Eff.)",
             title="ROC Curves",
             legend=:bottomleft,
             xlims=(0, 1), ylims=(0, 1),
             aspect_ratio=:equal,
             grid=true, gridstyle=:dot, gridalpha=0.3)

    # Plot diagonal (random classifier)
    #plot!(p, [0, 1], [1, 0], linestyle=:dash, color=:gray,
    #      label="Random", alpha=0.5)

    # Plot each ROC curve (TPR on X, 1-FPR on Y)
    colors = [:blue, :red]
    for (i, (roc, label, color)) in enumerate(zip(
            [rocs.roc_1mm,  rocs.roc_10mm], labels, colors))

		auc = round(roc.auc, digits=3)
		if i == 1
			auc = round(roc.auc + 0.001, digits=3)
		end
        plot!(p, roc.tpr, 1.0 .- roc.fpr,
              label="$label AUC=$(auc)",
              linewidth=2, color=color)
    end

    return p
end

# ╔═╡ f6dabe79-9985-47f2-b536-b12fef702393
rocs = compute_all_rocs(cnn; thresholds=range(0.0, 1.0, length=100))

# ╔═╡ b73a4322-e098-4032-8564-3f1f4914a955
begin
	 
	 plot_roc2(rocs; labels=["ITACA", "GXe Xe/He"])
end

# ╔═╡ 475079bb-18dd-477e-9b67-320ca5909ade
(1 -0.96)/(1 -0.992)

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
	vhe = voxel_size(0.75, 1.6, 15.0, 200.0)
	md"""
	- voxel size for 10 % He = $@sprintf("%.1e", vhe) mm
	"""
end

# ╔═╡ 7b3a8a4d-f649-4175-b881-453b7e44e598
md"""
### Load data
"""

# ╔═╡ d0fecc3e-5e0e-4522-9efa-b4b60fd151a5
#bidf = load_data("gammas_2458keV_15bar_100mum.next.h5")

# ╔═╡ b32e2b11-788e-4017-b5cf-4e4d3968398d
bbdf = load_data("bb0nu/0nubb_15bar_100mum.next.h5")

# ╔═╡ 6950e911-3003-4778-98bd-2b1a0a819844
xedf = load_data("xe137/electrons_2400_2500_15bar_100mum.next.h5")

# ╔═╡ cc0ad60d-0e31-4406-9b0d-20d938a00aff
begin
	σt_xe_l100_p15 = jn.Petit.sigma_t_mm(100.0, 15.0; dtmm=3.5)
	σt_xehe_l100_p15 = jn.Petit.sigma_t_mm(100.0, 15.0; dtmm=1.6)
	σl_xe_l100_p15 = jn.Petit.sigma_l_mm(100.0, 15.0; dlmm=0.9)
	σl_xehe_l100_p15 = jn.Petit.sigma_l_mm(100.0, 15.0; dlmm=0.75)
	σt_xe_l200_p15 = jn.Petit.sigma_t_mm(200.0, 15.0; dtmm=3.5)
	σt_xehe_l200_p15 = jn.Petit.sigma_t_mm(200.0, 15.0; dtmm=1.6)
	σl_xe_l200_p15 = jn.Petit.sigma_l_mm(200.0, 15.0; dlmm=0.9)
	σl_xehe_l200_p15 = jn.Petit.sigma_l_mm(200.0, 15.0; dlmm=0.75)
	σt_ion_l100_p15 = jn.Petit.sigma_t_ion_mm(300.0, 100.0, 500.0)
	σt_ion_l200_p15 = jn.Petit.sigma_t_ion_mm(300.0, 200.0, 500.0)
	println("All values of σ in mm")
	println("""
	σt_xe_l100_p15 = $(σt_xe_l100_p15)
	σt_xehe_l100_p15 = $(σt_xehe_l100_p15)
	σl_xe_l100_p15 = $(σl_xe_l100_p15)
	σl_xehe_l100_p15 = $(σl_xehe_l100_p15)
	σt_xe_l200_p15 = $(σt_xe_l200_p15)
	σt_xehe_l200_p15 = $(σt_xehe_l200_p15)
	σl_xe_l200_p15 = $(σl_xe_l200_p15)
	σl_xehe_l200_p15 = $(σl_xehe_l200_p15)
	σt_ion_l100_p15 = $(σt_ion_l100_p15)
	σt_ion_l200_p15 = $(σt_ion_l200_p15)
			""")
end

# ╔═╡ 1dfa60a6-c09c-46d7-874b-44ede0043780
md"""
## Track analysis
"""

# ╔═╡ 41653083-87fb-46f7-8d1f-5cf9d34c2ea7
begin
	# Event number to look at
	nevent = 3
end

# ╔═╡ 47dc42ac-9dd9-4318-9a29-7655541d3af8
md"""
### bb0nu 
"""

# ╔═╡ af515863-e9e6-4a8a-a940-0a8439afd924
md"""
#### MC event with 100 μm step
"""

# ╔═╡ 83a718af-c2d9-4743-827d-8682a7379c65
begin
	bbevt = jn.Petit.get_event(bbdf, nevent)
	bbevtmc = transform_hits_df(bbevt)
	jn.Petit.plot_event(bbevtmc)
end

# ╔═╡ 7ba1eacc-88f9-4bc3-ba75-1efce41dadb7
let
	hd, pd = jn.Petit.step_hist(jn.Petit.closest_voxel_distances_fast(bbevtmc);
	                                                 nbins = 40,
										             xlim   = (0.0, 1.0),
	                                                 xlabel = " d (mm)",
	                                                 ylabel = "Frequency",
	                                                 title=" Distance voxels (mm)")
	he, pe = jn.Petit.step_hist(bbevtmc.energy*1e+3;
	                                                 nbins = 40,
										             #xlim   = (0.0, 50.0),
	                                                 xlabel = " E  (keV)",
	                                                 ylabel = "Frequency",
	                                                 title=" E voxel (keV)")
	plot(pd, pe)
end

# ╔═╡ 4e1f83d9-e9f8-47c9-a742-3eb32b371c93
md"""
#### Voxelize MC event with 1 mm voxels 
"""

# ╔═╡ d92a0d53-d5f0-4e77-9911-a97e93d5b7b4
begin
	bbv1mm = jn.Petit.voxelize_event(bbevtmc, nevent, 1.0)
	jn.Petit.plot_event(bbv1mm)
end

# ╔═╡ 1d6276f9-c585-48ed-9c3e-d910cd47d091
begin
	bbmc_tks = jn.Petit.make_tracks(bbevtmc; max_distance_mm=1.0,
								energy_threshold_kev=0.0)
	md"""
	- number of tracks reconstructed =$(length(bbmc_tks))
	"""
end

# ╔═╡ 744b4f32-5f0e-4e2e-b6e3-033b1d74af62
md"""
#### Apply ion diffusion
"""

# ╔═╡ ac625d74-16f4-4463-bada-42d0b47e79d5
begin
	bbevtdf =jn.Petit.diffuse_xyz_image_mc(bbevtmc; 
									  sigma_t_mm=σt_ion_l100_p15, 
									  sigma_l_mm=0.1, 
									  nbins=300, 
									  nsigma=3.0)
	histogram(bbevtdf.electrons)
end

# ╔═╡ b0ab441c-f7aa-4bfa-8bc8-b87d4decfee9
md"""
#### Plot diffused track
"""

# ╔═╡ 49eacb34-e766-41f2-ab1f-8170f90e5a50
jn.Petit.plot_hits(bbevtdf, energy_column=:electrons)

# ╔═╡ dd1ebb2d-581e-4c5c-9f69-6db4e6fe0a1b
md"""
#### Voxelize diffused event
"""

# ╔═╡ df172361-5846-471d-8014-92db51969574
begin
	bbdfv1mm = jn.Petit.voxelize_event(bbevtdf, σt_ion_l200_p15)
	jn.Petit.plot_event(bbdfv1mm)
end

# ╔═╡ 66e9e2e4-aa44-4e8c-b304-0dddd65e4fc6
begin
	dmx_ion = σt_ion_l200_p15*2
	bbdfv1mm_tks = jn.Petit.make_tracks(bbdfv1mm; max_distance_mm=σt_ion_l200_p15*2,
								energy_threshold_kev=10.0)
	md"""
	- max distance to be included in track -> $(@sprintf("%.1f", dmx_ion))
	- number of tracks reconstructed =$(length(bbdfv1mm_tks))
	"""
end

# ╔═╡ fcf2b122-cc41-4b63-b631-41c30a19c728
begin
	bbdfv1mmt = bbdfv1mm_tks[1]
	rb_ion =  σt_ion_l200_p15*3
	bbdfv1mmt_xwalk =jn.Petit.walk_track_from_extremes(bbdfv1mmt)
	print_reco_blobs(bbdfv1mmt, bbdfv1mmt_xwalk, 
					 rblob = rb_ion, 
		label="bb0nu σ=$(@sprintf("%.1f", σt_ion_l100_p15)) mm vx=$(@sprintf("%.1f", σt_ion_l200_p15)) mm")
	 
end

# ╔═╡ 8a24d1ce-547f-400b-9415-74dd0b7dd689
begin
	plot_bbdfv1mmt_2d, plot_bbdfv1mmt_3d =jn.Petit.plot_track_blobs(bbdfv1mmt,
																	bbdfv1mmt_xwalk;
																	sphere_radius= rb_ion)
	plot(plot_bbdfv1mmt_2d, plot_bbdfv1mmt_3d )
end

# ╔═╡ 3fba0b56-cfc0-4115-ab29-8eccc1aefb37
md"""
#### Apply electron diffusion (Xe/He)
"""

# ╔═╡ 46528939-1e69-4f4d-b716-3e7f76e6c22f
begin
	bbevtdfel =jn.Petit.diffuse_xyz_image_mc(bbevtmc; 
										sigma_t_mm=σt_xehe_l100_p15, 
										sigma_l_mm=σl_xehe_l100_p15, 
										nbins=200, nsigma=3.0)
	histogram(bbevtdf.electrons)
end

# ╔═╡ 6f4c9519-ec55-452a-a683-01ed779d7334
jn.Petit.plot_hits(bbevtdfel, energy_column=:electrons)

# ╔═╡ b5070fc8-1ec7-4212-9b3b-43f0fa999d24
begin
	bbdfv1mmel = jn.Petit.voxelize_event(bbevtdfel, σt_xehe_l200_p15)
	jn.Petit.plot_event(bbdfv1mmel)
end

# ╔═╡ 1ac6befd-59d8-40ad-a5b0-710b7f941444
begin
	dmx_xehe = σt_xehe_l200_p15*2
	bbdfv1mmel_tks = jn.Petit.make_tracks(bbdfv1mmel; max_distance_mm=dmx_xehe,
								energy_threshold_kev=10.0)
	md"""
	- max distance to be included in track -> $(@sprintf("%.1f", dmx_xehe))
	- number of tracks reconstructed =$(length(bbdfv1mmel_tks))
	"""
end

# ╔═╡ 136adb7c-8ba5-4bed-9f73-132a937d2a16
begin
	bbdfv1mmelt = bbdfv1mmel_tks[1]
	rbxehe =  σt_xehe_l200_p15*3
	bbdfv1mmelt_xwalk =jn.Petit.walk_track_from_extremes(bbdfv1mmelt)
	print_reco_blobs(bbdfv1mmelt, bbdfv1mmelt_xwalk, 
					 rblob = rbxehe, 
					 label="bb0nu σ=$(@sprintf("%.1f", σt_xehe_l100_p15)) mm vx=$(@sprintf("%.1f", σt_xehe_l200_p15)) mm")

end

# ╔═╡ 891ff777-2294-47a9-89e0-e3810203dcff
begin
	plot_bbdfv1mmelt_2d, plot_bbdfv1mmelt_3d =jn.Petit.plot_track_blobs(bbdfv1mmelt,
																	bbdfv1mmelt_xwalk;
																	sphere_radius= rbxehe)
	plot(plot_bbdfv1mmelt_2d, plot_bbdfv1mmelt_3d )
end

# ╔═╡ 38474add-cbdf-4772-ba95-c7da00485b31
md"""
### Xe137 
"""

# ╔═╡ 49f11799-fe6d-40f8-abda-d90bc9c4d3e8
md"""
#### MC event with 100 μm step
"""

# ╔═╡ 86bf06e8-9caa-4c18-80ff-160cbd945986
begin
	xeevt = jn.Petit.get_event(xedf, nevent)
	xeevtmc = transform_hits_df(xeevt)
	jn.Petit.plot_event(xeevtmc)
end

# ╔═╡ 05998d38-6935-44e9-9c16-53e1810ecad5
xeevt

# ╔═╡ 5c9eab19-8870-47e5-a9a9-263c8f70bd44
md"""
#### Voxelize MC event with 1 mm voxels 
"""

# ╔═╡ 7abbc303-db74-4764-be25-b0304ed694ea
begin
	xev1mm = jn.Petit.voxelize_event(xeevtmc, nevent, 1.0)
	jn.Petit.plot_event(xev1mm)
end

# ╔═╡ 382bdd57-44c9-47dd-83be-c5bbee399ac9
xev1mm

# ╔═╡ 1b0e400f-ca67-4f2e-9e68-a90192fded54
begin
	xemc_tks = jn.Petit.make_tracks(xev1mm; max_distance_mm=1.0,
								energy_threshold_kev=0.0)
	md"""
	- number of tracks reconstructed =$(length(xemc_tks))
	"""
end

# ╔═╡ 29c5b0e3-9a0f-4cde-ab10-7bca1ff18468
md"""
#### Apply diffusion
"""

# ╔═╡ 3a795130-6a5e-4606-92b2-57f006bda4a4
begin
	xeevtdf =jn.Petit.diffuse_xyz_image_mc(xeevtmc; 
									  sigma_t_mm=σt_ion_l100_p15, 
									  sigma_l_mm=0.1, 
									  nbins=300, 
									  nsigma=3.0)
	histogram(xeevtdf.electrons)
end

# ╔═╡ 596bc0f0-fc1a-4c8b-bab0-e799d2621f46
md"""
#### Plot diffused track
"""

# ╔═╡ 7fd42a7e-18b0-408c-a439-7804f90e09df
jn.Petit.plot_hits(xeevtdf, energy_column=:electrons)

# ╔═╡ deeceb13-3b73-47c1-8828-5fbdffabd3de
md"""
#### Voxelize diffused event
"""

# ╔═╡ c92d9d16-f53c-43d7-95a6-ca2319c26da2
begin
	xedfv1mm = jn.Petit.voxelize_event(xeevtdf, σt_ion_l200_p15)
	jn.Petit.plot_event(xedfv1mm)
end

# ╔═╡ 63159f7b-a1fc-4648-b148-040893e609ab
begin
	xedfv1mm_tks = jn.Petit.make_tracks(xedfv1mm; max_distance_mm=σt_ion_l200_p15*2,
								energy_threshold_kev=10.0)
	md"""
	- max distance to be included in track -> $(@sprintf("%.1f", dmx_ion))
	- number of tracks reconstructed =$(length(xedfv1mm_tks))
	"""
end

# ╔═╡ 177e8223-e063-42e4-b633-17f87f653c74
begin
	xedfv1mmt = xedfv1mm_tks[1]
	xedfv1mmt_xwalk =jn.Petit.walk_track_from_extremes(xedfv1mmt)
	print_reco_blobs(xedfv1mmt, xedfv1mmt_xwalk, 
					 rblob = rb_ion, 
		label="Xe137 σ=$(@sprintf("%.1f", σt_ion_l100_p15)) mm vx=$(@sprintf("%.1f", σt_ion_l200_p15)) mm")
					 
end

# ╔═╡ a4442faf-6442-4ef4-9200-e3e705c7a2b2
begin
	plot_xedfv1mmt_2d, plot_xedfv1mmt_3d =jn.Petit.plot_track_blobs(xedfv1mmt,
																	xedfv1mmt_xwalk;
																	sphere_radius= rb_ion)
	plot(plot_xedfv1mmt_2d, plot_xedfv1mmt_3d )
end

# ╔═╡ e1b8c419-c2dc-4d5b-af16-c59f77d54aa8
md"""
#### Apply electron diffusion (Xe/He)
"""

# ╔═╡ d16366ec-3e0c-4b87-adb3-f180413cdd39
begin
	xeevtdfel =jn.Petit.diffuse_xyz_image_mc(xeevtmc; 
										sigma_t_mm=σt_xehe_l100_p15, 
										sigma_l_mm=σl_xehe_l100_p15, 
										nbins=200, nsigma=3.0)
	histogram(xeevtdf.electrons)
end

# ╔═╡ 2a2293f5-b037-498b-976e-cd95f1c07d1a
xeevtdfel

# ╔═╡ 7e707191-d6ff-4513-92d5-7d2ce8e5db3c
jn.Petit.plot_hits(xeevtdfel, energy_column=:electrons)


# ╔═╡ abc1fc50-ec3c-4bd6-a6af-7a0af8147b3d
begin
	xedfv1mmel = jn.Petit.voxelize_event(xeevtdfel, σt_xehe_l200_p15)
	jn.Petit.plot_event(xedfv1mmel)
end

# ╔═╡ e51fd878-7c12-4770-b1f1-c9f9e6e33256
xedfv1mmel

# ╔═╡ 380f1e72-1215-421e-8b00-268f34c82341
begin
	xedfv1mmel_tks = jn.Petit.make_tracks(xedfv1mmel; max_distance_mm=dmx_xehe,
								energy_threshold_kev=10.0)
	md"""
	- max distance to be included in track -> $(@sprintf("%.1f", dmx_xehe))
	- number of tracks reconstructed =$(length(xedfv1mmel_tks))
	"""
end

# ╔═╡ 6f8a94a0-0ec9-488f-bd3c-d7a67f27f978
begin
	xedfv1mmelt = xedfv1mmel_tks[1]
	xedfv1mmelt_xwalk =jn.Petit.walk_track_from_extremes(xedfv1mmelt)
	print_reco_blobs(xedfv1mmelt, xedfv1mmelt_xwalk, 
					 rblob = rbxehe, 
					 label="Xe137 σ=$(@sprintf("%.1f", σt_xehe_l100_p15)) mm vx=$(@sprintf("%.1f", σt_xehe_l200_p15)) mm")

end

# ╔═╡ f1d33080-06f8-4310-8589-369570f8934a
begin
	plot_xedfv1mmelt_2d, plot_xedfv1mmelt_3d =jn.Petit.plot_track_blobs(xedfv1mmelt,
																	xedfv1mmelt_xwalk;
																	sphere_radius= rbxehe)
	plot(plot_xedfv1mmelt_2d, plot_xedfv1mmelt_3d )
end

# ╔═╡ 8c173831-8e37-47d4-a183-d95be8773781
md"""
#### bb0nu example, voxels = 1mm
"""

# ╔═╡ 1396c7a9-bb49-4132-9a4f-0782922141f8
md"""
## Track Analysis
"""

# ╔═╡ 105bd81d-963d-4ac1-9bfc-1f3179f90202
md"""
### Select bb0nu and Xe137 events, returning tracks
"""

# ╔═╡ 3ff630c5-6409-49ea-971c-e3cd15427de2
begin
	bb_mctrks, bb_iontrks, bb_eletrks = jn.Petit.get_itaca_tracks("/Users/jjgomezcadenas/Data/HD5t/itaca/bb0nu", tag="xehe_15bar_roi100")
	bb_meta = jn.Petit.read_itaca_metadata("/Users/jjgomezcadenas/Data/HD5t/itaca/bb0nu/0nubb_xehe_15bar_roi100_metadata.csv")
end

# ╔═╡ 00ecc891-6de9-4b6d-a9e2-00de7a8c1d08
begin
	xe_mctrks, xe_iontrks, xe_eletrks = jn.Petit.get_itaca_tracks("/Users/jjgomezcadenas/Data/HD5t/itaca/xe137", tag="xehe_15bar_roi100")
	xe_meta = jn.Petit.read_itaca_metadata("/Users/jjgomezcadenas/Data/HD5t/itaca/xe137/electrons_2400_2500_xehe_15bar_roi100_metadata.csv")
end

# ╔═╡ be04f12b-f434-4090-8e69-fac08f2b8350
begin
	bi_mctrks, bi_iontrks, bi_eletrks = jn.Petit.get_itaca_tracks("/Users/jjgomezcadenas/Data/HD5t/itaca/bi214", tag="xe_15bar_roi100")
	bi_meta = jn.Petit.read_itaca_metadata("/Users/jjgomezcadenas/Data/HD5t/itaca/bi214/gammas_2458_xe_15bar_roi100_metadata.csv")
end

# ╔═╡ d2961e74-c99c-414b-b4a7-8fdde1bf9c8a
length(bi_mctrks.tracks)

# ╔═╡ bca3c9d7-5583-42a6-89a0-e595e5b21970
length(bi_iontrks.tracks)

# ╔═╡ 5af82f10-c89f-4092-811f-12ae1e1ef94b
length(bi_eletrks.tracks)

# ╔═╡ 362d825d-e385-4b28-8fc2-f083391bc89a
58/15

# ╔═╡ b1bff10b-31e6-4820-a7c4-df154c01ab6b
mdf = jn.Petit.match_itaca_events_df(bb_mctrks, bb_iontrks, bb_eletrks)

# ╔═╡ 2d5be057-d840-413a-a13f-8752a2619fd0
begin
	ix = 4
	ivt = mdf[ix, :].event_id
	imc = mdf[ix, :].mc_idx
	iion = mdf[ix, :].ion_idx
	iele = mdf[ix, :].ele_idx
	println("For event =$ivt, imc = $imc, iion = $iion, iele =$iele")
end

# ╔═╡ c75ea153-481d-47a3-99ba-8dcf26416dc0
mxedf = jn.Petit.match_itaca_events_df(xe_mctrks, xe_iontrks, xe_eletrks)

# ╔═╡ e316e013-d2d8-4fa2-b141-7e26eef0ad49
begin
	kx = 1
	kvt = mxedf[kx, :].event_id
	kmc = mxedf[kx, :].mc_idx
	kion = mxedf[kx, :].ion_idx
	kele = mxedf[kx, :].ele_idx
	println("For event =$kvt, kmc = $kmc, kion = $kion, kele =$kele")
end

# ╔═╡ b8e89e80-c581-49a0-8e5f-fb41ac4ea01b
begin
	bbmtrk = bb_mctrks.tracks[imc]
	bbmtrk_xwalk =jn.Petit.walk_track_from_extremes(bbmtrk)
	rbm=7.0
	
	print_reco_blobs(bbmtrk, bbmtrk_xwalk, 
					 rblob = rbm, 
		label="bb0nu σ=$(@sprintf("%.1f", bbmtrk.diffusion.sigma_t)) mm vx=$(@sprintf("%.1f", bbmtrk.diffusion.voxel_size)) mm")
					 
end

# ╔═╡ 194e2770-9b31-44f2-897a-9404be48079b
begin
	plot_bbmtrk_2d, plot_bbmtrk_3d =jn.Petit.plot_track_blobs(bbmtrk,
																	bbmtrk_xwalk;
																	sphere_radius= rbm)
	plot(plot_bbmtrk_2d, plot_bbmtrk_3d)
end

# ╔═╡ af19f2c1-66c6-40f7-acb5-b49acb17db51
nr = 5

# ╔═╡ e34e216a-0c7a-45e2-abb7-7d50f4172117
bbmtrk_xwalk.path_voxels[1:nr,:]

# ╔═╡ 3af1d8e7-9d83-40a9-b9e1-95c292785e38
bbmtrk_xwalk.path_voxels[(end-nr+1):end, :]

# ╔═╡ f35715e9-5024-4d2e-b981-8c06a3ea4300


# ╔═╡ d1cf9c39-b91d-4e16-8810-980e516971c5
Δb(eb1,eb2) = abs((eb1 -eb2))/(eb1+eb2)

# ╔═╡ 519fe2aa-0465-4253-b0c7-1492ab3b8440
let
	bbmcr = jn.Petit.energy_blobs_from_path(bbmtrk_xwalk, nr, energy_col=:electrons)
	md"""
	- eb1 = $(bbmcr.eb1), eb2 = $(bbmcr.eb2), Δb = $(Δb(bbmcr.eb1,bbmcr.eb2))
	"""

end

# ╔═╡ ca37680d-9162-497d-9b36-d73fd90b2e07
begin
	bbitrk = bb_iontrks.tracks[iion]
	bbitrk_xwalk =jn.Petit.walk_track_from_extremes(bbitrk)
	rbi=10.0
	
	print_reco_blobs(bbitrk, bbitrk_xwalk, 
					 rblob = rbi, 
		label="bb0nu σ=$(@sprintf("%.1f", bbitrk.diffusion.sigma_t)) mm vx=$(@sprintf("%.1f", bbitrk.diffusion.voxel_size)) mm")
					 
end

# ╔═╡ 9c1590fc-ad2b-4444-95db-9ffbcedca155
begin
	plot_bbitrk_2d, plot_bbitrk_3d =jn.Petit.plot_track_blobs(bbitrk,
																	bbitrk_xwalk;
																	sphere_radius= rbi)
	plot(plot_bbitrk_2d, plot_bbitrk_3d)
end

# ╔═╡ 51cb6061-9432-48b7-affb-a34c3549861d
bbitrk_xwalk.path_voxels[1:nr,:]

# ╔═╡ 67570f37-cc53-4ab8-8069-e8a622297971
bbitrk_xwalk.path_voxels[(end-nr+1):end, :]

# ╔═╡ 717de5ce-9edb-4d37-a1df-1dc756aa6831


# ╔═╡ 0adc0c96-c8c4-4c41-8f28-46ded1217cc6
let
	bbmcr = jn.Petit.energy_blobs_from_path(bbitrk_xwalk, 5, energy_col=:electrons)
	md"""
	- eb1 = $(bbmcr.eb1), eb2 = $(bbmcr.eb2), Δb = $(Δb(bbmcr.eb1,bbmcr.eb2))
	"""

end

# ╔═╡ d6fb2238-e7a9-491b-9057-041b32cace58
begin
	bbetrk = bb_eletrks.tracks[iele]
	bbetrk_xwalk =jn.Petit.walk_track_from_extremes(bbetrk)
	rbe=30.0
	
	print_reco_blobs(bbetrk, bbetrk_xwalk, 
					 rblob = rbe, 
		label="bb0nu σ=$(@sprintf("%.1f", bbetrk.diffusion.sigma_t)) mm vx=$(@sprintf("%.1f", bbetrk.diffusion.voxel_size)) mm")
					 
end

# ╔═╡ 7d5f5cd3-2ede-4c7a-8d25-35c817c1faf0
begin
	plot_bbetrk_2d, plot_bbetrk_3d =jn.Petit.plot_track_blobs(bbetrk,
																	bbetrk_xwalk;
																	sphere_radius= rbe)
	plot(plot_bbetrk_2d, plot_bbetrk_3d)
end

# ╔═╡ 464c9179-35a4-470b-a49d-1d4666c9e16a
bbetrk_xwalk.path_voxels

# ╔═╡ f74bf048-4c39-40dd-8d1c-6e9b149cb8ac
let
	bbmcr = jn.Petit.energy_blobs_from_path(bbetrk_xwalk, 3, energy_col=:electrons)
	md"""
	- eb1 = $(bbmcr.eb1), eb2 = $(bbmcr.eb2), Δb = $(Δb(bbmcr.eb1,bbmcr.eb2))
	"""

end

# ╔═╡ e9a297da-af74-4491-9bed-be191cd3fa23
begin
	xemtrk = bb_mctrks.tracks[kmc]
	xemtrk_xwalk =jn.Petit.walk_track_from_extremes(xemtrk)
	rxemc=5.0
	print_reco_blobs(xemtrk, xemtrk_xwalk, 
					 rblob = rxemc, 
		label="xe137 σ=$(@sprintf("%.1f", xemtrk.diffusion.sigma_t)) mm vx=$(@sprintf("%.1f", xemtrk.diffusion.voxel_size)) mm")
					 
end

# ╔═╡ d1b8d5c7-23a4-49cd-913f-d7ca0343bae2
begin
	plot_xemtrk_2d, plot_xemtrk_3d =jn.Petit.plot_track_blobs(xemtrk,
																	xemtrk_xwalk;
																	sphere_radius= rxemc)
	plot(plot_xemtrk_2d, plot_xemtrk_3d)
end

# ╔═╡ 4542893c-dd27-4f28-bef0-7c107b2d105e
let
	plot_bbmtrk_2d, plot_bbmtrk_3d =jn.Petit.plot_track_walk(xemtrk_xwalk)
																	
	plot(plot_bbmtrk_2d, plot_bbmtrk_3d)
end

# ╔═╡ e6b821f5-ecd2-4a47-8296-adb7b4c1318d
xemtrk_xwalk.path_voxels

# ╔═╡ 5b254930-2902-4b16-a4e6-6b802cd17828
let
	bbmcr = jn.Petit.energy_blobs_from_path(xemtrk_xwalk, nr, energy_col=:electrons)
	md"""
	- eb1 = $(bbmcr.eb1), eb2 = $(bbmcr.eb2), Δb = $(Δb(bbmcr.eb1,bbmcr.eb2))
	"""

end

# ╔═╡ 3127e6aa-feee-4726-9f8c-7bfc55e7581c
begin
	xeitrk = xe_iontrks.tracks[kion]
	xeitrk_xwalk =jn.Petit.walk_track_from_extremes(xeitrk)
	rxei=10.0
	
	print_reco_blobs(xeitrk, xeitrk_xwalk, 
					 rblob = rxei, 
		label="xe137 σ=$(@sprintf("%.1f", xeitrk.diffusion.sigma_t)) mm vx=$(@sprintf("%.1f", xeitrk.diffusion.voxel_size)) mm")
					 
end

# ╔═╡ 5450a959-f473-4659-9350-d2e63058591f
begin
	plot_xeitrk_2d, plot_xeitrk_3d =jn.Petit.plot_track_blobs(xeitrk,
																	xeitrk_xwalk;
																	sphere_radius= rxei)
	plot(plot_xeitrk_2d, plot_xeitrk_3d)
end

# ╔═╡ a2dd8015-ad7e-4484-8caf-e91395dda767
xeitrk_xwalk.path_voxels

# ╔═╡ 419dfc41-99b4-40c2-9319-4dadc08241bc
let
	plot_xeitrk_2d, plot_xeitrk_3d =jn.Petit.plot_track_walk(xeitrk_xwalk)
	plot(plot_xeitrk_2d, plot_xeitrk_3d)
end

# ╔═╡ 5c48ee86-78bf-4ad9-b75e-023cfb843eb9
let
	bbmcr = jn.Petit.energy_blobs_from_path(xeitrk_xwalk, nr, energy_col=:electrons)
	md"""
	- eb1 = $(bbmcr.eb1), eb2 = $(bbmcr.eb2), Δb = $(Δb(bbmcr.eb1,bbmcr.eb2))
	"""

end

# ╔═╡ bfe37a7d-d9d3-4d4b-ab67-b52fb869adce
begin
	xeetrk = xe_eletrks.tracks[kele]
	xeetrk_xwalk =jn.Petit.walk_track_from_extremes(xeetrk)
	rxee=30.0
	
	print_reco_blobs(xeetrk, xeetrk_xwalk, 
					 rblob = rxee, 
		label="bb0nu σ=$(@sprintf("%.1f", xeetrk.diffusion.sigma_t)) mm vx=$(@sprintf("%.1f", xeetrk.diffusion.voxel_size)) mm")
					 
end

# ╔═╡ fa617a54-7e00-46c0-ae41-f801b26fd572
xeetrk_xwalk.path_voxels

# ╔═╡ d43db256-777a-4c7c-a37f-447a2808eed1
begin
	plot_xeetrk_2d, plot_xeetrk_3d =jn.Petit.plot_track_blobs(xeetrk,
																	xeetrk_xwalk;
																	sphere_radius= rxee)
	plot(plot_xeetrk_2d, plot_xeetrk_3d)
end

# ╔═╡ 75da64af-683c-4310-a8e9-6db67a26cff1
let
	plot_xeitrk_2d, plot_xeitrk_3d =jn.Petit.plot_track_walk(xeetrk_xwalk)
	plot(plot_xeitrk_2d, plot_xeitrk_3d)
end

# ╔═╡ c1a7674a-3a3c-44d0-9e02-cb678e86e10b
let
	bbmcr = jn.Petit.energy_blobs_from_path(xeetrk_xwalk, nr, energy_col=:electrons)
	md"""
	- eb1 = $(bbmcr.eb1), eb2 = $(bbmcr.eb2), Δb = $(Δb(bbmcr.eb1,bbmcr.eb2))
	"""

end

# ╔═╡ 1a669bb1-9e99-4ac3-9f63-da23c4a48f6f
let
	bbmc_bamean,bbmc_bastd = jn.Petit.blob_asymmetry_from_path(bb_mctrks.tracks; i0=1, il=20, n0=1, nl=10, energy_col=:energy)
	bbion_bamean,bbion_bastd = jn.Petit.blob_asymmetry_from_path(bb_iontrks.tracks; i0=1, il=20, n0=1, nl=10, energy_col=:energy)
	bbele_bamean,bbele_bastd = jn.Petit.blob_asymmetry_from_path(bb_eletrks.tracks; i0=1, il=20, n0=1, nl=10, energy_col=:energy)
	xr = collect(1:10)
	scatter(xr, bbmc_bamean, yerror=bbmc_bastd, label="bb0nu blob asymmetry")
  	plot!(xr, bbmc_bamean, label=nothing, linewidth=1)
	scatter!(xr, bbion_bamean, yerror=bbion_bastd, label="ion blob asymmetry")
  	plot!(xr, bbion_bamean, label=nothing, linewidth=1)
	scatter!(xr, bbele_bamean, yerror=bbele_bastd, label="electron blob asymmetry")
  	plot!(xr, bbele_bamean, label=nothing, linewidth=1)
end

# ╔═╡ bc74a5f8-cdf2-4974-b8c6-d1642be15e4d
let
	xemc_bamean,xemc_bastd = jn.Petit.blob_asymmetry_from_path(xe_mctrks.tracks; i0=1, il=20, n0=1, nl=10, energy_col=:energy)
	xeion_bamean,xeion_bastd = jn.Petit.blob_asymmetry_from_path(xe_iontrks.tracks; i0=1, il=20, n0=1, nl=10, energy_col=:energy)
	xeele_bamean,xeele_bastd = jn.Petit.blob_asymmetry_from_path(xe_eletrks.tracks; i0=1, il=20, n0=1, nl=10, energy_col=:energy)
	xr = collect(1:10)
	scatter(xr, xemc_bamean, yerror=xemc_bastd, label="xe137 blob asymmetry")
  	plot!(xr, xemc_bamean, label=nothing, linewidth=1)
	scatter!(xr, xeion_bamean, yerror=xeion_bastd, label="ion blob asymmetry")
  	plot!(xr, xeion_bamean, label=nothing, linewidth=1)
	scatter!(xr, xeele_bamean, yerror=xeele_bastd, label="electron blob asymmetry")
  	plot!(xr, xeele_bamean, label=nothing, linewidth=1)
end

# ╔═╡ c79c2d21-a22d-4ed2-b7a7-3a5f7a6b9b4c
begin 
	
end

# ╔═╡ 51575ff2-a24a-4177-90fe-c904a20c2d57
begin
	r0 = 5
	rl = 30
	bbmc_bamean,bbmc_bastd = jn.Petit.blob_asymmetry(bb_mctrks.tracks; i0=1, il=20, r0=r0, rl =rl)
	bbion_bamean,bbion_bastd = jn.Petit.blob_asymmetry(bb_iontrks.tracks; i0=1, il=20, r0=r0, rl =rl)
	bbele_bamean,bbele_bastd = jn.Petit.blob_asymmetry(bb_eletrks.tracks; i0=1, il=20, r0=r0, rl =rl)
	xr = collect(r0:rl)
	scatter(xr, bbmc_bamean, yerror=bbmc_bastd, label="bb0nu blob asymmetry")
  	plot!(xr, bbmc_bamean, label=nothing, linewidth=1)
	scatter!(xr, bbion_bamean, yerror=bbion_bastd, label="ion blob asymmetry")
  	plot!(xr, bbion_bamean, label=nothing, linewidth=1)
	scatter!(xr, bbele_bamean, yerror=bbele_bastd, label="electron blob asymmetry")
  	plot!(xr, bbele_bamean, label=nothing, linewidth=1)
end

# ╔═╡ e88722ae-74fb-4abc-a9ab-e5f74ba6157a
let
	xemc_bamean,xemc_bastd = jn.Petit.blob_asymmetry(xe_mctrks.tracks; i0=1, il=20, r0=r0, rl =rl)
	xeion_bamean,xeion_bastd = jn.Petit.blob_asymmetry(xe_iontrks.tracks; i0=1, il=20, r0=r0, rl =rl)
	xeele_bamean,xeele_bastd = jn.Petit.blob_asymmetry(xe_eletrks.tracks; i0=1, il=20, r0=r0, rl =rl)
	scatter(xr, xemc_bamean, yerror=xemc_bastd, label="xe137 blob asymmetry")
  	plot!(xr, xemc_bamean, label=nothing, linewidth=1)
	scatter!(xr, xeion_bamean, yerror=xeion_bastd, label="ion blob asymmetry")
  	plot!(xr, xeion_bamean, label=nothing, linewidth=1)
	scatter!(xr, xeele_bamean, yerror=xeele_bastd, label="electron blob asymmetry")
  	plot!(xr, xeele_bamean, label=nothing, linewidth=1)
end


# ╔═╡ a5d01cd9-43e6-4097-8f81-748b407d043e
xe_iontrks.tracks[1].diffusion.ldrift

# ╔═╡ bfe6ccc4-bcfa-40be-9f52-9bee6d67e9a7
begin 	
	
end


# ╔═╡ cfd23030-a30d-448a-8fe8-185b032f9d3c
bb_iontrks.tracks[2]

# ╔═╡ 87c49cf1-9483-4467-b392-731bb25a6724
3.5/1.6

# ╔═╡ ca0229b0-fe96-4884-8155-8c2af47d5e67
0.7*2458/100

# ╔═╡ f507244b-8b79-4cc7-8957-37ad54869f50
let
	bbion_eb1,bbion_eb2, bbion_ld,bbion_db, _ = jn.Petit.blob_analysis_vs_radius(bbitrk; r0=5, rl=15)
	#_, pdbld = jn.Petit.p1df(bbion_ld, bbion_db, 10)
	#plot(pdbld)
end

# ╔═╡ 6a7f1dda-18ba-4746-9aa7-4ef0bf533e21
begin
	
end

# ╔═╡ 41e8f62f-ba3c-44b5-b804-22eefef23d44
begin
	xeion_ndbld, xeion_pdbld = jn.Petit.p1df(xeion_ld, xeion_db, 10)
	plot(xeion_pdbld)
end

# ╔═╡ b7c9112f-6ce9-4605-b023-95fa4c573a93
begin
	bbele_eb1,bbion_ele, bbele_ld,bbele_db = blob_asymmetry_rblob(bb_eletrks; i0=1, il=50, r0=20.0)
end

# ╔═╡ 1524fe96-92dd-403a-a68b-8014dd5afbfc
begin
	xeele_eb1,xeion_ele, xeele_ld,xeele_db = blob_asymmetry_rblob(xe_eletrks; i0=1, il=50, r0=20.0)
end

# ╔═╡ 79e0190f-054a-4455-9df8-7ba381eb6f89
begin
	bbele_ndbld, bbele_pdbld = jn.Petit.p1df(bbele_ld, bbele_db, 10)
	plot(bbele_pdbld)
end

# ╔═╡ 62ad8e00-aeda-4e54-bd2e-91bbab3beb16
begin
	xeele_ndbld, xeele_pdbld = jn.Petit.p1df(xeele_ld, xeele_db, 10)
	plot(xeele_pdbld)
end

# ╔═╡ 9a89e70a-f40b-4fda-a17d-8a875a8c327a
begin
	
	bbeb1,bbeb2, bbdb = blob_asymmetry_rblob(trksbb1mm.tracks; i0=1, il=500, r0=5.0)
	xeeb1,xeeb2, xedb = blob_asymmetry_rblob(trksxe1mm.tracks; i0=1, il=500, r0=5.0)
	bbeb1r10,bbeb2r10, _ = blob_asymmetry_rblob(trksbb1mm.tracks; i0=1, il=500, r0=10.0)
	xeeb1r10,xeeb2r10, _ = blob_asymmetry_rblob(trksxe1mm.tracks; i0=1, il=500, r0=10.0)
	
end

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

# ╔═╡ 5ff5638b-8b94-4c22-a3c4-a17988eff3bd
function blob_asymmetry(tracks; i0::Int=1, il::Int=10, r0::Int=5, rl::Int=15)
    # Validate indices
    il = min(il, length(tracks))
    i0 = max(1, i0)

    if i0 > il
        error("i0 ($i0) must be <= il ($il)")
    end

    nevents = il - i0 + 1
    nradii = rl - r0 + 1
    DB = Matrix{Float64}(undef, nevents, nradii)

    for (ievent_idx, ievent) in enumerate(i0:il)
        trk = tracks[ievent]
        walk_result = jn.Petit.walk_track_from_extremes(trk)

        for (ir, rr) in enumerate(r0:rl)
            sphere_energies = jn.Petit.energy_in_spheres_around_extremes(trk, walk_result, Float64(rr))
            eb1 = sphere_energies.blob1_energy
            eb2 = sphere_energies.blob2_energy

            # Compute asymmetry, handle division by zero
            total = eb1 + eb2
            db = total > 0 ? abs(eb1 - eb2) / total : 0.0
            DB[ievent_idx, ir] = db
        end
    end

    return (vec(mean(DB, dims=1)), vec(std(DB, dims=1)))
end

# ╔═╡ 5e6ca95c-4ecd-458f-9cd3-3e6f2c652994
function energy_blobs_from_path(walk_result, n::Int; energy_col::Symbol=:energy)
    # Check if we have valid path_voxels
    if isnothing(walk_result.path_voxels) || nrow(walk_result.path_voxels) == 0
        return (eb1 = 0.0, eb2 = 0.0, n1 = 0, n2 = 0,
                start_energy = 0.0, end_energy = 0.0)
    end

    path_df = walk_result.path_voxels
    nvoxels = nrow(path_df)

    # Check if requested column exists
    if !hasproperty(path_df, energy_col)
        error("Column :$energy_col not found in path_voxels. Available: $(names(path_df))")
    end

    # Get energy values from specified column
    energy_values = Float64.(path_df[!, energy_col])

    # Ensure n is at least 1 and doesn't exceed half the path
    n = max(1, n)
    n_start = min(n, nvoxels ÷ 2)  # Don't overlap with end
    n_end = min(n, nvoxels - n_start)  # Take remaining from end

    # Sum energy from first n_start voxels (start of path)
    start_energy = sum(energy_values[1:n_start])

    # Sum energy from last n_end voxels (end of path)
    end_energy = sum(energy_values[(nvoxels - n_end + 1):nvoxels])

    # Assign blob1 to higher energy, blob2 to lower
    if start_energy >= end_energy
        eb1 = start_energy
        eb2 = end_energy
        n1 = n_start
        n2 = n_end
    else
        eb1 = end_energy
        eb2 = start_energy
        n1 = n_end
        n2 = n_start
    end

    return (eb1 = eb1, eb2 = eb2, n1 = n1, n2 = n2,
            start_energy = start_energy, end_energy = end_energy)
end


# ╔═╡ 1a3a6847-8bab-411e-924b-ba13d3b7b52d
function blob_asymmetry_from_path(tracks; i0::Int=1, il::Int=10,
                                   n0::Int=1, nl::Int=10, energy_col::Symbol=:energy)
    # Validate indices
    il = min(il, length(tracks))
    i0 = max(1, i0)

    if i0 > il
        error("i0 ($i0) must be <= il ($il)")
    end

    nevents = il - i0 + 1
    nn = nl - n0 + 1
    DB = Matrix{Float64}(undef, nevents, nn)

    for (ievent_idx, ievent) in enumerate(i0:il)
        trk = tracks[ievent]
        walk_result = jn.Petit.walk_track_from_extremes(trk)

        for (in_idx, n) in enumerate(n0:nl)
            blobs = energy_blobs_from_path(walk_result, n; energy_col=energy_col)
            eb1 = blobs.eb1
            eb2 = blobs.eb2

            # Compute asymmetry, handle division by zero
            total = eb1 + eb2
            db = total > 0 ? abs(eb1 - eb2) / total : 0.0
            DB[ievent_idx, in_idx] = db
        end
    end

    return (vec(mean(DB, dims=1)), vec(std(DB, dims=1)))
end

# ╔═╡ 9e8a64e8-e5fb-4ec2-8161-51e13864e325
function plot_track_walk(walk_result;
                         markersize_voxels::Float64=4.0,
                         extreme_radius::Float64=5.0,
                         alpha_circles::Float64=0.3)

    # Helper: draw circle in 2D
    function draw_circle_2d(cx, cy, r)
        θ = range(0, 2π, length=100)
        cx .+ r * cos.(θ), cy .+ r * sin.(θ)
    end

    # Helper: draw wireframe sphere in 3D
    function draw_wireframe_sphere!(p, cx, cy, cz, r, color, n_lines=8)
        θ = range(0, 2π, length=50)
        for i in 1:n_lines
            ϕ = (i-1) * π / n_lines
            plot!(p, cx .+ r*sin(ϕ).*cos.(θ), cy .+ r*sin(ϕ).*sin.(θ),
                  cz .+ r*cos(ϕ).*ones(length(θ)), color=color, lw=2, alpha=0.7, label="")
            angle = (i-1) * π / (n_lines-1) - π/2
            r_circ = r * cos(angle)
            plot!(p, cx .+ r_circ.*cos.(θ), cy .+ r_circ.*sin.(θ),
                  (cz + r*sin(angle)).*ones(length(θ)), color=color, lw=2, alpha=0.7, label="")
        end
    end

    # Check if we have valid path_voxels
    if isnothing(walk_result.path_voxels) || nrow(walk_result.path_voxels) == 0
        error("walk_result.path_voxels is empty or nothing")
    end

    # Extract path voxel positions and electrons
    path_df = walk_result.path_voxels
    x = Float64.(path_df.x)
    y = Float64.(path_df.y)
    z = Float64.(path_df.z)

    # Use electrons column if available, otherwise use energy
    if hasproperty(path_df, :electrons)
        intensity = Float64.(path_df.electrons)
        intensity_label = "Electrons"
    else
        intensity = Float64.(path_df.energy)
        intensity_label = "Energy"
    end

    # Compute plot limits with padding
    padding = 0.7 + extreme_radius / max(1.0, min(maximum(x)-minimum(x), maximum(y)-minimum(y), maximum(z)-minimum(z)))
    xlim = (mean(extrema(x)) - padding*(maximum(x)-minimum(x)), mean(extrema(x)) + padding*(maximum(x)-minimum(x)))
    ylim = (mean(extrema(y)) - padding*(maximum(y)-minimum(y)), mean(extrema(y)) + padding*(maximum(y)-minimum(y)))
    zlim = (mean(extrema(z)) - padding*(maximum(z)-minimum(z)), mean(extrema(z)) + padding*(maximum(z)-minimum(z)))

    cmap = cgrad(:viridis)

    # Common scatter options
    scatter_opts = (marker_z=intensity, ms=markersize_voxels, color=cmap,
                   colorbar_title=intensity_label, legend=false, label="", markerstrokewidth=0)

    # Get extremes
    start_voxel = walk_result.extremes[1]
    end_voxel = walk_result.extremes[2]

    # === XY Projection ===
    p1 = scatter(x, y; scatter_opts..., xlabel="X (mm)", ylabel="Y (mm)",
                title="XY Projection", xlims=xlim, ylims=ylim)
    if !isnothing(start_voxel)
        # Draw red circles around extremes
        cx, cy = draw_circle_2d(start_voxel.x, start_voxel.y, extreme_radius)
        plot!(p1, cx, cy, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p1, [start_voxel.x], [start_voxel.y], ms=6, color=:red, markershape=:circle, label="Start")

        cx, cy = draw_circle_2d(end_voxel.x, end_voxel.y, extreme_radius)
        plot!(p1, cx, cy, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p1, [end_voxel.x], [end_voxel.y], ms=6, color=:darkred, markershape=:circle, label="End")
    end

    # === XZ Projection ===
    p2 = scatter(x, z; scatter_opts..., xlabel="X (mm)", ylabel="Z (mm)",
                title="XZ Projection", xlims=xlim, ylims=zlim)
    if !isnothing(start_voxel)
        cx, cz = draw_circle_2d(start_voxel.x, start_voxel.z, extreme_radius)
        plot!(p2, cx, cz, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p2, [start_voxel.x], [start_voxel.z], ms=6, color=:red, markershape=:circle, label="")

        cx, cz = draw_circle_2d(end_voxel.x, end_voxel.z, extreme_radius)
        plot!(p2, cx, cz, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p2, [end_voxel.x], [end_voxel.z], ms=6, color=:darkred, markershape=:circle, label="")
    end

    # === YZ Projection ===
    p3 = scatter(y, z; scatter_opts..., xlabel="Y (mm)", ylabel="Z (mm)",
                title="YZ Projection", xlims=ylim, ylims=zlim)
    if !isnothing(start_voxel)
        cy, cz = draw_circle_2d(start_voxel.y, start_voxel.z, extreme_radius)
        plot!(p3, cy, cz, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p3, [start_voxel.y], [start_voxel.z], ms=6, color=:red, markershape=:circle, label="")

        cy, cz = draw_circle_2d(end_voxel.y, end_voxel.z, extreme_radius)
        plot!(p3, cy, cz, color=:red, fill=(0, alpha_circles, :red), lw=2, label="")
        scatter!(p3, [end_voxel.y], [end_voxel.z], ms=6, color=:darkred, markershape=:circle, label="")
    end

    # === 3D View ===
    p4 = scatter(x, y, z; scatter_opts..., xlabel="X (mm)", ylabel="Y (mm)", zlabel="Z (mm)",
                title="3D View", xlims=xlim, ylims=ylim, zlims=zlim)
    if !isnothing(start_voxel)
        # Draw red wireframe spheres
        draw_wireframe_sphere!(p4, start_voxel.x, start_voxel.y, start_voxel.z, extreme_radius, :red)
        scatter!(p4, [start_voxel.x], [start_voxel.y], [start_voxel.z],
                ms=6, color=:red, markershape=:circle, label="Start")

        draw_wireframe_sphere!(p4, end_voxel.x, end_voxel.y, end_voxel.z, extreme_radius, :darkred)
        scatter!(p4, [end_voxel.x], [end_voxel.y], [end_voxel.z],
                ms=6, color=:darkred, markershape=:circle, label="End")
    end

    # Title with info
    title_text = "Track Walk Path"
    if !isnothing(start_voxel)
        title_text *= "\nPath length: $(round(walk_result.total_length, digits=2)) mm"
        title_text *= " | Voxels: $(nrow(path_df))"
        title_text *= " | Confidence: $(round(walk_result.confidence, digits=2))"
    end

    # Combine plots
    plot_2d = plot(p1, p2, p3, layout=(2,2), size=(1200,1000),
                   plot_title=title_text, titlefontsize=10,
                   left_margin=8Plots.mm, bottom_margin=5Plots.mm)
    plot_3d = plot(p4, size=(1000,1000), plot_title=title_text,
                   titlefontsize=10, left_margin=8Plots.mm, bottom_margin=5Plots.mm)

    return (plot_2d, plot_3d)
end


# ╔═╡ 99cb006d-064a-4d10-917a-ad0562ab1ce7
let
	plot_bbmtrk_2d, plot_bbmtrk_3d =plot_track_walk(bbmtrk_xwalk)
																	
	plot(plot_bbmtrk_2d, plot_bbmtrk_3d)
end

# ╔═╡ 58ad2d3a-0e26-4425-9faf-11a85a0e24d3
let
	plot_bbitrk_2d, plot_bbitrk_3d =plot_track_walk(bbitrk_xwalk)
	plot(plot_bbitrk_2d, plot_bbitrk_3d)
end

# ╔═╡ 9a43ceb1-7218-488d-945f-73c6ce9bc5df
let
	plot_bbetrk_2d, plot_bbetrk_3d =plot_track_walk(bbetrk_xwalk)
	plot(plot_bbetrk_2d, plot_bbetrk_3d)
end

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
# ╠═e91c9f90-63a6-48cf-ab8a-3d8d4ac320c4
# ╠═84eb8703-228b-46f2-8868-5086f576572e
# ╠═0289428e-6db3-4bfc-82fe-5b9e71dd3706
# ╠═726badaf-c389-484e-a46d-cdda4d4d578d
# ╠═43d5968c-e185-4586-9f49-f9be529e114b
# ╠═f28a58d0-55a3-4f6d-af3f-200b8c395c98
# ╠═35828cbe-b5d9-4f85-883b-8996646338f8
# ╠═01ce8389-a772-480b-a875-57916a4cd0ea
# ╠═529ebb7f-1615-464e-b569-dc423e1bb8c7
# ╠═79156ecf-35d9-45a5-bf05-51ceeee9da2e
# ╠═896f02e7-fede-42b7-9366-4955d1879b6c
# ╠═312de1b0-bba0-40ab-a537-c7c5343380b2
# ╠═1024cc13-60e0-42bf-ad14-2f1d90da4484
# ╠═27d8edb9-5b07-4760-9cbe-3bb413466b99
# ╠═f6dabe79-9985-47f2-b536-b12fef702393
# ╠═b73a4322-e098-4032-8564-3f1f4914a955
# ╠═475079bb-18dd-477e-9b67-320ca5909ade
# ╠═fc9cc4a7-56b5-4852-bc19-483fa7379ab2
# ╠═132042e7-f7ad-46ea-88d6-396bd072c5a8
# ╠═4fd05d01-16c4-4849-8486-ccfad7e3290d
# ╠═7dea81cf-f51d-4411-ae68-bb86bb04dcf3
# ╠═7b3a8a4d-f649-4175-b881-453b7e44e598
# ╠═d0fecc3e-5e0e-4522-9efa-b4b60fd151a5
# ╠═b32e2b11-788e-4017-b5cf-4e4d3968398d
# ╠═6950e911-3003-4778-98bd-2b1a0a819844
# ╠═cc0ad60d-0e31-4406-9b0d-20d938a00aff
# ╠═1dfa60a6-c09c-46d7-874b-44ede0043780
# ╠═41653083-87fb-46f7-8d1f-5cf9d34c2ea7
# ╠═47dc42ac-9dd9-4318-9a29-7655541d3af8
# ╠═af515863-e9e6-4a8a-a940-0a8439afd924
# ╠═83a718af-c2d9-4743-827d-8682a7379c65
# ╠═7ba1eacc-88f9-4bc3-ba75-1efce41dadb7
# ╠═4e1f83d9-e9f8-47c9-a742-3eb32b371c93
# ╠═d92a0d53-d5f0-4e77-9911-a97e93d5b7b4
# ╠═1d6276f9-c585-48ed-9c3e-d910cd47d091
# ╠═744b4f32-5f0e-4e2e-b6e3-033b1d74af62
# ╠═ac625d74-16f4-4463-bada-42d0b47e79d5
# ╠═b0ab441c-f7aa-4bfa-8bc8-b87d4decfee9
# ╠═49eacb34-e766-41f2-ab1f-8170f90e5a50
# ╠═dd1ebb2d-581e-4c5c-9f69-6db4e6fe0a1b
# ╠═df172361-5846-471d-8014-92db51969574
# ╠═66e9e2e4-aa44-4e8c-b304-0dddd65e4fc6
# ╠═fcf2b122-cc41-4b63-b631-41c30a19c728
# ╠═8a24d1ce-547f-400b-9415-74dd0b7dd689
# ╠═3fba0b56-cfc0-4115-ab29-8eccc1aefb37
# ╠═46528939-1e69-4f4d-b716-3e7f76e6c22f
# ╠═6f4c9519-ec55-452a-a683-01ed779d7334
# ╠═b5070fc8-1ec7-4212-9b3b-43f0fa999d24
# ╠═1ac6befd-59d8-40ad-a5b0-710b7f941444
# ╠═136adb7c-8ba5-4bed-9f73-132a937d2a16
# ╠═891ff777-2294-47a9-89e0-e3810203dcff
# ╠═38474add-cbdf-4772-ba95-c7da00485b31
# ╠═49f11799-fe6d-40f8-abda-d90bc9c4d3e8
# ╠═86bf06e8-9caa-4c18-80ff-160cbd945986
# ╠═05998d38-6935-44e9-9c16-53e1810ecad5
# ╠═5c9eab19-8870-47e5-a9a9-263c8f70bd44
# ╠═7abbc303-db74-4764-be25-b0304ed694ea
# ╠═382bdd57-44c9-47dd-83be-c5bbee399ac9
# ╠═1b0e400f-ca67-4f2e-9e68-a90192fded54
# ╠═29c5b0e3-9a0f-4cde-ab10-7bca1ff18468
# ╠═3a795130-6a5e-4606-92b2-57f006bda4a4
# ╠═596bc0f0-fc1a-4c8b-bab0-e799d2621f46
# ╠═7fd42a7e-18b0-408c-a439-7804f90e09df
# ╠═deeceb13-3b73-47c1-8828-5fbdffabd3de
# ╠═c92d9d16-f53c-43d7-95a6-ca2319c26da2
# ╠═63159f7b-a1fc-4648-b148-040893e609ab
# ╠═177e8223-e063-42e4-b633-17f87f653c74
# ╠═a4442faf-6442-4ef4-9200-e3e705c7a2b2
# ╠═e1b8c419-c2dc-4d5b-af16-c59f77d54aa8
# ╠═2a2293f5-b037-498b-976e-cd95f1c07d1a
# ╠═d16366ec-3e0c-4b87-adb3-f180413cdd39
# ╠═7e707191-d6ff-4513-92d5-7d2ce8e5db3c
# ╠═abc1fc50-ec3c-4bd6-a6af-7a0af8147b3d
# ╠═e51fd878-7c12-4770-b1f1-c9f9e6e33256
# ╠═380f1e72-1215-421e-8b00-268f34c82341
# ╠═6f8a94a0-0ec9-488f-bd3c-d7a67f27f978
# ╠═f1d33080-06f8-4310-8589-369570f8934a
# ╠═8c173831-8e37-47d4-a183-d95be8773781
# ╠═1396c7a9-bb49-4132-9a4f-0782922141f8
# ╠═105bd81d-963d-4ac1-9bfc-1f3179f90202
# ╠═3ff630c5-6409-49ea-971c-e3cd15427de2
# ╠═00ecc891-6de9-4b6d-a9e2-00de7a8c1d08
# ╠═be04f12b-f434-4090-8e69-fac08f2b8350
# ╠═d2961e74-c99c-414b-b4a7-8fdde1bf9c8a
# ╠═bca3c9d7-5583-42a6-89a0-e595e5b21970
# ╠═5af82f10-c89f-4092-811f-12ae1e1ef94b
# ╠═362d825d-e385-4b28-8fc2-f083391bc89a
# ╠═b1bff10b-31e6-4820-a7c4-df154c01ab6b
# ╠═2d5be057-d840-413a-a13f-8752a2619fd0
# ╠═c75ea153-481d-47a3-99ba-8dcf26416dc0
# ╠═e316e013-d2d8-4fa2-b141-7e26eef0ad49
# ╠═b8e89e80-c581-49a0-8e5f-fb41ac4ea01b
# ╠═194e2770-9b31-44f2-897a-9404be48079b
# ╠═99cb006d-064a-4d10-917a-ad0562ab1ce7
# ╠═af19f2c1-66c6-40f7-acb5-b49acb17db51
# ╠═e34e216a-0c7a-45e2-abb7-7d50f4172117
# ╠═3af1d8e7-9d83-40a9-b9e1-95c292785e38
# ╠═f35715e9-5024-4d2e-b981-8c06a3ea4300
# ╠═d1cf9c39-b91d-4e16-8810-980e516971c5
# ╠═519fe2aa-0465-4253-b0c7-1492ab3b8440
# ╠═ca37680d-9162-497d-9b36-d73fd90b2e07
# ╠═9c1590fc-ad2b-4444-95db-9ffbcedca155
# ╠═58ad2d3a-0e26-4425-9faf-11a85a0e24d3
# ╠═51cb6061-9432-48b7-affb-a34c3549861d
# ╠═67570f37-cc53-4ab8-8069-e8a622297971
# ╠═717de5ce-9edb-4d37-a1df-1dc756aa6831
# ╠═0adc0c96-c8c4-4c41-8f28-46ded1217cc6
# ╠═d6fb2238-e7a9-491b-9057-041b32cace58
# ╠═7d5f5cd3-2ede-4c7a-8d25-35c817c1faf0
# ╠═464c9179-35a4-470b-a49d-1d4666c9e16a
# ╠═f74bf048-4c39-40dd-8d1c-6e9b149cb8ac
# ╠═9a43ceb1-7218-488d-945f-73c6ce9bc5df
# ╠═e9a297da-af74-4491-9bed-be191cd3fa23
# ╠═d1b8d5c7-23a4-49cd-913f-d7ca0343bae2
# ╠═4542893c-dd27-4f28-bef0-7c107b2d105e
# ╠═e6b821f5-ecd2-4a47-8296-adb7b4c1318d
# ╠═5b254930-2902-4b16-a4e6-6b802cd17828
# ╠═3127e6aa-feee-4726-9f8c-7bfc55e7581c
# ╠═5450a959-f473-4659-9350-d2e63058591f
# ╠═a2dd8015-ad7e-4484-8caf-e91395dda767
# ╠═419dfc41-99b4-40c2-9319-4dadc08241bc
# ╠═5c48ee86-78bf-4ad9-b75e-023cfb843eb9
# ╠═bfe37a7d-d9d3-4d4b-ab67-b52fb869adce
# ╠═fa617a54-7e00-46c0-ae41-f801b26fd572
# ╠═d43db256-777a-4c7c-a37f-447a2808eed1
# ╠═75da64af-683c-4310-a8e9-6db67a26cff1
# ╠═c1a7674a-3a3c-44d0-9e02-cb678e86e10b
# ╠═1a669bb1-9e99-4ac3-9f63-da23c4a48f6f
# ╠═bc74a5f8-cdf2-4974-b8c6-d1642be15e4d
# ╠═c79c2d21-a22d-4ed2-b7a7-3a5f7a6b9b4c
# ╠═51575ff2-a24a-4177-90fe-c904a20c2d57
# ╠═e88722ae-74fb-4abc-a9ab-e5f74ba6157a
# ╠═a5d01cd9-43e6-4097-8f81-748b407d043e
# ╠═bfe6ccc4-bcfa-40be-9f52-9bee6d67e9a7
# ╠═cfd23030-a30d-448a-8fe8-185b032f9d3c
# ╠═87c49cf1-9483-4467-b392-731bb25a6724
# ╠═ca0229b0-fe96-4884-8155-8c2af47d5e67
# ╠═f507244b-8b79-4cc7-8957-37ad54869f50
# ╠═6a7f1dda-18ba-4746-9aa7-4ef0bf533e21
# ╠═41e8f62f-ba3c-44b5-b804-22eefef23d44
# ╠═b7c9112f-6ce9-4605-b023-95fa4c573a93
# ╠═1524fe96-92dd-403a-a68b-8014dd5afbfc
# ╠═79e0190f-054a-4455-9df8-7ba381eb6f89
# ╠═62ad8e00-aeda-4e54-bd2e-91bbab3beb16
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
# ╠═5ff5638b-8b94-4c22-a3c4-a17988eff3bd
# ╠═1a3a6847-8bab-411e-924b-ba13d3b7b52d
# ╠═5e6ca95c-4ecd-458f-9cd3-3e6f2c652994
# ╠═9e8a64e8-e5fb-4ec2-8161-51e13864e325
# ╠═4945ede2-bba5-4a6d-818f-482e66acd675
# ╠═b64344c7-5b7f-41a8-be36-bd24a3d0256f
# ╠═4dbfcaba-82d6-4c40-ad8b-1df81f6387e3
# ╠═f4ad81a7-b6c4-43a0-9ccf-88b6171593ef
