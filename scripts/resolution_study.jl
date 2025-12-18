#!/usr/bin/env julia


# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using Printf
using DataFrames
using CSV
using Glob
using ArgParse

# Load Petit module
include(joinpath(pdir, "src", "Petit.jl"))
import .Petit
# Helper functions (get_sigma, get_voxel_size_and_distance, etc.)
# are now part of Petit module via itaca_functions.jl

cmdir=joinpath(ENV["DATA"], "HD5t/itaca")



function continue_on_enter(plt)
        display(plt)
        println("Press Enter to continue...")
        readline()
    end
#=============================================================================
# Main Script
=============================================================================#



function print_blobs(blobs)
    println("Blob Analysis:")
    println("  Blob1 (high energy):")
    println("    position   = ($(round(blobs.blob1.x, digits=2)), $(round(blobs.blob1.y, digits=2)), $(round(blobs.blob1.z, digits=2))) mm")
    println("    energy     = $(round(blobs.Eb1*1e3, digits=1)) keV")
    println("  Blob2 (low energy):")
    println("    position   = ($(round(blobs.blob2.x, digits=2)), $(round(blobs.blob2.y, digits=2)), $(round(blobs.blob2.z, digits=2))) mm")
    println("    energy     = $(round(blobs.Eb2*1e3, digits=1)) keV")
    println("  Asymmetry    = $(round(blobs.asymmetry, digits=3))")
end

function process_event(ievent::Int, bbdf::DataFrame, ion_dfpars::Petit.DiffusionParams,
                       σt_ion_p15::Float64, ion_voxel::Float64, max_distance::Float64,
                       energy_threshold_keV::Float64, Rb::Float64)

    println("\n========== Processing event $ievent ==========")

    bbevt = Petit.get_event(bbdf, ievent)
    bbevtmc = transform_hits_df(bbevt)
    plt = Petit.plot_event(bbevtmc)
    continue_on_enter(plt)

    bbion_df = Petit.diffuse_xyz_image_mc(bbevtmc;
                                          sigma_t_mm=σt_ion_p15,
                                          sigma_l_mm=σl_ion,
                                          nbins=nbins,
                                          nsigma=nsigma)

    bbion_vx = Petit.voxelize_event(bbion_df, ion_voxel)

    bbdfion_tracks = Petit.make_tracks(bbion_vx;
                                       max_distance_mm=max_distance,
                                       energy_threshold_kev=energy_threshold_keV,
                                       diffusion=ion_dfpars)

    println("- number of tracks found = $(length(bbdfion_tracks))")

    if length(bbdfion_tracks) != 1
        println("Number of tracks found not one, skipping event $ievent")
        return
    end

    bbdfiont1 = bbdfion_tracks[1]
    bbdfiont1_walk = Petit.walk_track_from_extremes(bbdfiont1)

    plt = Petit.plot_track_with_extremes(bbdfiont1, bbdfiont1_walk)
    continue_on_enter(plt)

    bbdfionCP = Petit.reconstruct_central_path(bbdfiont1,
                                               bbdfiont1_walk.path_indices;
                                               filter_radius = nsigma*σt_ion_p15)

    plt = Petit.plot_reco_track_with_voxels(bbdfiont1, bbdfionCP)
    continue_on_enter(plt)

    bbdfion_blobs = Petit.find_blob_energies(bbdfiont1,
                                             bbdfionCP;
                                             radius=Rb)

    plt = Petit.plot_reco_track_with_voxels_and_spheres(bbdfiont1, bbdfionCP,
                                                        bbdfion_blobs, Rb)
    continue_on_enter(plt)
    print_blobs(bbdfion_blobs)
end


function main(; ievent::Int=1, levent::Int=1, ldrft::Float64=100.0,
               voxel_scale::Float64=1.0, voxel_distance_scale::Float64=1.5)
    
    # main parameters
    input_file="bb0nu/bb0nu_15bar_p1.h5" # get from command line
    particle_type = "ion"  # get from command line
    ldrft = 100.0 # get fro command line
    dt = 3.5 # get from command line.
    dl = 0.9 # get from command line.
    tK = 297.0 # get from command line. 
    edrift = 500.0 # get from command line. 
    Pbar=15.0 # get from command line.
    energy_threshold_ions =  10.0 # get from command line.
    energy_threshold_keV =  10.0 # get from command line.
    nbins = 100 # get from command line.
    nsigma = 3.0  # get from command line.
    n_kde_eval = 200  # get from command line.
    Rb = 10.0 # get from command line.
    nevent = 5 # get from command line.

    σt,  σl=  Petit.get_sigma(particle_type, ldrft; 
                        dt = dt, dl = dl, 
                        tK = tk, edrift = edrift, Pbar=Pbar)

    voxel_size, mcvox_size, max_distance = Petit.get_voxel_size_and_distance(ldrft, σt)

    energy_threshold_keV  = Petit.get_energy_threshold(particle_type; 
                                                    energy_threshold_ions =  energy_threshold_ions,
                                                    energy_threshold_keV =  energy_threshold_keV)

    kde_bandwidth =2 * voxel_size
    mc_kde_bandwidth=2 * voxel_size
    

    dfpars = Petit.DiffusionParams(ldrft, σt, σl, voxel_size, max_distance, 
                                        energy_threshold_keV, nbins, nsigma)


    diffusion_params_print(dfpars)

    # Load event
    md"- Examine event $(nevent)" # change by println 

  	hitsdf = load_data(input_file, cmdir)

    event_df = jn.Petit.get_event(hitsdf, nevent)
	plt = Petit.plot_event(event_df)
    continue_on_enter(plt) 

    # replace by println 
        md"""### Compute MC path 
    - MC path (voxelized primary particle trajectory with arc-length)
    - First/last rows give MC truth extremes
    """

    mc_path = Petit.compute_mc_path(event_df, mcvox_size)
	plt = Petit.plot_event(event_df; mc_path=mc_path)
    continue_on_enter(plt) 

    md"### Transform hits (removes time, particle_id, etc.)" #println

    event_mc = Petit.transform_hits_df(event_df)
    
    md"### Diffuse event" #println 

    diffused_df = Petit.diffuse_xyz_image_mc(event_mc;
                                             sigma_t_mm=σt,
                                             sigma_l_mm=σl,
                                             nbins=nbins,
                                             nsigma=nsigma)

    #
    plt = Petit.plot_hits(diffused_df, energy_column=:electrons)
    continue_on_enter(plt)

    md"### Voxelize event" # println

    voxels = Petit.voxelize_event(diffused_df, voxel_size)
	plt = Petit.plot_event(voxels)
    continue_on_enter(plt) 

    md"### Make tracks" #println 

    tracks = Petit.make_tracks(voxels;
                               max_distance_mm=max_distance,
                               energy_threshold_kev=energy_threshold_keV,
                               diffusion=dfpars)
    md"""
    - number of tracks found = $(length(tracks)) 
    """

    if length(tracks >1)
        println("Number of tracks = $(length(tracks)), stop.")
        # save plots and metadata including event number, run conditions (all parameters)
        # and diffusion parameter. Create a dir --outdir, by default in the directory where
        # we run (scripts), save plots as pngs with meaningful names in directory, as well as metadata in csv format.
        # then exit
    end

    track = tracks[1]
	walk_result = Petit.walk_track_from_extremes(track)
	walk_result_md(walk_result) # replace by println

    path = Petit.get_raw_path(track, walk_result.path_indices)
	path_md(path::DataFrame) # replace by println

    extreme_dists = Petit.compute_extreme_distances(path, mc_path)
	extreme_distances_md(extreme_dists) #println
    # As part of ours results, we want to save distances d1 and d2.

    plt = Petit.plot_track_with_paths(track, path, mc_path;
                              show_distances=false)
    continue_on_enter(plt) 

    md"""
    ### KDE
    """

    reco_kde = Petit.get_reco_kde(track, path; bandwidth=kde_bandwidth, 
									  n_eval=n_kde_eval)
	mc_kde = Petit.get_mc_kde(mc_path; bandwidth=mc_kde_bandwidth, 
								  n_eval=n_kde_eval)

    # replace below by function
    p1 = plot(reco_kde.kde_s, reco_kde.kde_f,
              xlabel="Arc length s (mm)",
              ylabel=" KDE Energy density f(s)",
              title="Longitudinal Energy Density (Event $(nevent))",
              label="RECO KDE (h=$(round(kde_bandwidth, digits=1))mm)",
              linewidth=2, 
              color=:blue)
	p2 = plot!(p1, mc_kde.kde_s, mc_kde.kde_f,
              xlabel="Arc length s (mm)",
              ylabel=" KDE Energy density f(s)",
              title="Longitudinal Energy Density (Event $(nevent))",
              label="MC KDE (h=$(round(kde_bandwidth, digits=1))mm)",
              linewidth=2, 
              color=:red)
	p3 = histogram(reco_kde.s, weights=reco_kde.E .* 1e3,
                   bins=range(0, reco_kde.track_length, length=n_kde_eval+1),
                   xlabel="Arc length s (mm)",
                   ylabel="Energy (keV)",
                   title="Energy vs S (Event $(nevent)",
                   label="RECO",
                   fillalpha=0.7,
                   color=:gray)
	p4 = histogram!(p3, mc_kde.s, weights=mc_kde.E .* 1e3,
                   bins=range(0, mc_kde.track_length, length=n_kde_eval+1),
                   xlabel="Arc length s (mm)",
                   ylabel="Energy (keV)",
                   title="Energy vs S (MC)",
                   label="MC",
                   fillalpha=0.7,
                   color=:green)
	plot(p2,p4, layout=(1, 2), size=(1200, 1000))

    # enter to continue

    md"### Find peaks " 

    peaks = Petit.find_peaks(reco_kde.kde_f, reco_kde.kde_s; prom_scale=0.2)
	plt = Petit.plot_kde_peaks(reco_kde.kde_s, reco_kde.kde_f, peaks; title="RECO KDE")
    continue_on_enter(plt) 

    # need to save the information on peaks. proms, positions, width, lefts and rights
    # for the peak whith the leftmost left and the peak with the rigtmost right. 

        md"""
    ### Blob analysis
    """

    blobs = Petit.find_blob_energies(track,
                                    path;
                                    radius=Rb) 

    plt = Petit.plot_reco_track_with_voxels_and_spheres(track, path, blobs, Rb)
    continue_on_enter(plt)

    # need to save blobs information. Eblo1, eblob2, asymmetry. 
    
end



# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--ievent", "-i"
            help = "Initial event number"
            arg_type = Int
            default = 1
        "--levent", "-e"
            help = "Last event number"
            arg_type = Int
            default = 1
        "--ldrft", "-l"
            help = "Drift length in cm"
            arg_type = Float64
            default = 100.0
        "--voxel_scale", "-v"
            help = "Voxel scale factor"
            arg_type = Float64
            default = 1.0
        "--voxel_distance_scale", "-d"
            help = "Voxel distance scale factor"
            arg_type = Float64
            default = 1.5
    end
    args = parse_args(s)
    main(; ievent=args["ievent"],
           levent=args["levent"],
           ldrft=args["ldrft"],
           voxel_scale=args["voxel_scale"],
           voxel_distance_scale=args["voxel_distance_scale"])
end
