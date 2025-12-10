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

# Import Petit module functions
import .Petit

cmdir=joinpath(ENV["DATA"], "HD5t/itaca")


function load_data(input_file)
	input_path = joinpath(cmdir, input_file)
	dfs = Petit.get_dataset_dfs(input_path)
	hitsdf = dfs["hits"]

    # Count events from loaded data
    println("Counting events from loaded data...")
    ntot = length(unique(hitsdf.event_id))
    println("Number of events with hits: $ntot")
	return hitsdf
end


function transform_hits_df(df::DataFrame; energy_to_electrons::Float64=1e5/2.5)
      df2 = select(df, Not([:time, :label, :particle_id, :hit_id]))
      df2.electrons = round.(Int, df2.energy .* energy_to_electrons)
      return df2
  end


function continue_on_enter(plt)
        display(plt)
        println("Press Enter to continue...")
        readline()
    end
#=============================================================================
# Main Script
=============================================================================#

f = 1e+5/2.5 # ions per MeV
fkeV = f*1e-3 # ions per keV
	
nbins = 100
nsigma = 3.0 
tK = 297.0 # K
edrift = 500.0  #V/cm

σl_ion  = 0.0
energy_threshold_ions =  10.0

function print_diffusion_params(dp::Petit.DiffusionParams)
    println("DiffusionParams:")
    println("  ldrift           = $(dp.ldrift) cm")
    println("  σ_t              = $(round(dp.sigma_t, digits=3)) mm")
    println("  σ_l              = $(round(dp.sigma_l, digits=3)) mm")
    println("  voxel_size       = $(round(dp.voxel_size, digits=3)) mm")
    println("  max_distance     = $(round(dp.max_distance, digits=3)) mm")
    println("  energy_threshold = $(round(dp.energy_threshold, digits=3)) keV")
    println("  nbins_df         = $(dp.nbins_df)")
    println("  nsigma_df        = $(dp.nsigma_df)")
end

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

    σt_ion_p15 = Petit.sigma_t_ion_mm(tK, ldrft, edrift)
    ion_voxel = σt_ion_p15 * voxel_scale
    max_distance = ion_voxel * voxel_distance_scale

    energy_threshold_keV = energy_threshold_ions/fkeV
    ion_dfpars = Petit.DiffusionParams(ldrft, σt_ion_p15, σl_ion,
                                       ion_voxel, max_distance, energy_threshold_keV,
                                       nbins, nsigma)
    print_diffusion_params(ion_dfpars)

    Rb = 10.0 #mm

    bbdf = load_data("bb0nu/bb0nu_15bar_p1.h5")

    for evt in ievent:levent
        process_event(evt, bbdf, ion_dfpars, σt_ion_p15, ion_voxel,
                      max_distance, energy_threshold_keV, Rb)
    end
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
