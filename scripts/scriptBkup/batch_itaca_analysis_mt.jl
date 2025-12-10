#!/usr/bin/env julia

"""
Multi-threaded event range processing script for ITACA analysis.

This script processes events from ievt to levt using multiple threads,
simulating ion and electron diffusion with ITACA model.
Each thread writes to separate output files for ions and electrons.
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using HDF5
using DataFrames
using Statistics
using Graphs
using CSV

# Load Petit module (itaca_functions.jl is included in Petit.jl)
include(joinpath(pdir, "src", "Petit.jl"))

# Import Petit module functions
import Petit: event_loop_itaca_mt, get_dataset_dfs, Tracks
import Petit: nof_events, count_events
import Petit: save_tracks_to_hdf5, read_tracks_from_hdf5

"""
    parse_bool(s::String)

Parse a string to boolean value.
"""
function parse_bool(s::String)
    s_lower = lowercase(strip(s))
    if s_lower in ["true", "t", "yes", "y", "1"]
        return true
    elseif s_lower in ["false", "f", "no", "n", "0"]
        return false
    else
        error("Cannot parse '$s' as boolean")
    end
end

"""
    save_run_metadata(filepath; kwargs...)

Save run configuration as CSV file for later analysis.
"""
function save_run_metadata(filepath::String;
                           cmdir::String,
                           input_file::String,
                           output_base::String,
                           ievt::Int,
                           levt::Int,
                           nthreads::Int,
                           lmin::Float64,
                           lmax::Float64,
                           lbuff::Float64,
                           pbar::Float64,
                           Dt::Float64,
                           Dl::Float64,
                           nbins_df::Int,
                           nsigma_df::Float64,
                           voxel_scale::Float64,
                           voxel_dd::Float64,
                           energy_threshold_kev::Float64,
                           emin::Float64,
                           emax::Float64)

    df = DataFrame(
        parameter = [
            "directory",
            "input_file",
            "output_base",
            "first_event",
            "last_event",
            "nthreads",
            "lmin_mm",
            "lmax_mm",
            "lbuff_mm",
            "pbar_bar",
            "Dt_mm",
            "Dl_mm",
            "nbins_df",
            "nsigma_df",
            "voxel_scale",
            "voxel_dd",
            "energy_threshold_kev",
            "emin_kev",
            "emax_kev"
        ],
        value = [
            cmdir,
            input_file,
            output_base,
            string(ievt),
            string(levt),
            string(nthreads),
            string(lmin),
            string(lmax),
            string(lbuff),
            string(pbar),
            string(Dt),
            string(Dl),
            string(nbins_df),
            string(nsigma_df),
            string(voxel_scale),
            string(voxel_dd),
            string(energy_threshold_kev),
            string(emin),
            string(emax)
        ]
    )

    CSV.write(filepath, df)
    println("Metadata saved to: $filepath")
end

"""
    main()

Main function to run from command line.

Usage:
    julia batch_itaca_analysis_mt.jl <cmdir> <input_file> <output_base> [options]

Required arguments:
    cmdir           Directory containing the input file
    input_file      Name of the HDF5 input file
    output_base     Base name for output files (without .h5 extension)
                    Each thread will create {output_base}_ion_th_{i}.h5 and {output_base}_ele_th_{i}.h5

Optional arguments:
    --ievt              First event to process (default: 1)
    --levt              Last event to process (default: -1, all events)
    --nthreads          Number of threads to use (default: 1)
    --lmin              Minimum drift length in mm (default: 0.0)
    --lmax              Maximum drift length in mm (default: 200.0)
    --lbuff             Buffer distance from lmin/lmax in mm (default: 10.0)
    --pbar              Pressure in bar (default: 15.0)
    --Dt                Transverse diffusion coefficient in mm (default: 1.6)
    --Dl                Longitudinal diffusion coefficient in mm (default: 0.75)
    --nbins-df          Number of bins for diffusion histogram (default: 300)
    --nsigma-df         Number of sigmas for histogram padding (default: 3.0)
    --voxel-scale       Voxel size as multiple of sigma_t (default: 2.0)
    --voxel-dd          Max distance as multiple of sigma_t (default: 3.0)
    --energy-threshold  Energy threshold in keV (default: 10.0)
    --emin              Minimum event energy in keV (default: -Inf)
    --emax              Maximum event energy in keV (default: Inf)

Example:
    # Process events 1000 to 2000 using 4 threads
    julia batch_itaca_analysis_mt.jl /path/to/data/ input.h5 output_base --ievt=1000 --levt=2000 --nthreads=4

    # Process all events with 15 bar pressure using 8 threads
    julia batch_itaca_analysis_mt.jl /path/to/data/ input.h5 output_base --nthreads=8 --pbar=15.0

Note:
    - You must start Julia with multiple threads: julia -t auto or julia -t 8
    - The script will check and cap threads at the system maximum
    - Each thread writes separate files for ions and electrons
"""
function main()
    # Default values
    ievt = 1
    levt = -1
    nthreads = 1

    # ITACA diffusion parameters
    lmin = 0.0              # Minimum drift length (mm)
    lmax = 200.0            # Maximum drift length (mm)
    lbuff = 10.0            # Buffer from lmin/lmax (mm)
    pbar = 15.0             # Pressure (bar)
    Dt = 1.6                # Transverse diffusion coefficient (mm)
    Dl = 0.75               # Longitudinal diffusion coefficient (mm)
    nbins_df = 300          # Number of bins for diffusion histogram
    nsigma_df = 3.0         # Number of sigmas for histogram padding
    voxel_scale = 2.0       # Voxel size as multiple of sigma_t
    voxel_dd = 3.0          # Max distance as multiple of sigma_t
    energy_threshold_kev = 10.0
    emin = -Inf
    emax = Inf



    # Check minimum required arguments
    if length(ARGS) < 3
        println("Error: Missing required arguments")
        println("\nUsage: julia -t <nthreads> batch_itaca_analysis_mt.jl <cmdir> <input_file> <output_base> [options]")
        println("\nRequired arguments:")
        println("  cmdir           Directory containing the input file")
        println("  input_file      Name of the HDF5 input file")
        println("  output_base     Base name for output files (no .h5 extension)")
        println("\nOptional arguments:")
        println("  --ievt=N               First event to process (default: 1)")
        println("  --levt=N               Last event to process (default: -1, all)")
        println("  --nthreads=N           Number of threads to use (default: 1)")
        println("  --lmin=X               Minimum drift length in mm (default: 0.0)")
        println("  --lmax=X               Maximum drift length in mm (default: 200.0)")
        println("  --lbuff=X              Buffer distance from lmin/lmax in mm (default: 10.0)")
        println("  --pbar=X               Pressure in bar (default: 15.0)")
        println("  --Dt=X                 Transverse diffusion coefficient in mm (default: 1.6)")
        println("  --Dl=X                 Longitudinal diffusion coefficient in mm (default: 0.75)")
        println("  --nbins-df=N           Number of bins for diffusion histogram (default: 300)")
        println("  --nsigma-df=X          Number of sigmas for histogram padding (default: 3.0)")
        println("  --voxel-scale=X        Voxel size as multiple of sigma_t (default: 2.0)")
        println("  --voxel-dd=X           Max distance as multiple of sigma_t (default: 3.0)")
        println("  --energy-threshold=X   Energy threshold in keV (default: 10.0)")
        println("  --emin=X               Minimum event energy in keV (default: -Inf)")
        println("  --emax=X               Maximum event energy in keV (default: Inf)")
        println("\nExample:")
        println("  # Start Julia with 4 threads and process events 1000 to 2000")
        println("  julia -t 4 batch_itaca_analysis_mt.jl /data/ input.h5 output --ievt=1000 --levt=2000 --nthreads=4")
        println("\nImportant:")
        println("  - You MUST start Julia with multiple threads: julia -t auto or julia -t <N>")
        println("  - Each thread creates separate files for ions and electrons")
        exit(1)
    end

    # Parse required arguments
    cmdir = ARGS[1]
    input_file = ARGS[2]
    output_base = ARGS[3]

    # Parse optional arguments
    for arg in ARGS[4:end]
        if startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) != 2
                println("Warning: Ignoring malformed argument: $arg")
                continue
            end
            key, value = parts

            try
                if key == "ievt"
                    ievt = parse(Int, value)
                elseif key == "levt"
                    levt = parse(Int, value)
                elseif key == "nthreads"
                    nthreads = parse(Int, value)
                elseif key == "lmin"
                    lmin = parse(Float64, value)
                elseif key == "lmax"
                    lmax = parse(Float64, value)
                elseif key == "lbuff"
                    lbuff = parse(Float64, value)
                elseif key == "pbar"
                    pbar = parse(Float64, value)
                elseif key == "Dt"
                    Dt = parse(Float64, value)
                elseif key == "Dl"
                    Dl = parse(Float64, value)
                elseif key == "nbins-df"
                    nbins_df = parse(Int, value)
                elseif key == "nsigma-df"
                    nsigma_df = parse(Float64, value)
                elseif key == "voxel-scale"
                    voxel_scale = parse(Float64, value)
                elseif key == "voxel-dd"
                    voxel_dd = parse(Float64, value)
                elseif key == "energy-threshold"
                    energy_threshold_kev = parse(Float64, value)
                elseif key == "emin"
                    emin = parse(Float64, value)
                elseif key == "emax"
                    emax = parse(Float64, value)
                else
                    println("Warning: Unknown argument: --$key")
                end
            catch e
                println("Error parsing argument --$key=$value: $e")
                exit(1)
            end
        end
    end

    # Display configuration
    println("="^60)
    println("MULTI-THREADED ITACA ANALYSIS")
    println("="^60)
    println("Configuration:")
    println("  Directory:          $cmdir")
    println("  Input file:         $input_file")
    println("  Output base:        $output_base")
    println("  First event:        $ievt")
    println("  Last event:         $(levt < 0 ? "all" : levt)")
    println("  Threads requested:  $nthreads")
    println("  Threads available:  $(Threads.nthreads())")
    println("\nITACA parameters:")
    println("  Drift range:        [$lmin, $lmax] mm (buffer: $lbuff mm)")
    println("  Pressure:           $pbar bar")
    println("  Diffusion (Dt/Dl):  $Dt / $Dl mm")
    println("  Histogram bins:     $nbins_df (nsigma: $nsigma_df)")
    println("  Voxel scale:        $voxel_scale x sigma_t")
    println("  Max distance:       $voxel_dd x sigma_t")
    println("  Energy threshold:   $energy_threshold_kev keV")
    println("  Energy range:       [$emin, $emax] keV")
    println("="^60)

    # Save metadata to CSV
    metadata_file = joinpath(cmdir, "$(output_base)_metadata.csv")
    save_run_metadata(metadata_file;
                      cmdir=cmdir,
                      input_file=input_file,
                      output_base=output_base,
                      ievt=ievt,
                      levt=levt,
                      nthreads=nthreads,
                      lmin=lmin,
                      lmax=lmax,
                      lbuff=lbuff,
                      pbar=pbar,
                      Dt=Dt,
                      Dl=Dl,
                      nbins_df=nbins_df,
                      nsigma_df=nsigma_df,
                      voxel_scale=voxel_scale,
                      voxel_dd=voxel_dd,
                      energy_threshold_kev=energy_threshold_kev,
                      emin=emin,
                      emax=emax)

    # Check if Julia was started with enough threads
    if Threads.nthreads() < nthreads
        println("\nWarning: Julia started with only $(Threads.nthreads()) threads")
        println("Restart Julia with: julia -t $nthreads <script>")
        println("Proceeding with $(Threads.nthreads()) threads...")
    end

    # Run multi-threaded processing
    println("\nStarting multi-threaded ITACA event processing...")
    result = Petit.event_loop_itaca_mt(cmdir, output_base;
                             input_file=input_file,
                             ievt=ievt,
                             levt=levt,
                             nthreads=nthreads,
                             lmin=lmin,
                             lmax=lmax,
                             lbuff=lbuff,
                             pbar=pbar,
                             Dt=Dt,
                             Dl=Dl,
                             nbins_df=nbins_df,
                             nsigma_df=nsigma_df,
                             voxel_scale=voxel_scale,
                             voxel_dd=voxel_dd,
                             energy_threshold_kev=energy_threshold_kev,
                             emin=emin,
                             emax=emax)
                   

    println("\n" * "="^60)
    println("PROCESSING COMPLETE")
    println("="^60)
    println("\nOutput files created:")
    println("  MC tracks:")
    for thread_result in result["thread_results"]
        println("    $(thread_result["mc_output_file"])")
    end
    println("  Ion tracks:")
    for thread_result in result["thread_results"]
        println("    $(thread_result["ion_output_file"])")
    end
    println("  Electron tracks:")
    for thread_result in result["thread_results"]
        println("    $(thread_result["ele_output_file"])")
    end

    println("\n" * "="^60)
    println("DONE")
    println("="^60)
end

# Run main if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

