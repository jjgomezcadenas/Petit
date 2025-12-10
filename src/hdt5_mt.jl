"""
Multi-threaded analysis functions for HD5t event processing.
"""

using DataFrames
using HDF5
using Base.Threads


"""
    get_events_to_process(hitsdf, events_to_run, initial_event)

Get list of event IDs to process based on parameters. Validates inputs.
"""
function get_events_to_process(hitsdf::DataFrame, events_to_run, initial_event::Int)
    unique_event_ids = sort(unique(hitsdf.event_id))

    if initial_event < 1
        throw(ArgumentError("initial_event must be >= 1, got $initial_event"))
    end
    if initial_event > length(unique_event_ids)
        throw(ArgumentError("initial_event ($initial_event) exceeds total events ($(length(unique_event_ids)))"))
    end

    if events_to_run === nothing
        return unique_event_ids[initial_event:end]
    else
        end_idx = min(initial_event + events_to_run - 1, length(unique_event_ids))
        return unique_event_ids[initial_event:end_idx]
    end
end


"""
    analysis_loop_single_track_mt(hitsdf, thread_id; kwargs...)

Multi-threaded single-track analysis. Returns vector of single-track events.
"""
function analysis_loop_single_track_mt(hitsdf::DataFrame, thread_id::Int;
                                       events_to_run=nothing,
                                       initial_event::Int=1,
                                       voxel_size_mm::Float64=2.0,
                                       max_distance_mm::Float64=5.0,
                                       energy_threshold_kev::Float64=1.0,
                                       emin::Float64=-Inf,
                                       emax::Float64=Inf)

    event_ids = get_events_to_process(hitsdf, events_to_run, initial_event)
    n_single_track = 0
    n_events_processed = 0
    n_failed = 0
    TRACKS = Tracks[]

    for nevent in event_ids
        n_events_processed += 1
        try
            tracks = select_events(hitsdf, nevent;
                                   voxel_size_mm=voxel_size_mm,
                                   max_distance_mm=max_distance_mm,
                                   energy_threshold_kev=energy_threshold_kev,
                                   emin=emin, emax=emax)
            if length(tracks) == 1
                push!(TRACKS, tracks[1])
                n_single_track += 1
            end
        catch e
            println("Thread $thread_id: Warning: Error processing event $nevent: $e")
            n_failed += 1
        end
    end

    println("Thread $thread_id: Analysis completed! Processed $n_events_processed events.")
    println("Thread $thread_id: Number of single track events $n_single_track")
    return TRACKS
end


"""
    analysis_loop_itaca_mt(hitsdf, thread_id; kwargs...)

Multi-threaded ITACA analysis. Returns tuple (ion_tracks, ele_tracks).
"""
function analysis_loop_itaca_mt(hitsdf::DataFrame, thread_id::Int;
                                events_to_run=nothing,
                                initial_event::Int=1,
                                lmin::Float64=0.0,
                                lmax::Float64=200.0,
                                lbuff::Float64=10.0,
                                pbar::Float64=15.0,
                                Dt::Float64=1.6,
                                Dl::Float64=0.75,
                                nbins_df::Int=300,
                                nsigma_df::Float64=3.0,
                                voxel_scale::Float64=2.0,
                                voxel_dd::Float64=3.0,
                                energy_threshold_kev::Float64=1.0,
                                emin::Float64=-Inf,
                                emax::Float64=Inf)

    event_ids = get_events_to_process(hitsdf, events_to_run, initial_event)
    n_events_processed = 0
    n_failed = 0
    ION_TRACKS = Tracks[]
    ELE_TRACKS = Tracks[]
    MC_TRACKS = Tracks[]

    for nevent in event_ids
        n_events_processed += 1
        try
            mc_tracks, ion_tracks, ele_tracks = select_events_itaca(hitsdf, nevent;
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
                                                         emin=emin, emax=emax)

            if length(mc_tracks) == 1
                push!(MC_TRACKS, mc_tracks[1])
            end
             if length(ion_tracks) == 1
                push!(ION_TRACKS, ion_tracks[1])
            end
            if length(ele_tracks) == 1
                push!(ELE_TRACKS, ele_tracks[1])
            end
        catch e
            println("Thread $thread_id: Warning: Error processing event $nevent: $e")
            n_failed += 1
        end
    end

    println("Thread $thread_id: Processed $n_events_processed events.")
    println("Thread $thread_id: Ion single-tracks: $(length(ION_TRACKS)), Electron single-tracks: $(length(ELE_TRACKS))")
    return (MC_TRACKS, ION_TRACKS, ELE_TRACKS)
end


"""
    get_optimal_threads(requested_threads::Int)

Determine the optimal number of threads to use based on system capabilities
and user request.

# Arguments
- `requested_threads::Int`: Number of threads requested by user

# Returns
- `Int`: Optimal number of threads to use (capped at system limit)
"""
function get_optimal_threads(requested_threads::Int)
    # Get the number of available CPU threads
    available_threads = Threads.nthreads()

    # Get the number of physical cores (more conservative estimate)
    # This is system-dependent, so we'll use available threads as upper bound
    max_recommended = available_threads

    if requested_threads > max_recommended
        println("Warning: Requested $requested_threads threads, but only $max_recommended available.")
        println("Setting nthreads to $max_recommended")
        return max_recommended
    elseif requested_threads < 1
        println("Warning: Requested $requested_threads threads (invalid). Setting to 1.")
        return 1
    else
        return requested_threads
    end
end


"""
    split_events_for_threads(total_events::Int, nthreads::Int, initial_event::Int=1)

Split the event range among threads as evenly as possible.

# Arguments
- `total_events::Int`: Total number of events to process
- `nthreads::Int`: Number of threads to use
- `initial_event::Int=1`: First event to start from

# Returns
- `Vector{Tuple{Int, Int}}`: Vector of (start_event, num_events) for each thread
"""
function split_events_for_threads(total_events::Int, nthreads::Int, initial_event::Int=1)
    events_per_thread = div(total_events, nthreads)
    remainder = total_events % nthreads

    thread_ranges = Tuple{Int, Int}[]
    current_start = initial_event

    for i in 1:nthreads
        # Distribute remainder events to first threads
        thread_events = events_per_thread + (i <= remainder ? 1 : 0)

        if thread_events > 0
            push!(thread_ranges, (current_start, thread_events))
            current_start += thread_events
        end
    end

    return thread_ranges
end


"""
    load_and_validate_input(cmdir, input_file, ievt, levt)

Load HDF5 data and validate event range. Returns (hitsdf, ntot, nevents_config, last_event).
"""
function load_and_validate_input(cmdir::String, input_file::String, ievt::Int, levt::Int)
    input_path = joinpath(cmdir, input_file)

    nevents_config = nof_events(input_path)
    println("Number of events from config: $nevents_config")

    println("Loading data from: $input_path")
    dfs = get_dataset_dfs(input_path)
    hitsdf = dfs["hits"]

    println("Counting events from loaded data...")
    ntot = length(unique(hitsdf.event_id))
    println("Number of events with hits: $ntot")

    if nevents_config != ntot
        println("Warning: Config events ($nevents_config) != events with hits ($ntot)")
    end

    if ievt < 1
        error("ievt must be >= 1, got $ievt")
    end
    if ievt > ntot
        error("ievt ($ievt) exceeds total events ($ntot)")
    end

    last_event = levt < 0 ? ntot : min(levt, ntot)
    if last_event < ievt
        error("levt ($last_event) must be >= ievt ($ievt)")
    end

    return hitsdf, ntot, nevents_config, last_event
end


"""
    print_processing_config(ievt, last_event, optimal_nthreads, thread_ranges)

Print event processing configuration and thread distribution.
"""
function print_processing_config(ievt::Int, last_event::Int, optimal_nthreads::Int,
                                 thread_ranges::Vector{Tuple{Int,Int}})
    nevents_to_process = last_event - ievt + 1
    println("\nEvent processing configuration:")
    println("  First event:        $ievt")
    println("  Last event:         $last_event")
    println("  Events to skip:     $(ievt - 1)")
    println("  Events to process:  $nevents_to_process")
    println("  Threads to use:     $optimal_nthreads")

    println("\nThread event distribution:")
    for (i, (start_evt, num_evts)) in enumerate(thread_ranges)
        println("  Thread $i: events $start_evt to $(start_evt + num_evts - 1) ($num_evts events)")
    end
end


"""
    save_thread_results(output_path, tracks, metadata)

Save tracks and metadata to HDF5 file.
"""
function save_thread_results(output_path::String, tracks::Vector{Tracks}, metadata::Dict)
    h5open(output_path, "w") do fid
        for (key, val) in metadata
            attrs(fid)[key] = val
        end
        attrs(fid)["total_tracks_saved"] = length(tracks)
        if !isempty(tracks)
            save_tracks_to_hdf5(tracks, fid, 1)
        end
    end
end


"""
    save_thread_results_itaca(output_path, mc_tracks, ion_tracks, ele_tracks, metadata)

Save MC, ion and electron tracks to separate HDF5 files.
"""
function save_thread_results_itaca(output_base::String, mc_tracks::Vector{Tracks},
                                   ion_tracks::Vector{Tracks}, ele_tracks::Vector{Tracks},
                                   metadata::Dict)
    # Save MC tracks
    mc_path = output_base * "_mc.h5"
    h5open(mc_path, "w") do fid
        for (key, val) in metadata
            attrs(fid)[key] = val
        end
        attrs(fid)["track_type"] = "mc"
        attrs(fid)["total_tracks_saved"] = length(mc_tracks)
        if !isempty(mc_tracks)
            save_tracks_to_hdf5(mc_tracks, fid, 1)
        end
    end

    # Save ion tracks
    ion_path = output_base * "_ion.h5"
    h5open(ion_path, "w") do fid
        for (key, val) in metadata
            attrs(fid)[key] = val
        end
        attrs(fid)["track_type"] = "ion"
        attrs(fid)["total_tracks_saved"] = length(ion_tracks)
        if !isempty(ion_tracks)
            save_tracks_to_hdf5(ion_tracks, fid, 1)
        end
    end

    # Save electron tracks
    ele_path = output_base * "_ele.h5"
    h5open(ele_path, "w") do fid
        for (key, val) in metadata
            attrs(fid)[key] = val
        end
        attrs(fid)["track_type"] = "electron"
        attrs(fid)["total_tracks_saved"] = length(ele_tracks)
        if !isempty(ele_tracks)
            save_tracks_to_hdf5(ele_tracks, fid, 1)
        end
    end

    return (mc_path, ion_path, ele_path)
end


"""
    print_and_return_summary(thread_ranges, results, nevents_to_process, ntot, output_base)

Print summary and return results dictionary.
"""
function print_and_return_summary(thread_ranges, results, nevents_to_process, ntot, output_base)
    println("\n" * "="^60)
    println("MULTI-THREADED PROCESSING COMPLETE")
    println("="^60)

    total_tracks = sum(r["tracks_saved"] for r in results)
    println("\nSummary:")
    println("  Total threads used:     $(length(thread_ranges))")
    println("  Total events processed: $nevents_to_process")
    println("  Total tracks saved:     $total_tracks")
    println("\nPer-thread results:")
    for r in results
        println("  Thread $(r["thread_id"]): $(r["tracks_saved"]) tracks, events $(r["first_event"])-$(r["last_event"])")
    end

    return Dict(
        "total_events" => ntot,
        "events_processed" => nevents_to_process,
        "threads_used" => length(thread_ranges),
        "total_tracks_saved" => total_tracks,
        "thread_results" => results,
        "output_base" => output_base
    )
end


"""
    print_and_return_summary_itaca(thread_ranges, results, nevents_to_process, ntot, output_base)

Print summary for ITACA (mc + ion + electron tracks) and return results dictionary.
"""
function print_and_return_summary_itaca(thread_ranges, results, nevents_to_process, ntot, output_base)
    println("\n" * "="^60)
    println("MULTI-THREADED ITACA PROCESSING COMPLETE")
    println("="^60)

    total_mc = sum(r["mc_tracks_saved"] for r in results)
    total_ion = sum(r["ion_tracks_saved"] for r in results)
    total_ele = sum(r["ele_tracks_saved"] for r in results)
    println("\nSummary:")
    println("  Total threads used:       $(length(thread_ranges))")
    println("  Total events processed:   $nevents_to_process")
    println("  Total MC tracks saved:    $total_mc")
    println("  Total ion tracks saved:   $total_ion")
    println("  Total electron tracks:    $total_ele")
    println("\nPer-thread results:")
    for r in results
        println("  Thread $(r["thread_id"]): $(r["mc_tracks_saved"]) mc, $(r["ion_tracks_saved"]) ion, $(r["ele_tracks_saved"]) electron, events $(r["first_event"])-$(r["last_event"])")
    end

    return Dict(
        "total_events" => ntot,
        "events_processed" => nevents_to_process,
        "threads_used" => length(thread_ranges),
        "total_mc_tracks" => total_mc,
        "total_ion_tracks" => total_ion,
        "total_ele_tracks" => total_ele,
        "thread_results" => results,
        "output_base" => output_base
    )
end


"""
    event_loop_single_track_mt(cmdir, output_base; kwargs...)

Multi-threaded event processing. Loads data once and distributes events among threads.
Each thread writes to a separate output file named {output_base}_th_{i}.h5
"""
function event_loop_single_track_mt(cmdir::String, output_base::String;
                                    input_file::String="0nubb.next.h5",
                                    ievt::Int=1,
                                    levt::Int=-1,
                                    nthreads::Int=1,
                                    voxel_size_mm::Float64=5.0,
                                    max_distance_mm::Float64=10.0,
                                    energy_threshold_kev::Float64=10.0,
                                    emin::Float64=-Inf,
                                    emax::Float64=Inf,
                                    xyc::Float64=1950.0,
                                    zc::Float64=10.0)

    hitsdf, ntot, nevents_config, last_event = load_and_validate_input(cmdir, input_file, ievt, levt)
    nevents_to_process = last_event - ievt + 1
    optimal_nthreads = get_optimal_threads(nthreads)
    thread_ranges = split_events_for_threads(nevents_to_process, optimal_nthreads, ievt)

    print_processing_config(ievt, last_event, optimal_nthreads, thread_ranges)

    println("\n" * "="^60)
    println("STARTING MULTI-THREADED PROCESSING")
    println("="^60)

    results = Vector{Any}(undef, length(thread_ranges))

    Threads.@threads for i in 1:length(thread_ranges)
        start_evt, num_evts = thread_ranges[i]

        tracks = analysis_loop_single_track_mt(hitsdf, i;
                                               events_to_run=num_evts,
                                               initial_event=start_evt,
                                               voxel_size_mm=voxel_size_mm,
                                               max_distance_mm=max_distance_mm,
                                               energy_threshold_kev=energy_threshold_kev,
                                               emin=emin, emax=emax)

        output_file = "$(output_base)_th_$i.h5"
        output_path = joinpath(cmdir, output_file)

        metadata = Dict(
            "input_file" => input_file,
            "nevents_from_config" => nevents_config,
            "total_events_in_file" => ntot,
            "thread_id" => i,
            "first_event_processed" => start_evt,
            "last_event_processed" => start_evt + num_evts - 1,
            "events_skipped" => start_evt - 1,
            "events_processed" => num_evts,
            "voxel_size_mm" => voxel_size_mm,
            "max_distance_mm" => max_distance_mm,
            "energy_threshold_kev" => energy_threshold_kev,
            "emin" => emin, "emax" => emax,
            "xyc" => xyc, "zc" => zc
        )
        save_thread_results(output_path, tracks, metadata)

        results[i] = Dict(
            "thread_id" => i,
            "output_file" => output_path,
            "events_processed" => num_evts,
            "tracks_saved" => length(tracks),
            "first_event" => start_evt,
            "last_event" => start_evt + num_evts - 1
        )
        println("Thread $i: Saved $(length(tracks)) tracks to $output_file")
    end

    return print_and_return_summary(thread_ranges, results, nevents_to_process, ntot, output_base)
end


"""
    event_loop_itaca_mt(cmdir, output_base; kwargs...)

Multi-threaded ITACA event processing with diffusion simulation.
"""
function event_loop_itaca_mt(cmdir::String, output_base::String;
                             input_file::String="0nubb.next.h5",
                             ievt::Int=1,
                             levt::Int=-1,
                             nthreads::Int=1,
                             lmin::Float64=0.0,
                             lmax::Float64=200.0,
                             lbuff::Float64=10.0,
                             pbar::Float64=15.0,
                             Dt::Float64=1.6,
                             Dl::Float64=0.75,
                             nbins_df::Int=300,
                             nsigma_df::Float64=3.0,
                             voxel_scale::Float64=2.0,
                             voxel_dd::Float64=3.0,
                             energy_threshold_kev::Float64=10.0,
                             emin::Float64=-Inf,
                             emax::Float64=Inf)

    hitsdf, ntot, nevents_config, last_event = load_and_validate_input(cmdir, input_file, ievt, levt)
    nevents_to_process = last_event - ievt + 1
    optimal_nthreads = get_optimal_threads(nthreads)
    thread_ranges = split_events_for_threads(nevents_to_process, optimal_nthreads, ievt)

    print_processing_config(ievt, last_event, optimal_nthreads, thread_ranges)

    println("\n" * "="^60)
    println("STARTING MULTI-THREADED PROCESSING")
    println("="^60)

    results = Vector{Any}(undef, length(thread_ranges))

    Threads.@threads for i in 1:length(thread_ranges)
        start_evt, num_evts = thread_ranges[i]

        mc_tracks, ion_tracks, ele_tracks = analysis_loop_itaca_mt(hitsdf, i;
                                                        events_to_run=num_evts,
                                                        initial_event=start_evt,
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
                                                        emin=emin, emax=emax)

        output_base_th = joinpath(cmdir, "$(output_base)_th_$i")

        metadata = Dict(
            "input_file" => input_file,
            "nevents_from_config" => nevents_config,
            "total_events_in_file" => ntot,
            "thread_id" => i,
            "first_event_processed" => start_evt,
            "last_event_processed" => start_evt + num_evts - 1,
            "events_skipped" => start_evt - 1,
            "events_processed" => num_evts,
            "lmin" => lmin, "lmax" => lmax, "lbuff" => lbuff,
            "pbar" => pbar, "Dt" => Dt, "Dl" => Dl,
            "nbins_df" => nbins_df, "nsigma_df" => nsigma_df,
            "voxel_scale" => voxel_scale, "voxel_dd" => voxel_dd,
            "energy_threshold_kev" => energy_threshold_kev,
            "emin" => emin, "emax" => emax
        )
        mc_path, ion_path, ele_path = save_thread_results_itaca(output_base_th, mc_tracks, ion_tracks, ele_tracks, metadata)

        results[i] = Dict(
            "thread_id" => i,
            "mc_output_file" => mc_path,
            "ion_output_file" => ion_path,
            "ele_output_file" => ele_path,
            "events_processed" => num_evts,
            "mc_tracks_saved" => length(mc_tracks),
            "ion_tracks_saved" => length(ion_tracks),
            "ele_tracks_saved" => length(ele_tracks),
            "first_event" => start_evt,
            "last_event" => start_evt + num_evts - 1
        )
        println("Thread $i: Saved $(length(mc_tracks)) mc, $(length(ion_tracks)) ion, $(length(ele_tracks)) electron tracks")
    end

    return print_and_return_summary_itaca(thread_ranges, results, nevents_to_process, ntot, output_base)
end
