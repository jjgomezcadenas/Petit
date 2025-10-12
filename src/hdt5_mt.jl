"""
Multi-threaded analysis functions for HD5t event processing.
"""

using DataFrames
using HDF5
using ProgressMeter
using Base.Threads

"""
    analysis_loop_single_track_mt(hitsdf::DataFrame, thread_id::Int;
                                  events_to_run=nothing,
                                  initial_event=1,
                                  show_progress=false,
                                  voxel_size_mm::Float64=2.0,
                                  max_distance_mm::Float64=5.0,
                                  energy_threshold_kev::Float64=1.0,
                                  emin::Float64=-Inf,
                                  emax::Float64=Inf)

Multi-threaded version of analysis_loop_single_track2.
Processes events assigned to a specific thread.

# Arguments
- `hitsdf::DataFrame`: DataFrame containing hit data
- `thread_id::Int`: Thread identifier for this worker
- `events_to_run`: Number of events to process (integer)
- `initial_event::Int=1`: First event to process (skips events before this)
- `show_progress::Bool=false`: If true, display a progress bar
- `voxel_size_mm::Float64=2.0`: Voxel size in mm
- `max_distance_mm::Float64=5.0`: Maximum distance for track clustering
- `energy_threshold_kev::Float64=1.0`: Energy threshold in keV
- `emin::Float64=-Inf`: Minimum event energy in keV
- `emax::Float64=Inf`: Maximum event energy in keV

# Returns
- `Vector{Tracks}`: Vector containing only single-track events
"""
function analysis_loop_single_track_mt(hitsdf::DataFrame, thread_id::Int;
                      events_to_run=nothing,
                      initial_event::Int=1,
                      show_progress::Bool=false,
                      voxel_size_mm::Float64=2.0,
                      max_distance_mm::Float64=5.0,
                      energy_threshold_kev::Float64=1.0,
                      emin::Float64=-Inf,
                      emax::Float64=Inf)

    # Get the actual event IDs from the DataFrame
    unique_event_ids = sort(unique(hitsdf.event_id))

    # Validate initial_event
    if initial_event < 1
        throw(ArgumentError("initial_event must be >= 1, got $initial_event"))
    end
    if initial_event > length(unique_event_ids)
        throw(ArgumentError("initial_event ($initial_event) exceeds total events ($(length(unique_event_ids)))"))
    end

    # Determine which events to process
    if events_to_run === nothing
        # Process from initial_event to end
        event_ids_to_process = unique_event_ids[initial_event:end]
    else
        # Process events_to_run events starting from initial_event
        end_idx = min(initial_event + events_to_run - 1, length(unique_event_ids))
        event_ids_to_process = unique_event_ids[initial_event:end_idx]
    end

    # Initialize counters
    n_single_track = 0
    n_events_processed = 0
    n_failed = 0

    TRACKS = Tracks[]

    # Create progress bar if requested
    if show_progress
        prog = Progress(length(event_ids_to_process), 1, "Thread $thread_id processing events: ")
    end

    for nevent in event_ids_to_process
        n_events_processed += 1

        try
            tracks = select_events(hitsdf, nevent;
                                 voxel_size_mm=voxel_size_mm,
                                 max_distance_mm=max_distance_mm,
                                 energy_threshold_kev=energy_threshold_kev,
                                 emin=emin,
                                 emax=emax)

            if length(tracks) == 1
                push!(TRACKS, tracks[1])
                n_single_track += 1
            end

        catch e
            println("Thread $thread_id: Warning: Error processing event $nevent: $e")
            n_failed += 1
        end

        # Update progress bar if enabled
        if show_progress
            next!(prog)
        end
    end

    # Finish progress bar if enabled
    if show_progress
        finish!(prog)
    end

    println("Thread $thread_id: Analysis completed! Processed $n_events_processed events.")
    println("Thread $thread_id: Number of single track events $n_single_track")

    return TRACKS
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
    event_loop_single_track_mt(cmdir::String, output_base::String;
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

Multi-threaded event processing. Loads data once and distributes events among threads.
Each thread writes to a separate output file named {output_base}_th_{i}.h5

# Arguments
- `cmdir::String`: Directory containing the input file
- `output_base::String`: Base name for output files (without extension)
- `input_file::String="0nubb.next.h5"`: Name of the HDF5 input file
- `ievt::Int=1`: First event to process
- `levt::Int=-1`: Last event to process (-1 means all events)
- `nthreads::Int=1`: Number of threads to use
- `voxel_size_mm::Float64=5.0`: Voxel size in mm
- `max_distance_mm::Float64=10.0`: Maximum distance for track clustering
- `energy_threshold_kev::Float64=10.0`: Energy threshold in keV
- `emin::Float64=-Inf`: Minimum event energy in keV
- `emax::Float64=Inf`: Maximum event energy in keV
- `xyc::Float64=1950.0`: XY fiducial cut radius
- `zc::Float64=10.0`: Z fiducial cut distance

# Returns
- `Dict`: Summary statistics including per-thread results
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

    input_path = joinpath(cmdir, input_file)

    # Get number of events from configuration
    nevents_config = nof_events(input_path)
    println("Number of events from config: $nevents_config")

    # Load data once (shared across all threads)
    println("Loading data from: $input_path")
    dfs = get_dataset_dfs(input_path)
    hitsdf = dfs["hits"]

    # Count events from loaded data
    println("Counting events from loaded data...")
    ntot = length(unique(hitsdf.event_id))
    println("Number of events with hits: $ntot")

    if nevents_config != ntot
        println("Warning: Config events ($nevents_config) != events with hits ($ntot)")
    end

    # Validate event range
    if ievt < 1
        error("ievt must be >= 1, got $ievt")
    end
    if ievt > ntot
        error("ievt ($ievt) exceeds total events ($ntot)")
    end

    # Determine last event to process
    last_event = levt < 0 ? ntot : min(levt, ntot)

    if last_event < ievt
        error("levt ($last_event) must be >= ievt ($ievt)")
    end

    # Calculate number of events to process
    nevents_to_process = last_event - ievt + 1

    # Determine optimal number of threads
    optimal_nthreads = get_optimal_threads(nthreads)

    println("\nEvent processing configuration:")
    println("  First event:        $ievt")
    println("  Last event:         $last_event")
    println("  Events to skip:     $(ievt - 1)")
    println("  Events to process:  $nevents_to_process")
    println("  Threads to use:     $optimal_nthreads")

    # Split events among threads
    thread_ranges = split_events_for_threads(nevents_to_process, optimal_nthreads, ievt)

    println("\nThread event distribution:")
    for (i, (start_evt, num_evts)) in enumerate(thread_ranges)
        println("  Thread $i: events $start_evt to $(start_evt + num_evts - 1) ($num_evts events)")
    end

    # Process events in parallel
    println("\n" * "="^60)
    println("STARTING MULTI-THREADED PROCESSING")
    println("="^60)

    results = Vector{Any}(undef, length(thread_ranges))

    Threads.@threads for i in 1:length(thread_ranges)
        start_evt, num_evts = thread_ranges[i]

        # Process events for this thread
        tracks = analysis_loop_single_track_mt(hitsdf, i;
                                              events_to_run=num_evts,
                                              initial_event=start_evt,
                                              show_progress=false,  # Disable for cleaner output
                                              voxel_size_mm=voxel_size_mm,
                                              max_distance_mm=max_distance_mm,
                                              energy_threshold_kev=energy_threshold_kev,
                                              emin=emin,
                                              emax=emax)

        # Save tracks to thread-specific file
        output_file = "$(output_base)_th_$i.h5"
        output_path = joinpath(cmdir, output_file)

        h5open(output_path, "w") do fid
            # Store metadata
            attrs(fid)["input_file"] = input_file
            attrs(fid)["nevents_from_config"] = nevents_config
            attrs(fid)["total_events_in_file"] = ntot
            attrs(fid)["thread_id"] = i
            attrs(fid)["first_event_processed"] = start_evt
            attrs(fid)["last_event_processed"] = start_evt + num_evts - 1
            attrs(fid)["events_skipped"] = start_evt - 1
            attrs(fid)["events_processed"] = num_evts
            attrs(fid)["voxel_size_mm"] = voxel_size_mm
            attrs(fid)["max_distance_mm"] = max_distance_mm
            attrs(fid)["energy_threshold_kev"] = energy_threshold_kev
            attrs(fid)["emin"] = emin
            attrs(fid)["emax"] = emax
            attrs(fid)["xyc"] = xyc
            attrs(fid)["zc"] = zc
            attrs(fid)["total_tracks_saved"] = length(tracks)

            # Save tracks
            if !isempty(tracks)
                save_tracks_to_hdf5(tracks, fid, 1)
            end
        end

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

    println("\n" * "="^60)
    println("MULTI-THREADED PROCESSING COMPLETE")
    println("="^60)

    # Summarize results
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
