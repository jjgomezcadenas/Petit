struct AnalysisResults
    single_track::DataFrame
    two_track_primary::DataFrame
    two_track_secondary::DataFrame
    three_track_primary::DataFrame
    three_track_secondary::DataFrame
    n_events_processed::Int
    n_single_track::Int
    n_two_track::Int
    n_three_plus_track::Int
    n_failed::Int
end

function determine_events_to_process(events_to_run, unique_event_ids)
    # Determine which events to process
    if events_to_run === nothing
        # Process all events
        event_ids_to_process = unique_event_ids
    elseif isa(events_to_run, Integer)
        # Process first N events
        n_available = length(unique_event_ids)
        n_to_process = min(events_to_run, n_available)
        event_ids_to_process = unique_event_ids[1:n_to_process]
    elseif isa(events_to_run, AbstractRange) || isa(events_to_run, AbstractVector)
        # Process specific event IDs if they exist in the DataFrame
        event_ids_to_process = filter(id -> id in unique_event_ids, collect(events_to_run))
        if isempty(event_ids_to_process) && !isempty(events_to_run)
            # Only warn if we requested specific IDs that don't exist (not for empty input)
            println("Warning: None of the requested event IDs exist in the DataFrame")
            if !isempty(unique_event_ids)
                println("Available event IDs range from $(minimum(unique_event_ids)) to $(maximum(unique_event_ids))")
            else
                println("No events available in the DataFrame")
            end
        end
    else
        throw(ArgumentError("events_to_run must be nothing, an Integer, a Range, or a Vector"))
    end
    return event_ids_to_process
end

"""
    analysis_loop(hitsdf::DataFrame; events_to_run=nothing, voxel_size_mm=2.0, max_distance_mm=5.0, energy_threshold_kev=1.0)

Analyze multiple events to extract track energy distributions.

# Arguments
- `hitsdf::DataFrame`: DataFrame containing hit data
- `events_to_run`: Either:
  - `nothing` (default): Process all unique event IDs in the DataFrame
  - An integer: Process the first N unique event IDs
  - A range (e.g., 1:100): Process event IDs in this range if they exist
  - A vector of specific event IDs to process
- `voxel_size_mm::Float64=2.0`: Voxel size in mm
- `max_distance_mm::Float64=5.0`: Maximum distance for track clustering
- `energy_threshold_kev::Float64=1.0`: Energy threshold in keV

# Returns
- `AnalysisResults`: Struct containing energy arrays and statistics
"""
function analysis_loop(hitsdf::DataFrame; 
                      events_to_run=nothing, 
                      voxel_size_mm::Float64=2.0,
                      max_distance_mm::Float64=5.0, 
                      energy_threshold_kev::Float64=1.0)
    
    # Get the actual event IDs from the DataFrame
    unique_event_ids = sort(unique(hitsdf.event_id))

    # Determine which events to process using helper function
    event_ids_to_process = determine_events_to_process(events_to_run, unique_event_ids)
    
    # Initialize arrays to collect data for DataFrames
    single_track_data = (event_id=Int[], energy=Float64[], x=Float64[], y=Float64[], z=Float64[])
    two_track_primary_data = (event_id=Int[], energy=Float64[], x=Float64[], y=Float64[], z=Float64[])
    two_track_secondary_data = (event_id=Int[], energy=Float64[], x=Float64[], y=Float64[], z=Float64[])
    three_track_primary_data = (event_id=Int[], energy=Float64[], x=Float64[], y=Float64[], z=Float64[])
    three_track_secondary_data = (event_id=Int[], energy=Float64[], x=Float64[], y=Float64[], z=Float64[])
    
    # Initialize counters
    n_single_track = 0
    n_two_track = 0
    n_three_plus_track = 0
    n_failed = 0
    n_events_processed = 0
    
    # Create progress bar
    total_events = length(event_ids_to_process)
    if total_events == 0
        println("No events to process")
        # Return empty DataFrames
        empty_df = DataFrame(event_id=Int[], energy=Float64[], x=Float64[], y=Float64[], z=Float64[])
        return AnalysisResults(
            empty_df, empty_df, empty_df, empty_df, empty_df,
            0, 0, 0, 0, 0
        )
    end
    
    progress = Progress(total_events, dt=0.5, 
                       desc="Processing events: ",
                       barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|'),
                       barlen=50)
    
    for nevent in event_ids_to_process
        n_events_processed += 1
        
        try
            tracks = select_events(hitsdf, nevent; 
                                 voxel_size_mm=voxel_size_mm, 
                                 max_distance_mm=max_distance_mm, 
                                 energy_threshold_kev=energy_threshold_kev)

            if length(tracks) == 1
                # Single track event
                energy_kev = 1e+3 * sum(tracks[1].voxels.energy)
                # Add data for each voxel in the track
                for i in 1:nrow(tracks[1].voxels)
                    push!(single_track_data.event_id, nevent)
                    push!(single_track_data.energy, energy_kev)
                    push!(single_track_data.x, tracks[1].voxels.x[i])
                    push!(single_track_data.y, tracks[1].voxels.y[i])
                    push!(single_track_data.z, tracks[1].voxels.z[i])
                end
                n_single_track += 1
                
            elseif length(tracks) == 2 
                # Two track event
                primary_energy = 1e+3 * sum(tracks[1].voxels.energy)
                secondary_energy = 1e+3 * sum(tracks[2].voxels.energy)
                
                # Add primary track data
                for i in 1:nrow(tracks[1].voxels)
                    push!(two_track_primary_data.event_id, nevent)
                    push!(two_track_primary_data.energy, primary_energy)
                    push!(two_track_primary_data.x, tracks[1].voxels.x[i])
                    push!(two_track_primary_data.y, tracks[1].voxels.y[i])
                    push!(two_track_primary_data.z, tracks[1].voxels.z[i])
                end
                
                # Add secondary track data
                for i in 1:nrow(tracks[2].voxels)
                    push!(two_track_secondary_data.event_id, nevent)
                    push!(two_track_secondary_data.energy, secondary_energy)
                    push!(two_track_secondary_data.x, tracks[2].voxels.x[i])
                    push!(two_track_secondary_data.y, tracks[2].voxels.y[i])
                    push!(two_track_secondary_data.z, tracks[2].voxels.z[i])
                end
                
                n_two_track += 1
                
            elseif length(tracks) >= 3 
                # Three or more track event
                primary_energy = 1e+3 * sum(tracks[1].voxels.energy)
                
                # Add primary track data
                for i in 1:nrow(tracks[1].voxels)
                    push!(three_track_primary_data.event_id, nevent)
                    push!(three_track_primary_data.energy, primary_energy)
                    push!(three_track_primary_data.x, tracks[1].voxels.x[i])
                    push!(three_track_primary_data.y, tracks[1].voxels.y[i])
                    push!(three_track_primary_data.z, tracks[1].voxels.z[i])
                end
                
                # Add secondary tracks data
                for n in 2:length(tracks)
                    secondary_energy = 1e+3 * sum(tracks[n].voxels.energy)
                    for i in 1:nrow(tracks[n].voxels)
                        push!(three_track_secondary_data.event_id, nevent)
                        push!(three_track_secondary_data.energy, secondary_energy)
                        push!(three_track_secondary_data.x, tracks[n].voxels.x[i])
                        push!(three_track_secondary_data.y, tracks[n].voxels.y[i])
                        push!(three_track_secondary_data.z, tracks[n].voxels.z[i])
                    end
                end
                n_three_plus_track += 1
            end
            
        catch e
            println("Warning: Error processing event $nevent: $e")
            n_failed += 1
        finally
            # Update progress bar with current statistics
            next!(progress, showvalues = [
                (:Event, nevent),
                (:Single_tracks, n_single_track),
                (:Two_tracks, n_two_track), 
                (:Three_plus_tracks, n_three_plus_track),
                (:Failed, n_failed)
            ])
        end
    end
    
    # Finish progress bar
    finish!(progress)
    println(" Analysis completed! Processed $n_events_processed events.")
    
    # Create DataFrames from collected data
    single_track_df = DataFrame(
        event_id=single_track_data.event_id,
        energy=single_track_data.energy,
        x=single_track_data.x,
        y=single_track_data.y,
        z=single_track_data.z
    )
    
    two_track_primary_df = DataFrame(
        event_id=two_track_primary_data.event_id,
        energy=two_track_primary_data.energy,
        x=two_track_primary_data.x,
        y=two_track_primary_data.y,
        z=two_track_primary_data.z
    )
    
    two_track_secondary_df = DataFrame(
        event_id=two_track_secondary_data.event_id,
        energy=two_track_secondary_data.energy,
        x=two_track_secondary_data.x,
        y=two_track_secondary_data.y,
        z=two_track_secondary_data.z
    )
    
    three_track_primary_df = DataFrame(
        event_id=three_track_primary_data.event_id,
        energy=three_track_primary_data.energy,
        x=three_track_primary_data.x,
        y=three_track_primary_data.y,
        z=three_track_primary_data.z
    )
    
    three_track_secondary_df = DataFrame(
        event_id=three_track_secondary_data.event_id,
        energy=three_track_secondary_data.energy,
        x=three_track_secondary_data.x,
        y=three_track_secondary_data.y,
        z=three_track_secondary_data.z
    )
    
    return AnalysisResults(
        single_track_df,
        two_track_primary_df,
        two_track_secondary_df,
        three_track_primary_df,
        three_track_secondary_df,
        n_events_processed,
        n_single_track,
        n_two_track,
        n_three_plus_track,
        n_failed
    )

end

"""
    analysis_loop_single_track(hitsdf::DataFrame; events_to_run=nothing, voxel_size_mm=2.0, max_distance_mm=5.0, energy_threshold_kev=1.0)

Analyze multiple events to extract only single-track events and return them as a vector of Tracks objects.

# Arguments
- `hitsdf::DataFrame`: DataFrame containing hit data
- `events_to_run`: Same options as analysis_loop
- `voxel_size_mm::Float64=2.0`: Voxel size in mm
- `max_distance_mm::Float64=5.0`: Maximum distance for track clustering
- `energy_threshold_kev::Float64=1.0`: Energy threshold in keV

# Returns
- `Vector{Tracks}`: Vector containing only single-track events
"""
function analysis_loop_single_track(hitsdf::DataFrame; 
                      events_to_run=nothing, 
                      voxel_size_mm::Float64=2.0,
                      max_distance_mm::Float64=5.0, 
                      energy_threshold_kev::Float64=1.0)
    
    # Get the actual event IDs from the DataFrame
    unique_event_ids = sort(unique(hitsdf.event_id))
    
    event_ids_to_process = determine_events_to_process(events_to_run, unique_event_ids)
    # Initialize counters
    n_single_track = 0
    n_events_processed = 0
    n_failed = 0
    
    TRACKS = Tracks[]
    for nevent in event_ids_to_process
        n_events_processed += 1
        
        try
            tracks = select_events(hitsdf, nevent; 
                                 voxel_size_mm=voxel_size_mm, 
                                 max_distance_mm=max_distance_mm, 
                                 energy_threshold_kev=energy_threshold_kev)

            if length(tracks) == 1
                push!(TRACKS, tracks[1])

                n_single_track += 1
            else
                continue
            end

        catch e
            println("Warning: Error processing event $nevent: $e")
            n_failed += 1
        end
    end
    
    println(" Analysis completed! Processed $n_events_processed events.")
    println(" Number of single track events $n_single_track")
    
    return TRACKS

end

"""
    validate_event_loop_parameters(events_to_run, voxel_size_mm, max_distance_mm, energy_threshold_kev)

Helper function to validate common parameters for event loop functions.
"""
function validate_event_loop_parameters(events_to_run, voxel_size_mm, max_distance_mm, energy_threshold_kev)
    if events_to_run <= 0
        throw(ArgumentError("events_to_run must be positive, got $events_to_run"))
    end
    if voxel_size_mm <= 0
        throw(ArgumentError("voxel_size_mm must be positive, got $voxel_size_mm"))
    end
    if max_distance_mm <= 0
        throw(ArgumentError("max_distance_mm must be positive, got $max_distance_mm"))
    end
    if energy_threshold_kev < 0
        throw(ArgumentError("energy_threshold_kev must be non-negative, got $energy_threshold_kev"))
    end
end

"""
    load_and_prepare_data(cmdir, input_file, xyc, zc)

Helper function to load HDF5 data and apply fiducial cuts.

# Returns
- `DataFrame`: Processed hits DataFrame after fiducial cuts
"""
function load_and_prepare_data(cmdir, input_file, xyc, zc)
    # Load data file
    xfile = joinpath(cmdir, input_file)
    if !isfile(xfile)
        throw(ArgumentError("Input file not found: $xfile"))
    end

    println(" Loading data from: $xfile")
    dfs = get_dataset_dfs(xfile)
    hitsdf = dfs["hits"]

    ###
    ### Aplying fiducial cuts to large files is very slow
    ### better to apply them at analysis level
    ####

    #println(" Applying fiducial volume cuts (xyc=$xyc, zc=$zc)")
    #hitsdf_fiducial = filter_fiducial_events(hitsdf, xyc, zc)

    # Get actual number of unique events after fiducial cuts
    #n_events_available = length(unique(hitsdf_fiducial.event_id))
    #println(" Found $n_events_available events after fiducial cuts")

    #return hitsdf_fiducial
    return hitsdf
end

"""
    event_loop(cmdir; input_file="0nubb.next.h5", events_to_run=100, voxel_size_mm=5,
               max_distance_mm=10, energy_threshold_kev=10, xyc=1800.0, zc=100.0)

Process HDF5 data file through the full analysis pipeline with fiducial volume cuts.

# Arguments
- `cmdir`: Directory containing the input file
- `input_file="0nubb.next.h5"`: Name of the HDF5 input file
- `events_to_run=100`: Number of events to process
- `voxel_size_mm=5`: Voxel size in mm for spatial discretization
- `max_distance_mm=10`: Maximum distance in mm for connecting voxels into tracks
- `energy_threshold_kev=10`: Minimum energy threshold in keV for including voxels
- `xyc=1800.0`: Fiducial cut value for x and y coordinates (mm)
- `zc=100.0`: Fiducial cut value for z coordinate (mm)

# Returns
- `AnalysisResults`: Analysis results with track statistics and energy distributions
"""
function event_loop(cmdir; input_file="0nubb.next.h5",
                    events_to_run=100,
                    voxel_size_mm=5,
                    max_distance_mm=10,
                    energy_threshold_kev=10,
                    xyc::Float64=1950.0,
                    zc::Float64=10.0)

    # Validate input parameters
    validate_event_loop_parameters(events_to_run, voxel_size_mm, max_distance_mm, energy_threshold_kev)

    # Load and prepare data
    hitsdf_fiducial = load_and_prepare_data(cmdir, input_file, xyc, zc)

    println(" Starting analysis of up to $events_to_run events")
    results = analysis_loop(hitsdf_fiducial;
                           events_to_run=events_to_run,  # Pass integer to process first N events
                           voxel_size_mm,
                           max_distance_mm,
                           energy_threshold_kev)
    return results
end

"""
    event_loop_single_track(cmdir; input_file="0nubb.next.h5", events_to_run=100, voxel_size_mm=5,
                            max_distance_mm=10, energy_threshold_kev=10, xyc=1950.0, zc=10.0)

Process HDF5 data file to extract only single-track events.

# Arguments
- Same as `event_loop`

# Returns
- `Vector{Tracks}`: Vector containing only single-track events
"""
function event_loop_single_track(cmdir; input_file="0nubb.next.h5",
                    events_to_run=100,
                    voxel_size_mm=5,
                    max_distance_mm=10,
                    energy_threshold_kev=10,
                    xyc::Float64=1950.0,
                    zc::Float64=10.0)

    # Validate input parameters
    validate_event_loop_parameters(events_to_run, voxel_size_mm, max_distance_mm, energy_threshold_kev)

    # Load and prepare data
    hitsdf_fiducial = load_and_prepare_data(cmdir, input_file, xyc, zc)

    println(" Starting analysis of up to $events_to_run events")
    results = analysis_loop_single_track(hitsdf_fiducial;
                           events_to_run=events_to_run,  # Pass integer to process first N events
                           voxel_size_mm,
                           max_distance_mm,
                           energy_threshold_kev)
    return results
end


"""
    analysis_loop_single_track2(hitsdf::DataFrame;
                               events_to_run=nothing,
                               initial_event=1,
                               show_progress=false,
                               voxel_size_mm::Float64=2.0,
                               max_distance_mm::Float64=5.0,
                               energy_threshold_kev::Float64=1.0,
                               emin::Float64=-Inf,
                               emax::Float64=Inf)

Analyze events to extract single-track events with support for starting from a specific event
and showing progress.

# Arguments
- `hitsdf::DataFrame`: DataFrame containing hit data
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
function analysis_loop_single_track2(hitsdf::DataFrame;
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
        prog = Progress(length(event_ids_to_process), 1, "Processing events: ")
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
            #else
                #println("number of tracks found = $(length(tracks))")
            end

        catch e
            println("Warning: Error processing event $nevent: $e")
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

    println(" Analysis completed! Processed $n_events_processed events.")
    println(" Number of single track events $n_single_track")

    return TRACKS
end


"""
    event_loop_single_track2(cmdir; input_file="0nubb.next.h5",
                            events_to_run=100,
                            initial_event=1,
                            show_progress=false,
                            voxel_size_mm=5,
                            max_distance_mm=10,
                            energy_threshold_kev=10,
                            xyc::Float64=1950.0,
                            zc::Float64=10.0)

Process HDF5 data file to extract only single-track events, with support for starting
from a specific event and showing progress.

# Arguments
- `cmdir`: Directory containing the input file
- `input_file="0nubb.next.h5"`: Name of the HDF5 input file
- `events_to_run=100`: Number of events to process
- `initial_event=1`: First event to process (1-indexed, skips initial_event-1 events)
- `show_progress=false`: If true, display a progress bar
- `voxel_size_mm=5`: Voxel size in mm for spatial discretization
- `max_distance_mm=10`: Maximum distance in mm for connecting voxels into tracks
- `energy_threshold_kev=10`: Minimum energy threshold in keV for including voxels
- `xyc=1950.0`: Fiducial cut value for x and y coordinates (mm)
- `zc=10.0`: Fiducial cut value for z coordinate (mm)

# Returns
- `Vector{Tracks}`: Vector containing only single-track events

# Example
```julia
# Process events 1000 to 2000
tracks = event_loop_single_track2(cmdir;
                                 input_file="input.h5",
                                 events_to_run=1001,
                                 initial_event=1000,
                                 show_progress=true)
```
"""
function event_loop_single_track2(cmdir; input_file="0nubb.next.h5",
                    events_to_run=100,
                    initial_event::Int=1,
                    show_progress::Bool=false,
                    voxel_size_mm=5,
                    max_distance_mm=10,
                    energy_threshold_kev=10,
                    xyc::Float64=1950.0,
                    zc::Float64=10.0)

    # Validate input parameters
    validate_event_loop_parameters(events_to_run, voxel_size_mm, max_distance_mm, energy_threshold_kev)

    # Load and prepare data
    hitsdf_fiducial = load_and_prepare_data(cmdir, input_file, xyc, zc)

    println(" Starting analysis from event $initial_event, processing $events_to_run events")
    results = analysis_loop_single_track2(hitsdf_fiducial;
                           events_to_run=events_to_run,
                           initial_event=initial_event,
                           show_progress=show_progress,
                           voxel_size_mm,
                           max_distance_mm,
                           energy_threshold_kev)
    return results
end


"""
    histogram_results(df::DataFrame; nbins=50, xrange=(-2000.0, 2000.0), 
                     yrange=(-2000.0, 2000.0), zrange=(0.0, 4000.0), erange=(0.0, 2700.0))

Create histograms from a DataFrame containing track data.

# Arguments
- `df::DataFrame`: DataFrame with columns: event_id, energy, x, y, z
- `nbins`: Number of bins for histograms
- `xrange`, `yrange`, `zrange`, `erange`: Ranges for each histogram

# Returns
- Tuple of (h_x, h_y, h_z, h_e) histograms
"""
function histogram_results(df::DataFrame; 
                           nbins = 50,
                           xrange =  (-2000.0, 2000.0),
                           yrange =  (-2000.0, 2000.0),
                           zrange =  (0.0, 4000.0),
                           erange =  (0.0, 2700.0))
	
	# Extract unique energies per event (since energy is repeated for each voxel)
	unique_energies = unique(select(df, [:event_id, :energy]))[:, :energy]
	
	h_x = hist1d(df.x; nbins = nbins, xlim  = xrange)
	h_y = hist1d(df.y; nbins = nbins, xlim  = yrange)
	h_z = hist1d(df.z; nbins = nbins, xlim  = zrange)
	h_e = hist1d(unique_energies; nbins = nbins, xlim = erange )

	return(h_x, h_y, h_z, h_e)
end


"""
    smear_histogram(h::Histo1d, sigma::Float64)

Apply Gaussian smearing to a histogram. For each bin with center `xc` and weight `w`,
generate `w` random numbers from a Gaussian distribution with mean `xc` and standard 
deviation `sigma`.

# Arguments
- `h::Histo1d`: Input histogram to smear
- `sigma::Float64`: Standard deviation of the Gaussian smearing

# Returns
- `Vector{Float64}`: Vector containing all smeared values
"""
function smear_histogram(h::Histo1d, sigma::Float64)
    smeared_values = Float64[]
    
    for (i, weight) in enumerate(h.weights)
        if weight > 0  # Only process bins with non-zero weights
            # Get the center of the bin
            xc = h.centers[i]
            
            # Generate weight number of random Gaussian samples
            n_samples = round(Int, weight)
            if n_samples > 0
                # Generate Gaussian random numbers with mean=xc and std=sigma
                samples = randn(n_samples) .* sigma .+ xc
                append!(smeared_values, samples)
            end
        end
    end
    
    return smeared_values
end


"""
    example_smearing()

Demonstration of histogram smearing functionality.
Creates a sample histogram and applies Gaussian smearing with different sigma values.

# Returns
- Tuple of (original_histogram, smeared_values_sigma1, smeared_values_sigma5)
"""
function example_smearing()
    # Create example data: a peak at 2500 keV with some background
    energies = vcat(
        fill(2500.0, 100),    # Peak at 2500 keV
        fill(2400.0, 20),     # Some counts at 2400 keV
        fill(2600.0, 30)      # Some counts at 2600 keV
    )
    
    # Create histogram
    h = hist1d(energies; nbins=50, xlim=(2300.0, 2700.0))
    
    # Apply different amounts of smearing
    smeared_1keV = smear_histogram(h, 1.0)    # 1 keV resolution
    smeared_5keV = smear_histogram(h, 5.0)    # 5 keV resolution
    
    println("Original histogram: $(length(energies)) total counts")
    println("Smeared with σ=1 keV: $(length(smeared_1keV)) samples")
    println("Smeared with σ=5 keV: $(length(smeared_5keV)) samples")
    
    return h, smeared_1keV, smeared_5keV
end


"""
    counts_in_range(h::Histo1d, xlow::Float64, xup::Float64)

Count the total number of counts (sum of weights) in a histogram within a specified range.

This function sums the weights of all bins whose centers fall within the range [xlow, xup].
Bins are included if their center satisfies: xlow ≤ bin_center ≤ xup.

# Arguments
- `h::Histo1d`: Input histogram
- `xlow::Float64`: Lower bound of the range (inclusive)
- `xup::Float64`: Upper bound of the range (inclusive)

# Returns
- `Float64`: Total counts (sum of weights) in the specified range

# Example
```julia
# Create a histogram
h = hist1d(data; nbins=50, xlim=(0.0, 100.0))

# Count entries between 20 and 30
counts = counts_in_range(h, 20.0, 30.0)
```
"""
function counts_in_range(h::Histo1d, xlow::Float64, xup::Float64)
    if xlow > xup
        throw(ArgumentError("xlow ($xlow) must be less than or equal to xup ($xup)"))
    end
    
    total_counts = 0.0
    
    for (i, center) in enumerate(h.centers)
        if xlow <= center <= xup
            total_counts += h.weights[i]
        end
    end
    
    return total_counts
end


"""
    save_analysis_results(results::AnalysisResults, output_dir::String)

Save AnalysisResults to a directory containing:
- DataFrames as HDF5 files (single_track.h5, two_track_primary.h5, etc.)
- Metadata as JSON file (metadata.json)

# Arguments
- `results::AnalysisResults`: The results structure to save
- `output_dir::String`: Directory path where files will be created

# Output Structure
```
output_dir/
├── single_track.h5
├── two_track_primary.h5  
├── two_track_secondary.h5
├── three_track_primary.h5
├── three_track_secondary.h5
└── metadata.json
```
"""
function save_analysis_results(results::AnalysisResults, output_dir::String)
    # Create output directory if it doesn't exist
    if !isdir(output_dir)
        mkpath(output_dir)
        println("Created directory: $output_dir")
    end
    
    # Save DataFrames as HDF5 files
    track_types = [:single_track, :two_track_primary, :two_track_secondary, 
                   :three_track_primary, :three_track_secondary]
    
    for track_type in track_types
        df = getfield(results, track_type)
        h5_file = joinpath(output_dir, "$(track_type).h5")
        
        # Save DataFrame to HDF5
        h5open(h5_file, "w") do file
            # Save each column as a dataset
            if nrow(df) > 0
                for col_name in names(df)
                    file[string(col_name)] = df[!, col_name]
                end
                # Save column names and types as attributes
                attrs(file)["column_names"] = collect(names(df))
                attrs(file)["column_types"] = [string(eltype(df[!, col])) for col in names(df)]
            else
                # Handle empty DataFrames
                attrs(file)["empty"] = true
                attrs(file)["column_names"] = collect(names(df))
                attrs(file)["column_types"] = [string(eltype(df[!, col])) for col in names(df)]
            end
        end
    end
    
    # Save metadata as JSON
    metadata = Dict(
        "n_events_processed" => results.n_events_processed,
        "n_single_track" => results.n_single_track,
        "n_two_track" => results.n_two_track,
        "n_three_plus_track" => results.n_three_plus_track,
        "n_failed" => results.n_failed,
        "format_version" => "1.0",
        "created_at" => string(now())
    )
    
    json_file = joinpath(output_dir, "metadata.json")
    open(json_file, "w") do io
        JSON.print(io, metadata, 2)  # Pretty print with 2-space indent
    end
    
    println("Analysis results saved to directory: $output_dir")
    println("Files created:")
    for track_type in track_types
        println("  - $(track_type).h5")
    end
    println("  - metadata.json")
end


"""
    load_analysis_results(input_dir::String) -> AnalysisResults

Load AnalysisResults from a directory containing HDF5 and JSON files.

# Arguments
- `input_dir::String`: Directory path containing the analysis results files

# Expected Structure
```
input_dir/
├── single_track.h5
├── two_track_primary.h5  
├── two_track_secondary.h5
├── three_track_primary.h5
├── three_track_secondary.h5
└── metadata.json
```

# Returns
- `AnalysisResults`: The reconstructed results structure

# Example
```julia
results = load_analysis_results("path/to/results_directory")
println("Loaded results with \$(results.n_events_processed) events")
```
"""
function load_analysis_results(input_dir::String)
    if !isdir(input_dir)
        error("Directory not found: $input_dir")
    end
    
    # Load metadata
    metadata_file = joinpath(input_dir, "metadata.json")
    if !isfile(metadata_file)
        error("Metadata file not found: $metadata_file")
    end
    
    metadata = JSON.parsefile(metadata_file)
    
    # Load DataFrames from HDF5 files
    track_types = [:single_track, :two_track_primary, :two_track_secondary, 
                   :three_track_primary, :three_track_secondary]
    
    loaded_dfs = Dict{Symbol, DataFrame}()
    
    for track_type in track_types
        h5_file = joinpath(input_dir, "$(track_type).h5")
        if !isfile(h5_file)
            error("HDF5 file not found: $h5_file")
        end
        
        # Load DataFrame from HDF5
        h5open(h5_file, "r") do file
            if haskey(attrs(file), "empty") && attrs(file)["empty"]
                # Handle empty DataFrame
                col_names = attrs(file)["column_names"]
                col_types = attrs(file)["column_types"]
                
                # Create empty DataFrame with correct column types
                df = DataFrame()
                for (name, type_str) in zip(col_names, col_types)
                    col_type = eval(Symbol(type_str))  # Convert string back to type
                    df[!, Symbol(name)] = col_type[]
                end
                loaded_dfs[track_type] = df
            else
                # Load data from HDF5
                col_names = attrs(file)["column_names"]
                df_data = Dict()
                
                for col_name in col_names
                    df_data[Symbol(col_name)] = read(file, col_name)
                end
                
                loaded_dfs[track_type] = DataFrame(df_data)
            end
        end
    end
    
    # Construct AnalysisResults
    return AnalysisResults(
        loaded_dfs[:single_track],
        loaded_dfs[:two_track_primary],
        loaded_dfs[:two_track_secondary],
        loaded_dfs[:three_track_primary],
        loaded_dfs[:three_track_secondary],
        metadata["n_events_processed"],
        metadata["n_single_track"],
        metadata["n_two_track"],
        metadata["n_three_plus_track"],
        metadata["n_failed"]
    )
end


