struct TracksSummary
    energies::Vector{Float64}
    xs::Vector{Float64}
    ys::Vector{Float64}
    zs::Vector{Float64}
    ids::Vector{Int}
end

# Constructor for empty TracksSummary
TracksSummary() = TracksSummary(Float64[], Float64[], Float64[], Float64[], Int[])

struct AnalysisResults
    single_track::TracksSummary
    two_track_primary::TracksSummary
    two_track_secondary::TracksSummary
    three_track_primary::TracksSummary
    three_track_secondary::TracksSummary
    n_events_processed::Int
    n_single_track::Int
    n_two_track::Int
    n_three_plus_track::Int
    n_failed::Int
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
        if isempty(event_ids_to_process)
            println("Warning: None of the requested event IDs exist in the DataFrame")
            println("Available event IDs range from $(minimum(unique_event_ids)) to $(maximum(unique_event_ids))")
        end
    else
        throw(ArgumentError("events_to_run must be nothing, an Integer, a Range, or a Vector"))
    end
    
    # Initialize energy arrays with proper types
    single_tracks = TracksSummary()
    two_track_primary = TracksSummary()
    two_track_secondary = TracksSummary()
    three_track_primary = TracksSummary()
    three_track_secondary = TracksSummary()
    
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
        return AnalysisResults(
            single_tracks, two_track_primary, two_track_secondary,
            three_track_primary, three_track_secondary,
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
                push!(single_tracks.energies, energy_kev)
                push!(single_tracks.ids, nevent)
                append!(single_tracks.xs, tracks[1].voxels.x)
                append!(single_tracks.ys, tracks[1].voxels.y)
                append!(single_tracks.zs, tracks[1].voxels.z)

                n_single_track += 1
                
            elseif length(tracks) == 2 
                # Two track event
                primary_energy = 1e+3 * sum(tracks[1].voxels.energy)
                secondary_energy = 1e+3 * sum(tracks[2].voxels.energy)
                push!(two_track_primary.energies, primary_energy)
                push!(two_track_primary.ids, nevent)
                append!(two_track_primary.xs, tracks[1].voxels.x)
                append!(two_track_primary.ys, tracks[1].voxels.y)
                append!(two_track_primary.zs, tracks[1].voxels.z)
                push!(two_track_secondary.energies, secondary_energy)
                push!(two_track_secondary.ids, nevent)
                append!(two_track_secondary.xs, tracks[2].voxels.x)
                append!(two_track_secondary.ys, tracks[2].voxels.y)
                append!(two_track_secondary.zs, tracks[2].voxels.z)
                
                n_two_track += 1
                
            elseif length(tracks) >= 3 
                # Three or more track event
                primary_energy = 1e+3 * sum(tracks[1].voxels.energy)
                push!(three_track_primary.energies, primary_energy)
                push!(three_track_primary.ids, nevent)
                append!(three_track_primary.xs, tracks[1].voxels.x)
                append!(three_track_primary.ys, tracks[1].voxels.y)
                append!(three_track_primary.zs, tracks[1].voxels.z)
                for n in 2:length(tracks)
                    secondary_energy = 1e+3 * sum(tracks[n].voxels.energy)
                    push!(three_track_secondary.energies, secondary_energy)
                    push!(three_track_secondary.ids, nevent)
                    append!(three_track_secondary.xs, tracks[n].voxels.x)
                    append!(three_track_secondary.ys, tracks[n].voxels.y)
                    append!(three_track_secondary.zs, tracks[n].voxels.z)
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
    println("✅ Analysis completed! Processed $n_events_processed events.")
    
    return AnalysisResults(
        single_tracks,
        two_track_primary,
        two_track_secondary,
        three_track_primary,
        three_track_secondary,
        n_events_processed,
        n_single_track,
        n_two_track,
        n_three_plus_track,
        n_failed
    )

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
                    xyc::Float64=1800.0,
                    zc::Float64=100.0)
    
    # Validate input parameters
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
    
    # Load data file
    xfile = joinpath(cmdir, input_file)
    if !isfile(xfile)
        throw(ArgumentError("Input file not found: $xfile"))
    end
    
    println("📁 Loading data from: $xfile")
    dfs = get_dataset_dfs(xfile)
    hitsdf = dfs["hits"]
    
    println("🔍 Applying fiducial volume cuts (xyc=$xyc, zc=$zc)")
    hitsdf_fiducial = filter_fiducial_events(hitsdf, xyc, zc)
    
    # Get actual number of unique events after fiducial cuts
    n_events_available = length(unique(hitsdf_fiducial.event_id))
    println("📊 Found $n_events_available events after fiducial cuts")
    
    println("🔬 Starting analysis of up to $events_to_run events")
    results = analysis_loop(hitsdf_fiducial; 
                           events_to_run=events_to_run,  # Pass integer to process first N events
                           voxel_size_mm,
                           max_distance_mm, 
                           energy_threshold_kev)
    return results
end




function histogram_results(ts::TracksSummary; 
                           nbins = 50,
                           xrange =  (-2000.0, 2000.0),
                           yrange =  (-2000.0, 2000.0),
                           zrange =  (0.0, 4000.0),
                           erange =  (0.0, 2700.0))
	
	h_x = hist1d(ts.xs; nbins = nbins, xlim  = xrange)
	h_y = hist1d(ts.ys; nbins = nbins, xlim  = yrange)
	h_z = hist1d(ts.zs; nbins = nbins, xlim  = zrange)
	h_e = hist1d(ts.energies; nbins = nbins, xlim = erange )

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
    save_analysis_results(results::AnalysisResults, filename::String)

Save an AnalysisResults structure to disk using Julia's built-in serialization.

# Arguments
- `results::AnalysisResults`: The results structure to save
- `filename::String`: Path to the output file (should end with .jls)

# Example
```julia
results = event_loop(cmdir; events_to_run=1000)
save_analysis_results(results, "analysis_results.jls")
```
"""
function save_analysis_results(results::AnalysisResults, filename::String)
    open(filename, "w") do io
        Serialization.serialize(io, results)
    end
    println("Analysis results saved to: $filename")
end


"""
    load_analysis_results(filename::String) -> AnalysisResults

Load an AnalysisResults structure from disk using Julia's built-in serialization.

# Arguments
- `filename::String`: Path to the input file (should end with .jls)

# Returns
- `AnalysisResults`: The loaded results structure

# Example
```julia
results = load_analysis_results("analysis_results.jls")
println("Loaded results with \$(results.n_events_processed) events")
```
"""
function load_analysis_results(filename::String)
    if !isfile(filename)
        error("File not found: $filename")
    end
    
    return open(filename, "r") do io
        results = Serialization.deserialize(io)
        if !isa(results, AnalysisResults)
            error("File does not contain valid AnalysisResults data")
        end
        return results
    end
end


