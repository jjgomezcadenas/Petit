mutable struct Eff
	eff1tr::Float64
	efffid::Float64
	effroi::Float64
	effblb::Float64
end

mutable struct TRK
	tracks::Vector{Tracks}
	metas::Vector{Dict{String, Any}}
	ntot::Int
	n1trk::Int
end

mutable struct NTRKS
	ntot::Int64  # events with 1 track, Copper + Inner
	nfid::Int64  # events, pass cut on confidence, track length, eblob2
	nroi::Int64  # events in ROI
	nblobs::Int64  # events with eblob2 > ecut
end

struct Blobs
	confidence::Vector{Float64}
	trackLength::Vector{Float64}
	energyKeV::Vector{Float64}
	eB1::Vector{Float64}
	eB2::Vector{Float64}
end

struct BlobsVariable
	confidence::Vector{Float64}
	trackLength::Vector{Float64}
	energyKeV::Vector{Float64}
	eB1::Vector{Float64}
	eB2::Vector{Float64}
	rB1::Vector{Float64}  # Optimal radius for blob1
	rB2::Vector{Float64}  # Optimal radius for blob2
end


function get_bkgnd_tracks(cmdir;isotope::String="bi214", 
						  bkgnd::String ="copperbkg", 
						  tag::String ="*st3mm*")

	path = joinpath(cmdir, bkgnd, isotope)
	println(path)
	files = Glob.glob(tag, path)
	tracks, metadatas = chain_track_files(files)
	
	ntot = 0
	for i in 1:length(metadatas)
		ntot+= metadatas[i]["nevents_from_config"]
	end
	
	return TRK(tracks, metadatas, ntot, length(tracks))
end


function get_bb0nu_tracks(cmdir; tag::String ="*st3mm*")
	path = joinpath(cmdir, "bb0nu")
	files = Glob.glob(tag, path)
	tracks, metadatas = chain_track_files(files)
	ntot = 0
	for i in 1:length(metadatas)
		ntot += metadatas[i]["events_processed"] # NB only one file for bb
	end

	return TRK(tracks, metadatas, ntot, length(tracks))
end


function get_xe137_tracks(cmdir; xedir="xe137r2", tag="st3mm")
	path = joinpath(cmdir, xedir)
	# Construct glob pattern from tag
	pattern = "*$(tag)*"
	files = Glob.glob(pattern, path)
	tracks, metadatas = chain_track_files(files)
	ntot = 0
	for i in 1:length(metadatas)
		ntot += metadatas[i]["events_processed"] # NB only one file for bb
	end

	return TRK(tracks, metadatas, ntot, length(tracks))
end


"""
    get_itaca_tracks(cmdir; tag::String="*")

Read ITACA MC, ion and electron tracks from HDF5 files.

# Arguments
- `cmdir::String`: Directory containing the track files
- `tag::String`: Glob pattern to match files (default: "*")

# Returns
- `(mc_trk::TRK, ion_trk::TRK, ele_trk::TRK)`: Tuple of TRK structs for MC, ion and electron tracks
"""
function get_itaca_tracks(cmdir::String; tag::String="*")
    # Find mc, ion and electron files
    mc_pattern = "*$(tag)*_mc.h5"
    ion_pattern = "*$(tag)*_ion.h5"
    ele_pattern = "*$(tag)*_ele.h5"

    mc_files = Glob.glob(mc_pattern, cmdir)
    ion_files = Glob.glob(ion_pattern, cmdir)
    ele_files = Glob.glob(ele_pattern, cmdir)

    if isempty(mc_files)
        @warn "No MC track files found matching pattern: $mc_pattern in $cmdir"
    end
    if isempty(ion_files)
        @warn "No ion track files found matching pattern: $ion_pattern in $cmdir"
    end
    if isempty(ele_files)
        @warn "No electron track files found matching pattern: $ele_pattern in $cmdir"
    end

    # Chain MC tracks
    mc_tracks, mc_metas = if !isempty(mc_files)
        chain_track_files(mc_files)
    else
        (Tracks[], Dict{String, Any}[])
    end

    # Chain ion tracks
    ion_tracks, ion_metas = if !isempty(ion_files)
        chain_track_files(ion_files)
    else
        (Tracks[], Dict{String, Any}[])
    end

    # Chain electron tracks
    ele_tracks, ele_metas = if !isempty(ele_files)
        chain_track_files(ele_files)
    else
        (Tracks[], Dict{String, Any}[])
    end

    # Count total events processed
    mc_ntot = sum(m["events_processed"] for m in mc_metas; init=0)
    ion_ntot = sum(m["events_processed"] for m in ion_metas; init=0)
    ele_ntot = sum(m["events_processed"] for m in ele_metas; init=0)

    mc_trk = TRK(mc_tracks, mc_metas, mc_ntot, length(mc_tracks))
    ion_trk = TRK(ion_tracks, ion_metas, ion_ntot, length(ion_tracks))
    ele_trk = TRK(ele_tracks, ele_metas, ele_ntot, length(ele_tracks))

    return (mc_trk, ion_trk, ele_trk)
end


"""
    read_itaca_metadata(filepath::String)

Read ITACA run metadata from CSV file into a DataFrame.

# Arguments
- `filepath::String`: Path to the metadata CSV file

# Returns
- `DataFrame`: DataFrame with parameter and value columns
"""
function read_itaca_metadata(filepath::String)
    CSV.read(filepath, DataFrame)
end


"""
    get_itaca_metadata_value(df::DataFrame, param::String)

Get a specific parameter value from metadata DataFrame.

# Arguments
- `df::DataFrame`: Metadata DataFrame from read_itaca_metadata
- `param::String`: Parameter name to retrieve

# Returns
- Value as appropriate type (String, Int, or Float64)
"""
function get_itaca_metadata_value(df::DataFrame, param::String)
    row = filter(r -> r.parameter == param, df)
    if nrow(row) == 0
        error("Parameter '$param' not found in metadata")
    end
    val = row.value[1]

    # Try to convert to appropriate type
    if param in ["first_event", "last_event", "nthreads", "nbins_df"]
        return parse(Int, val)
    elseif param in ["directory", "input_file", "output_base"]
        return val
    else
        return parse(Float64, val)
    end
end


"""
    match_itaca_events(mc_trk::TRK, ion_trk::TRK, ele_trk::TRK)

Match event IDs across MC, ion, and electron tracks.

For each unique event ID found in any of the three track collections,
returns the indices (mc_idx, ion_idx, ele_idx) where that event appears.
If an event is not present in a collection, -1 is returned for that index.

# Arguments
- `mc_trk::TRK`: MC tracks
- `ion_trk::TRK`: Ion tracks
- `ele_trk::TRK`: Electron tracks

# Returns
- `Dict{Int, Tuple{Int, Int, Int}}`: Dictionary mapping event_id to (mc_idx, ion_idx, ele_idx)
"""
function match_itaca_events(mc_trk::TRK, ion_trk::TRK, ele_trk::TRK)
    # Build index maps: event_id -> track index
    mc_map = Dict{Int, Int}()
    for (idx, track) in enumerate(mc_trk.tracks)
        if nrow(track.voxels) > 0
            event_id = track.voxels.event_id[1]
            mc_map[event_id] = idx
        end
    end

    ion_map = Dict{Int, Int}()
    for (idx, track) in enumerate(ion_trk.tracks)
        if nrow(track.voxels) > 0
            event_id = track.voxels.event_id[1]
            ion_map[event_id] = idx
        end
    end

    ele_map = Dict{Int, Int}()
    for (idx, track) in enumerate(ele_trk.tracks)
        if nrow(track.voxels) > 0
            event_id = track.voxels.event_id[1]
            ele_map[event_id] = idx
        end
    end

    # Collect all unique event IDs
    all_event_ids = union(keys(mc_map), keys(ion_map), keys(ele_map))

    # Build result dictionary
    result = Dict{Int, Tuple{Int, Int, Int}}()
    for event_id in all_event_ids
        mc_idx = get(mc_map, event_id, -1)
        ion_idx = get(ion_map, event_id, -1)
        ele_idx = get(ele_map, event_id, -1)
        result[event_id] = (mc_idx, ion_idx, ele_idx)
    end

    return result
end


"""
    match_itaca_events_df(mc_trk::TRK, ion_trk::TRK, ele_trk::TRK)

Match event IDs across MC, ion, and electron tracks, returning a DataFrame.

# Returns
- `DataFrame`: With columns event_id, mc_idx, ion_idx, ele_idx, sorted by event_id
"""
function match_itaca_events_df(mc_trk::TRK, ion_trk::TRK, ele_trk::TRK)
    matches = match_itaca_events(mc_trk, ion_trk, ele_trk)

    df = DataFrame(
        event_id = collect(keys(matches)),
        mc_idx = [m[1] for m in values(matches)],
        ion_idx = [m[2] for m in values(matches)],
        ele_idx = [m[3] for m in values(matches)]
    )

    sort!(df, :event_id)
    return df
end


function blob_analysis(strks::Vector{Tracks}, r::Float64; nmax::Int=-1, nprint::Int=0)
	eB1=Float64[]
	eB2=Float64[]
	CON = Float64[]
	TL = Float64[]

	# Determine how many tracks to process
	ntracks_total = length(strks)
	ntracks_to_process = nmax > 0 ? min(nmax, ntracks_total) : ntracks_total

	println(stderr, "Processing $ntracks_to_process out of $ntracks_total tracks")
	flush(stderr)

	# Process only the first ntracks_to_process tracks
	for (i, strk) in enumerate(strks[1:ntracks_to_process])
		xresult = walk_track_from_extremes(strk)
	  	xtrack_length = xresult.total_length
		confidence = xresult.confidence

		blobs = energy_in_spheres_around_extremes(strk, xresult, r)
		eb1 = blobs.blob1_energy * 1e+3
		eb2 = blobs.blob2_energy * 1e+3
		push!(eB1, eb1)
		push!(eB2, eb2)
		push!(CON, confidence)
		push!(TL, xtrack_length)

		# Print progress every nprint tracks
		if nprint > 0 && i % nprint == 0
			println(stderr, "Processed $i/$ntracks_to_process tracks")
            flush(stderr)
		end
	end

	# Calculate energies only for processed tracks
	E = track_energies_keV(strks[1:ntracks_to_process])

	return Blobs(CON, TL, E, eB1, eB2)
end


function blob_analysis_variable(strks::Vector{Tracks};
                                seed_radius::Float64=3.0,
                                step::Float64=1.0,
                                max_radius::Float64=10.0,
                                threshold::Float64=0.05,
                                nmax::Int=-1,
                                nprint::Int=0)
	eB1 = Float64[]
	eB2 = Float64[]
	rB1 = Float64[]
	rB2 = Float64[]
	CON = Float64[]
	TL = Float64[]

	# Determine how many tracks to process
	ntracks_total = length(strks)
	ntracks_to_process = nmax > 0 ? min(nmax, ntracks_total) : ntracks_total

	println(stderr, "Processing $ntracks_to_process out of $ntracks_total tracks (variable radius mode)")
	flush(stderr)

	# Process only the first ntracks_to_process tracks
	for (i, strk) in enumerate(strks[1:ntracks_to_process])
		xresult = walk_track_from_extremes(strk)
		xtrack_length = xresult.total_length
		confidence = xresult.confidence

		blobs = energy_in_variable_spheres_around_extremes(
			strk, xresult;
			seed_radius=seed_radius,
			step=step,
			max_radius=max_radius,
			threshold=threshold
		)

		eb1 = blobs.blob1_energy * 1e+3
		eb2 = blobs.blob2_energy * 1e+3
		rb1 = blobs.blob1_radius
		rb2 = blobs.blob2_radius

		push!(eB1, eb1)
		push!(eB2, eb2)
		push!(rB1, rb1)
		push!(rB2, rb2)
		push!(CON, confidence)
		push!(TL, xtrack_length)

		# Print progress every nprint tracks
		if nprint > 0 && i % nprint == 0
			println(stderr, "Processed $i/$ntracks_to_process tracks")
			flush(stderr)
		end
	end

	# Calculate energies only for processed tracks
	E = track_energies_keV(strks[1:ntracks_to_process])

	return BlobsVariable(CON, TL, E, eB1, eB2, rB1, rB2)
end


function write_blobs_csv(blobs::Blobs, filename::String)
      df = DataFrame(
          confidence = blobs.confidence,
          trackLength = blobs.trackLength,
          energyKeV = blobs.energyKeV,
          eB1 = blobs.eB1,
          eB2 = blobs.eB2
      )
      CSV.write(filename, df)
      #println("Blobs written to $filename")
  end


  function fom_blobs(eblob1::Vector{Float64}, eblob2::Vector{Float64}; min_cut=100.0, max_cut=1000.0,step=50.0)

	function simple_binomial_error(k::Int, n::Int)
	    if n == 0
	        return (0.0, 0.0)
	    end
	
	    p = k / n
	
	    # Standard error for binomial proportion
	    # σ = sqrt(p(1-p)/n)
	    error = sqrt(p * (1 - p) / n)
	
	    return (p, error)
	end
    # Total number of events below max_cut
    n_total = length(eblob2)
	

    # Initialize results arrays
    cuts = Float64[]
    efficiencies = Float64[]
    errors = Float64[]

    # Iterate over cut values
    cut_value = min_cut
    while cut_value <= max_cut
        # Count events passing the cut
        n_pass = sum(eblob2 .>= cut_value)

        # Calculate efficiency and error
        
        eff, err = simple_binomial_error(n_pass, n_total)
        

        push!(cuts, cut_value)
        push!(efficiencies, eff)
        push!(errors, err)

        cut_value += step
    end

    # Create output dataframe
    results_df = DataFrame(
        eblob2 = cuts,
        eff = round.(efficiencies, digits=4),
        err = round.(errors, digits=4)
    )

    return results_df
end