"""
Auxiliary functions for itaca
"""

using DataFrames
using Statistics
#using Plots
#using Images
#using ImageFiltering
#using Interpolations
using Random
using Distributions

"""
sigma_t (mm) = dtmm /sqrt(Pbar) x sqrt(Lcm)
"""
function sigma_t_mm(Lcm, Pbar; dtmm=3.5) 
    dtmm /sqrt(Pbar) * sqrt(Lcm)
end

"""
sigma_l (mm) = dlmm /sqrt(Pbar) x sqrt(Lcm)
"""
function sigma_l_mm(Lcm, Pbar; dlmm=0.9) 
    dlmm /sqrt(Pbar) * sqrt(Lcm)
end


"""
sigma_t (mm) = sqrt(2Kb/e x T x L/E) * 10
where:
- 2Kb/e in V/K
- T in K
- L in cm
- E in V/cm
"""
function sigma_t_ion_mm(Tk, Lcm, Evcm; kB_ovr_e = 8.6173E-5)  
    # T -> K, L -> cm, E -> v/cm
    return sqrt(2*kB_ovr_e*Tk*(Lcm/Evcm)) * 10 #-> mm
end



"""
    transform_hits_df(df; energy_to_electrons=1e5/2.5)

Drop time/label/particle_id/hit_id columns and add electrons column from energy.
"""
function transform_hits_df(df::DataFrame; energy_to_electrons::Float64=1e5/2.5)
    df2 = select(df, Not([:time, :label, :particle_id, :hit_id]))
    df2.electrons = round.(Int, df2.energy .* energy_to_electrons)
    return df2
end


"""
    diffuse_xy_image_mc(df::DataFrame; sigma_mm=0.6, nbins=100, nsigma=3.0)

Simulate Gaussian diffusion using Monte Carlo method.
For each electron in the 'electrons' column, generate random positions according to
a 2D Gaussian distribution centered at (x, y) with standard deviation sigma_mm.
Then histogram the diffused positions.

# Arguments
- `df::DataFrame`: DataFrame containing x, y, electrons columns
- `sigma_mm::Float64`: Standard deviation of Gaussian diffusion in mm (default: 0.6)
- `nbins::Int`: Number of bins for the 2D histogram (default: 100)
- `nsigma::Float64`: Number of sigmas to pad the histogram range (default: 3.0)

# Returns
- `DataFrame`: DataFrame with columns x, y, intensity representing the diffused image
"""
function diffuse_xy_image_mc(df::DataFrame; sigma_mm::Float64=0.6, nbins::Int=100, nsigma::Float64=3.0,
                             ef::Float64=1e+5/2.5)

    # Extract hit positions and electron counts
    x_hits = Float64.(df.x)
    y_hits = Float64.(df.y)

    # Check if electrons column exists, otherwise use 1 electron per hit
    if hasproperty(df, :electrons)
        electrons = Int.(df.electrons)
    else
        @warn "DataFrame does not have 'electrons' column, using 1 electron per hit"
        electrons = ones(Int, length(x_hits))
    end

    # Calculate histogram range with padding
    x_min = minimum(x_hits) - nsigma * sigma_mm
    x_max = maximum(x_hits) + nsigma * sigma_mm
    y_min = minimum(y_hits) - nsigma * sigma_mm
    y_max = maximum(y_hits) + nsigma * sigma_mm

    # Create grid edges
    x_edges = range(x_min, x_max, length=nbins+1)
    y_edges = range(y_min, y_max, length=nbins+1)

    # Bin width for checking
    dx = x_edges[2] - x_edges[1]
    dy = y_edges[2] - y_edges[1]

    # Initialize histogram
    hist_counts = zeros(Int, nbins, nbins)

    # For each hit, generate diffused electron positions
    for i in 1:length(x_hits)
        x_center = x_hits[i]
        y_center = y_hits[i]
        n_electrons = electrons[i]

        # Generate random positions for each electron using 2D Gaussian
        # In 2D, x and y are independent Gaussians with same sigma
        for _ in 1:n_electrons
            # Sample from Normal(mu, sigma)
            x_diffused = x_center + randn() * sigma_mm
            y_diffused = y_center + randn() * sigma_mm

            # Find which bin - using integer division approach
            # This is more robust than searchsortedfirst
            ix = Int(floor((x_diffused - x_min) / dx)) + 1
            iy = Int(floor((y_diffused - y_min) / dy)) + 1

            # Check if within histogram range
            if 1 <= ix <= nbins && 1 <= iy <= nbins
                hist_counts[ix, iy] += 1
            end
        end
    end

    # Get bin centers
    x_centers = (x_edges[1:end-1] .+ x_edges[2:end]) ./ 2
    y_centers = (y_edges[1:end-1] .+ y_edges[2:end]) ./ 2

    # Create output DataFrame
    result_x = Float64[]
    result_y = Float64[]
    result_intensity = Float64[]

    for ix in 1:nbins
        for iy in 1:nbins
            push!(result_x, x_centers[ix])
            push!(result_y, y_centers[iy])
            push!(result_intensity, hist_counts[ix, iy])
        end
    end

    # Get event_id (scalar, not vector)
    event_id = nrow(df) > 0 ? df.event_id[1] : 0

    return DataFrame(
        event_id = event_id,
        x = result_x,
        y = result_y,
        energy = result_intensity / ef,  # in MeV
        electrons = Int.(result_intensity)
    )
end


"""
    diffuse_xyz_image_mc(df::DataFrame; sigma_t_mm=0.6, sigma_l_mm=0.6, nbins=100, nsigma=3.0)

Simulate 3D Gaussian diffusion using Monte Carlo method with separate transverse and longitudinal diffusion.
For each electron in the 'electrons' column, generate random positions according to:
- Transverse (x, y): 2D Gaussian with standard deviation sigma_t_mm
- Longitudinal (z): 1D Gaussian with standard deviation sigma_l_mm

Then create 3D histogram of the diffused positions.

# Arguments
- `df::DataFrame`: DataFrame containing x, y, z, electrons columns
- `sigma_t_mm::Float64`: Standard deviation of transverse diffusion (x, y) in mm (default: 0.6)
- `sigma_l_mm::Float64`: Standard deviation of longitudinal diffusion (z) in mm (default: 0.6)
- `nbins::Int`: Number of bins per dimension for the 3D histogram (default: 100)
- `nsigma::Float64`: Number of sigmas to pad the histogram range (default: 3.0)

# Returns
- `DataFrame`: DataFrame with columns x, y, z, electrons representing the diffused 3D image
"""
function diffuse_xyz_image_mc(df::DataFrame; sigma_t_mm=0.6, sigma_l_mm=0.6, nbins=100, nsigma=3.0,
                              ef::Float64=1e+5/2.5)

    # Extract hit positions and electron counts
    x_hits = Float64.(df.x)
    y_hits = Float64.(df.y)
    z_hits = Float64.(df.z)

    # Check if electrons column exists, otherwise use 1 electron per hit
    if hasproperty(df, :electrons)
        electrons = Int.(df.electrons)
    elseif hasproperty(df, :energy)
        electrons = round.(Int, df.energy .* ef)
    else
        @warn "DataFrame does not have 'electrons' column, using 1 electron per hit"
        electrons = ones(Int, length(x_hits))
    end

    # Calculate histogram range with padding
    x_min = minimum(x_hits) - nsigma * sigma_t_mm
    x_max = maximum(x_hits) + nsigma * sigma_t_mm
    y_min = minimum(y_hits) - nsigma * sigma_t_mm
    y_max = maximum(y_hits) + nsigma * sigma_t_mm
    z_min = minimum(z_hits) - nsigma * sigma_l_mm
    z_max = maximum(z_hits) + nsigma * sigma_l_mm

    # Create grid edges
    x_edges = range(x_min, x_max, length=nbins+1)
    y_edges = range(y_min, y_max, length=nbins+1)
    z_edges = range(z_min, z_max, length=nbins+1)

    # Bin widths for checking
    dx = x_edges[2] - x_edges[1]
    dy = y_edges[2] - y_edges[1]
    dz = z_edges[2] - z_edges[1]

    # Initialize 3D histogram
    hist_counts = zeros(Int, nbins, nbins, nbins)

    # For each hit, generate diffused electron positions
    for i in 1:length(x_hits)
        x_center = x_hits[i]
        y_center = y_hits[i]
        z_center = z_hits[i]
        n_electrons = electrons[i]

        # Generate random positions for each electron
        # Transverse: x and y are independent Gaussians with sigma_t
        # Longitudinal: z is Gaussian with sigma_l
        for _ in 1:n_electrons
            # Sample from Normal(mu, sigma)
            x_diffused = x_center + randn() * sigma_t_mm
            y_diffused = y_center + randn() * sigma_t_mm
            z_diffused = z_center + randn() * sigma_l_mm

            # Find which bin - using integer division approach
            ix = Int(floor((x_diffused - x_min) / dx)) + 1
            iy = Int(floor((y_diffused - y_min) / dy)) + 1
            iz = Int(floor((z_diffused - z_min) / dz)) + 1

            # Check if within histogram range
            if 1 <= ix <= nbins && 1 <= iy <= nbins && 1 <= iz <= nbins
                hist_counts[ix, iy, iz] += 1
            end
        end
    end

    # Get bin centers
    x_centers = (x_edges[1:end-1] .+ x_edges[2:end]) ./ 2
    y_centers = (y_edges[1:end-1] .+ y_edges[2:end]) ./ 2
    z_centers = (z_edges[1:end-1] .+ z_edges[2:end]) ./ 2

    # Create output DataFrame with all non-zero voxels
    result_x = Float64[]
    result_y = Float64[]
    result_z = Float64[]
    result_electrons = Int[]

    for ix in 1:nbins
        for iy in 1:nbins
            for iz in 1:nbins
                count = hist_counts[ix, iy, iz]
                if count > 0  # Only store non-zero voxels
                    push!(result_x, x_centers[ix])
                    push!(result_y, y_centers[iy])
                    push!(result_z, z_centers[iz])
                    push!(result_electrons, count)
                end
            end
        end
    end

    # Get event_id (scalar, not vector)
    event_id = nrow(df) > 0 ? df.event_id[1] : 0

    return DataFrame(
        event_id = event_id,
        x = result_x,
        y = result_y,
        z = result_z,
        energy = result_electrons / ef,  # in MeV
        electrons = result_electrons
    )
end


#=
# COMMENTED OUT: This function requires ImageFiltering which has heavy dependencies
# Use diffuse_xyz_image_mc() instead (Monte Carlo version)

"""
    diffuse_xyz_image_kernel(df::DataFrame; sigma_t_mm=0.6, sigma_l_mm=0.6, nbins=100, nsigma=3.0)

Simulate 3D Gaussian diffusion using ImageFiltering.jl separable convolution.
Bins hits into a 3D histogram and applies Gaussian blur using `imfilter`.

# Arguments
- `df::DataFrame`: DataFrame containing x, y, z, electrons columns
- `sigma_t_mm::Float64`: Standard deviation of transverse diffusion (x, y) in mm (default: 0.6)
- `sigma_l_mm::Float64`: Standard deviation of longitudinal diffusion (z) in mm (default: 0.6)
- `nbins::Int`: Number of bins per dimension for the 3D histogram (default: 100)
- `nsigma::Float64`: Number of sigmas for histogram padding (default: 3.0)

# Returns
- `DataFrame`: DataFrame with columns x, y, z, electrons representing the diffused 3D image

# Performance
Uses optimized separable Gaussian convolution from ImageFiltering.jl.
Deterministic output (no random number generation).
"""
function diffuse_xyz_image_kernel(df::DataFrame; sigma_t_mm=0.6, sigma_l_mm=0.6, nbins=100, nsigma=3.0,
                                   ef::Float64=1e+5/2.5)

    # Extract hit positions and electron counts
    x_hits = Float64.(df.x)
    y_hits = Float64.(df.y)
    z_hits = Float64.(df.z)

    # Check if electrons column exists, otherwise use energy
    if hasproperty(df, :electrons)
        electrons = Float64.(df.electrons)
    elseif hasproperty(df, :energy)
        electrons = df.energy .* ef
    else
        @warn "DataFrame does not have 'electrons' column, using 1 electron per hit"
        electrons = ones(Float64, length(x_hits))
    end

    n_hits = length(x_hits)
    if n_hits == 0
        event_id = nrow(df) > 0 ? df.event_id[1] : 0
        return DataFrame(
            event_id = Int[],
            x = Float64[],
            y = Float64[],
            z = Float64[],
            energy = Float64[],
            electrons = Int[]
        )
    end

    # Handle sigma_l = 0 case (no longitudinal diffusion, e.g., for ions)
    no_z_diffusion = (sigma_l_mm <= 0.0)
    sigma_l_eff = no_z_diffusion ? sigma_t_mm : sigma_l_mm  # Use sigma_t for range calculation

    # Calculate histogram range with padding
    x_min = minimum(x_hits) - nsigma * sigma_t_mm
    x_max = maximum(x_hits) + nsigma * sigma_t_mm
    y_min = minimum(y_hits) - nsigma * sigma_t_mm
    y_max = maximum(y_hits) + nsigma * sigma_t_mm
    z_min = minimum(z_hits) - nsigma * sigma_l_eff
    z_max = maximum(z_hits) + nsigma * sigma_l_eff

    # Create grid edges and compute bin widths
    x_edges = range(x_min, x_max, length=nbins+1)
    y_edges = range(y_min, y_max, length=nbins+1)
    z_edges = range(z_min, z_max, length=nbins+1)

    dx = x_edges[2] - x_edges[1]
    dy = y_edges[2] - y_edges[1]
    dz = z_edges[2] - z_edges[1]

    # Bin centers as vectors
    x_centers = collect((x_edges[1:end-1] .+ x_edges[2:end]) ./ 2)
    y_centers = collect((y_edges[1:end-1] .+ y_edges[2:end]) ./ 2)
    z_centers = collect((z_edges[1:end-1] .+ z_edges[2:end]) ./ 2)

    # Bin hits into 3D histogram
    hist = zeros(Float64, nbins, nbins, nbins)
    @inbounds for i in 1:n_hits
        ix = clamp(Int(floor((x_hits[i] - x_min) / dx)) + 1, 1, nbins)
        iy = clamp(Int(floor((y_hits[i] - y_min) / dy)) + 1, 1, nbins)
        iz = clamp(Int(floor((z_hits[i] - z_min) / dz)) + 1, 1, nbins)
        hist[ix, iy, iz] += electrons[i]
    end

    # Build Gaussian kernel in bin units and apply filtering
    σx_bins = sigma_t_mm / dx
    σy_bins = sigma_t_mm / dy

    if no_z_diffusion
        # 2D blur only in x,y - no z diffusion for ions
        kernel = KernelFactors.gaussian((σx_bins, σy_bins, 0.0))
    else
        σz_bins = sigma_l_mm / dz
        kernel = KernelFactors.gaussian((σx_bins, σy_bins, σz_bins))
    end

    # Apply 3D Gaussian filtering
    blurred = imfilter(hist, kernel, Fill(0.0))

    # Count non-zero voxels first to preallocate
    n_nonzero = 0
    @inbounds for ix in 1:nbins, iy in 1:nbins, iz in 1:nbins
        if blurred[ix, iy, iz] > 0.5
            n_nonzero += 1
        end
    end

    # Preallocate output arrays
    result_x = Vector{Float64}(undef, n_nonzero)
    result_y = Vector{Float64}(undef, n_nonzero)
    result_z = Vector{Float64}(undef, n_nonzero)
    result_electrons = Vector{Int}(undef, n_nonzero)

    # Fill output arrays
    idx = 1
    @inbounds for ix in 1:nbins
        for iy in 1:nbins
            for iz in 1:nbins
                count = blurred[ix, iy, iz]
                if count > 0.5
                    result_x[idx] = x_centers[ix]
                    result_y[idx] = y_centers[iy]
                    result_z[idx] = z_centers[iz]
                    result_electrons[idx] = round(Int, count)
                    idx += 1
                end
            end
        end
    end

    # Get event_id
    event_id = nrow(df) > 0 ? df.event_id[1] : 0

    return DataFrame(
        event_id = event_id,
        x = result_x,
        y = result_y,
        z = result_z,
        energy = result_electrons / ef,
        electrons = result_electrons
    )
end
=#


"""
    select_events_itaca(hitsdf, nevent; kwargs...)

Process event with ITACA diffusion: transform hits, diffuse, voxelize, build tracks.
Returns vector of Tracks (empty if energy outside [emin, emax]).
"""
function select_events_itaca(hitsdf::DataFrame, nevent::Int;
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

    tK = 293.0 # 20 C
    Ed = 500.0 # V/cm

    # get mc event
    event_data = get_event(hitsdf, nevent)
    evtmc = transform_hits_df(event_data)

    # Calculate total event energy in keV
    energy = 1e+3 * sum(event_data.energy)
   
    # Check if event energy is within range
    if energy < emin || energy > emax
        #println("Event $nevent rejected: energy=$energy keV not in range [$emin, $emax] keV")
        return (Tracks[], Tracks[], Tracks[]) # Return empty tracks arrays
    end
    
    # decide the L of the event
    levt = rand(Uniform(lmin + lbuff, lmax - lbuff))
    lion = lmax - levt # position ion
    lele = levt # position electron

    # MC tracks: voxelize with 1.0 mm, no diffusion
    mc_voxel_size = 1.0
    mc_max_distance = 3.0
    vmct = voxelize_event(evtmc, mc_voxel_size)
    mc_diffusion = DiffusionParams(0.0, 0.0, 0.0, mc_voxel_size, mc_max_distance,
                                   0.0, nbins_df, nsigma_df)
    mc_tracks = make_tracks(vmct;
                            max_distance_mm=mc_max_distance,
                            energy_threshold_kev=0.0,
                            diffusion=mc_diffusion)

    # Ion tracks
    sigma_t_ion = sigma_t_ion_mm(tK, lion, Ed) # 20 C 500 V/cm
    sigma_l_ion = sigma_t_ion # approximation!
    voxel_size_ion = voxel_scale * sigma_t_ion
    max_distance_ion = voxel_dd * sigma_t_ion

    iondf = diffuse_xyz_image_mc(evtmc;
                                 sigma_t_mm=sigma_t_ion,
                                 sigma_l_mm=sigma_l_ion,
                                 nbins=nbins_df,
                                 nsigma=nsigma_df)
    vion = voxelize_event(iondf, voxel_size_ion)

    ion_diffusion = DiffusionParams(lion, sigma_t_ion, sigma_l_ion, voxel_size_ion,
                                    max_distance_ion, energy_threshold_kev, nbins_df, nsigma_df)
    ion_tracks = make_tracks(vion;
                             max_distance_mm=max_distance_ion,
                             energy_threshold_kev=energy_threshold_kev,
                             diffusion=ion_diffusion)

    # Electron tracks
    sigma_t_ele = sigma_t_mm(lele, pbar; dtmm=Dt)
    sigma_l_ele = sigma_l_mm(lele, pbar; dlmm=Dl)
    voxel_size_ele = voxel_scale * sigma_t_ele
    max_distance_ele = voxel_dd * sigma_t_ele

    eledf = diffuse_xyz_image_mc(evtmc;
                                 sigma_t_mm=sigma_t_ele,
                                 sigma_l_mm=sigma_l_ele,
                                 nbins=nbins_df,
                                 nsigma=nsigma_df)
    vele = voxelize_event(eledf, voxel_size_ele)

    ele_diffusion = DiffusionParams(lele, sigma_t_ele, sigma_l_ele, voxel_size_ele,
                                    max_distance_ele, energy_threshold_kev, nbins_df, nsigma_df)
    ele_tracks = make_tracks(vele;
                             max_distance_mm=max_distance_ele,
                             energy_threshold_kev=energy_threshold_kev,
                             diffusion=ele_diffusion)

    return (mc_tracks, ion_tracks, ele_tracks)
end


# ============================================================
# Helper functions for ITACA analysis scripts
# (Previously in itaca_aux.jl)
# ============================================================

"""
    get_sigma(particle_type, ldrft; dt, dl, tK, edrift, Pbar)

Compute transverse and longitudinal sigma based on particle type.

# Arguments
- `particle_type`: "ion" or "electron"
- `ldrft`: Drift length in cm
- `dt`: Transverse diffusion coefficient mm/√cm (default: 3.5)
- `dl`: Longitudinal diffusion coefficient mm/√cm (default: 0.9)
- `tK`: Temperature in Kelvin (default: 297.0)
- `edrift`: Drift field in V/cm (default: 500.0)
- `Pbar`: Pressure in bar (default: 15.0)

# Returns
- `(σt, σl)`: Tuple of transverse and longitudinal sigma in mm
"""
function get_sigma(particle_type, ldrft;
                   dt=3.5, dl=0.9,
                   tK=297.0, edrift=500.0, Pbar=15.0)
    if particle_type == "ion"
        σt = sigma_t_ion_mm(tK, ldrft, edrift)
        σl = 0.0
    else
        σt = sigma_t_mm(ldrft, Pbar; dtmm=dt)
        σl = sigma_l_mm(ldrft, Pbar; dlmm=dl)
    end
    σt, σl
end


"""
    get_energy_threshold(particle_type; energy_threshold_ions, energy_threshold_keV)

Get energy threshold based on particle type.

# Arguments
- `particle_type`: "ion" or "electron"
- `energy_threshold_ions`: Threshold for ions (default: 10.0)
- `energy_threshold_keV`: Threshold in keV for electrons (default: 10.0)

# Returns
- Energy threshold in keV
"""
function get_energy_threshold(particle_type;
                              energy_threshold_ions=10.0,
                              energy_threshold_keV=10.0)
    f = 1e+5/2.5  # ions per MeV
    fkeV = f*1e-3  # ions per keV

    if particle_type == "ion"
        energy_threshold_keV = energy_threshold_ions/fkeV
    end
    energy_threshold_keV
end


"""
    get_voxel_size_and_distance(ldrft, σt)

Compute voxel size, MC voxel size, and max distance based on diffusion.

# Arguments
- `ldrft`: Drift length in cm
- `σt`: Transverse sigma in mm

# Returns
- `(voxel_size, mcvox_size, max_distance)`: Tuple of sizes in mm
"""
function get_voxel_size_and_distance(ldrft, σt)
    if ldrft > 50.0
        voxel_scale = 1.5
        voxel_distance_scale = 1.5
    else
        voxel_scale = 3.0
        voxel_distance_scale = 2.0
    end

    voxel_size = σt * voxel_scale
    mcvox_size = 0.5
    max_distance = voxel_size * voxel_distance_scale
    (voxel_size, mcvox_size, max_distance)
end


"""
    length_and_energy_of_tracks(tracks)

Get length (number of voxels) and energy (keV) for each track.

# Returns
- `(LT, E)`: Tuple of vectors with lengths and energies
"""
function length_and_energy_of_tracks(tracks)
    LT = [length(track.voxels.energy) for track in tracks]
    E = [sum(track.voxels.energy)*1e+3 for track in tracks]
    LT, E
end


"""
    kde_peaks(peaks, kde_f)

Extract peak1 (leftmost) and peak2 (rightmost) from KDE peaks.
Prominence is normalized by the KDE range.

# Arguments
- `peaks`: Result from find_peaks()
- `kde_f`: KDE values (for normalization)

# Returns
NamedTuple with:
- `peak1_left`, `peak1_right`, `peak1_prom`: Leftmost peak info
- `peak2_left`, `peak2_right`, `peak2_prom`: Rightmost peak info

Values are 0.0 if peak not found.
"""
function kde_peaks(peaks, kde_f)
    # Default values (0 = not found)
    peak1_left, peak1_right, peak1_prom = 0.0, 0.0, 0.0
    peak2_left, peak2_right, peak2_prom = 0.0, 0.0, 0.0

    f_range = maximum(kde_f) - minimum(kde_f)

    if length(peaks.indices) >= 1
        # Find leftmost peak (peak1)
        left_idx = argmin(peaks.positions)
        peak1_left = peaks.lefts[left_idx]
        peak1_right = peaks.rights[left_idx]
        peak1_prom = f_range > 0 ? peaks.proms[left_idx] / f_range : 0.0

        if length(peaks.indices) >= 2
            # Find rightmost peak (peak2)
            right_idx = argmax(peaks.positions)
            peak2_left = peaks.lefts[right_idx]
            peak2_right = peaks.rights[right_idx]
            peak2_prom = f_range > 0 ? peaks.proms[right_idx] / f_range : 0.0
        end
    end

    (peak1_left=peak1_left, peak1_right=peak1_right, peak1_prom=peak1_prom,
     peak2_left=peak2_left, peak2_right=peak2_right, peak2_prom=peak2_prom)
end
