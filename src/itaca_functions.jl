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

