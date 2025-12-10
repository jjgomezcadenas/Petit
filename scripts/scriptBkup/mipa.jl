#!/usr/bin/env julia
# focusing_solver.jl
# Solve Laplace in a periodic 2D unit cell (x–z), trace field lines from a hole to a disk,
# and report the capture fraction. Also saves a PNG with potential + sample streamlines.

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using Random
using Plots
using Plots.Measures
using Printf

# ----------------------------
# Simple named-argument parser (--key=value or --key value)
# ----------------------------
function parse_args(defaults::Dict{String,Any})
    args = copy(ARGS)
    conf = copy(defaults)
    i = 1
    while i <= length(args)
        tok = args[i]
        if startswith(tok, "--")
            kv = split(tok[3:end], "=", limit=2)
            key = kv[1]
            if length(kv) == 2
                val = kv[2]
                conf[key] = val
                i += 1
            else
                # expect next token is the value
                if i == length(args)
                    error("Missing value for --$key")
                end
                conf[key] = args[i+1]
                i += 2
            end
        else
            i += 1
        end
    end
    return conf
end

function parse_float!(conf::Dict{String,Any}, key::String)
    conf[key] = parse(Float64, string(conf[key]))
end
function parse_int!(conf::Dict{String,Any}, key::String)
    conf[key] = parse(Int, string(conf[key]))
end

# ----------------------------
# Geometry + solver
# ----------------------------
function main()
    # Defaults (SI)
    defaults = Dict{String,Any}(
        "pitch" => 4.0e-4,     # cell width s (m)
        "l1"    => 5.0e-4,     # plane1->plane2 gap (m)
        "l2"    => 5.0e-4,     # plane2->plane3 gap (m)
        "t"     => 5.0e-5,     # disk thickness (m)
        "d1"    => 1.5e-4,     # hole diameter in plane 2 (m)
        "d2"    => 1.0e-4,     # disk diameter (m)
        "V0"    => 500.0,      # plane 1 potential (V)
        "dV1"   => 200.0,      # drop to plane 2
        "dV2"   => 150.0,      # extra drop (disk is at V0-dV1-dV2)
        "dV3"   => 150.0,      # extra drop to plane 3
        "Nx"    => 401,        # grid x
        "Nz"    => 401,        # grid z
        "nseeds"=> 2000,       # number of field lines to trace
        "stream_plot" => 200,  # streamlines to draw (subset)
        "outfile" => "fieldmap.png",
        "omega" => 1.8,
        "tol"   => 1e-6,
        "maxit" => 50000
    )

    conf = parse_args(defaults)
    # Parse numerics
    for k in ["pitch","l1","l2","t","d1","d2","V0","dV1","dV2","dV3","omega","tol"]
        parse_float!(conf, k)
    end
    for k in ["Nx","Nz","nseeds","stream_plot","maxit"]
        parse_int!(conf, k)
    end
    outfile = string(conf["outfile"])

    # Unpack
    s   = conf["pitch"]
    l1  = conf["l1"]
    l2  = conf["l2"]
    t   = conf["t"]
    d1  = conf["d1"]
    d2  = conf["d2"]
    V0  = conf["V0"]
    dV1 = conf["dV1"]
    dV2 = conf["dV2"]
    dV3 = conf["dV3"]
    Nx  = conf["Nx"]
    Nz  = conf["Nz"]
    nseeds = conf["nseeds"]
    nplot  = conf["stream_plot"]
    ω    = conf["omega"]
    tol  = conf["tol"]
    maxit = conf["maxit"]

    # Potentials
    V_plane1 = V0
    V_plane2 = V0 - dV1
    V_disk   = V0 - dV1 - dV2
    V_plane3 = V0 - dV1 - dV2 - dV3

    # Domain
    Lx = s
    Lz = l1 + l2 + t
    dx = Lx/(Nx-1)
    dz = Lz/(Nz-1)

    x(i) = -Lx/2 + (i-1)*dx
    z(k) = (k-1)*dz

    r_hole = d1/2
    r_disk = d2/2
    xc = 0.0
    zc_plane1 = 0.0
    zc_plane2 = l1
    zc_plane3 = l1 + l2

    # Masks
    plane1 = falses(Nx, Nz)
    plane2 = falses(Nx, Nz)
    plane3 = falses(Nx, Nz)
    disk   = falses(Nx, Nz)

    # plane 1 at z~0 (bottom row)
    for i in 1:Nx
        plane1[i,1] = true
    end
    # plane 2 at z=l1, except central hole of d1
    k2 = clamp(round(Int, zc_plane2/dz) + 1, 1, Nz)
    for i in 1:Nx
        r = abs(x(i) - xc)
        if r > r_hole + dx/2
            plane2[i,k2] = true
        end
    end
    # plane 3 at z=l1+l2 (full strip)
    k3 = clamp(round(Int, zc_plane3/dz) + 1, 1, Nz)
    for i in 1:Nx
        plane3[i,k3] = true
    end
    # disk volume (cylinder) above plane 3 with height t
    k3_top = clamp(round(Int, (zc_plane3 + t)/dz) + 1, 1, Nz)
    for i in 1:Nx, k in k3:k3_top
        r = abs(x(i) - xc)
        if r ≤ r_disk + dx/2
            disk[i,k] = true
        end
    end

    is_dirichlet = plane1 .| plane2 .| plane3 .| disk

    # Potential init
    V = zeros(Float64, Nx, Nz)
    function enforce_dirichlet!(V)
        for i in 1:Nx, k in 1:Nz
            if plane1[i,k]; V[i,k] = V_plane1; end
            if plane2[i,k]; V[i,k] = V_plane2; end
            if plane3[i,k]; V[i,k] = V_plane3; end
            if disk[i,k];   V[i,k] = V_disk;   end
        end
    end
    enforce_dirichlet!(V)

    # Periodic helpers
    idxp(i,N) = (i==N) ? 1 : i+1
    idxm(i,N) = (i==1) ? N : i-1

    # SOR Laplace solver with progress reporting
    println("Starting SOR solver...")
    println("Grid: $(Nx)×$(Nz), Max iterations: $(maxit), Tolerance: $(tol)")

    # Precompute constants for efficiency
    dx2_inv = 1.0 / dx^2
    dz2_inv = 1.0 / dz^2
    denom_inv = 1.0 / (2*dx2_inv + 2*dz2_inv)

    report_interval = max(1, maxit ÷ 20)  # Report ~20 times during solve

    converged = false
    for it in 1:maxit
        maxΔ = 0.0

        @inbounds for k in 2:Nz-1
            for i in 1:Nx
                if is_dirichlet[i,k]; continue; end
                ip = idxp(i,Nx); im = idxm(i,Nx)
                kp = k+1; km = k-1
                Vnew = ((V[ip,k] + V[im,k]) * dx2_inv +
                        (V[i,kp] + V[i,km]) * dz2_inv) * denom_inv
                δ = Vnew - V[i,k]
                V[i,k] += ω*δ
                maxΔ = max(maxΔ, abs(δ))
            end
        end
        enforce_dirichlet!(V)

        # Progress reporting
        if it % report_interval == 0 || maxΔ < tol
            @printf("  Iteration %5d/%d: residual = %.2e\n", it, maxit, maxΔ)
        end

        if maxΔ < tol
            println("✓ Converged in $it iterations (residual = $(@sprintf("%.2e", maxΔ)))")
            converged = true
            break
        end
        if it == maxit
            @warn "Reached max iterations; residual ~ $(@sprintf("%.2e", maxΔ))"
        end
    end

    if !converged
        println("⚠ Did not converge to tolerance $(tol)")
    end

    # Field computation
    println("Computing electric field...")
    Ex = zeros(Float64, Nx, Nz)
    Ez = zeros(Float64, Nx, Nz)
    dx2_inv_field = 1.0 / (2dx)
    dz2_inv_field = 1.0 / (2dz)

    @inbounds for k in 2:Nz-1
        for i in 1:Nx
            ip = idxp(i,Nx); im = idxm(i,Nx)
            Ex[i,k] = -(V[ip,k] - V[im,k]) * dx2_inv_field   # periodic x
            Ez[i,k] = -(V[i,k+1] - V[i,k-1]) * dz2_inv_field
        end
    end

    # Interp E (bilinear), periodic in x
    function interp_E(xq::Float64, zq::Float64)
        ξ = (xq + Lx/2)/dx + 1
        ζ = zq/dz + 1
        i = clamp(floor(Int, ξ), 1, Nx-1)
        k = clamp(floor(Int, ζ), 1, Nz-1)
        tx = clamp(ξ - i, 0.0, 1.0)
        tz = clamp(ζ - k, 0.0, 1.0)
        ip = idxp(i,Nx)
        kp = k+1

        ex = (1-tx)*(1-tz)*Ex[i,k] + tx*(1-tz)*Ex[ip,k] +
             (1-tx)*tz*Ex[i,kp] + tx*tz*Ex[ip,kp]
        ez = (1-tx)*(1-tz)*Ez[i,k] + tx*(1-tz)*Ez[ip,k] +
             (1-tx)*tz*Ez[i,kp] + tx*tz*Ez[ip,kp]
        return ex, ez
    end

    in_disk(i,k) = disk[i,k]
    in_conductor(i,k) = is_dirichlet[i,k] && !plane2[i,k]  # plane2 hole not a solid

    # RK4 field-line tracer
    function trace_line(x0,z0; h=0.25*min(dx,dz), maxsteps=20000)
        x, z = x0, z0
        hit_disk = false
        xs = Float64[]; zs = Float64[]   # for plotting
        for _ in 1:maxsteps
            if z ≤ 0 || z ≥ Lz; break; end
            push!(xs, x); push!(zs, z)
            i = clamp(round(Int, (x + Lx/2)/dx) + 1, 1, Nx)
            k = clamp(round(Int,  z/dz) + 1, 1, Nz)
            if in_conductor(i,k)
                hit_disk = in_disk(i,k)
                break
            end
            ex1, ez1 = interp_E(x, z)
            n1 = hypot(ex1, ez1) + 1e-30
            vx1, vz1 = ex1/n1, ez1/n1

            ex2, ez2 = interp_E(x + 0.5h*vx1, z + 0.5h*vz1)
            n2 = hypot(ex2, ez2) + 1e-30
            vx2, vz2 = ex2/n2, ez2/n2

            ex3, ez3 = interp_E(x + 0.5h*vx2, z + 0.5h*vz2)
            n3 = hypot(ex3, ez3) + 1e-30
            vx3, vz3 = ex3/n3, ez3/n3

            ex4, ez4 = interp_E(x + h*vx3, z + h*vz3)
            n4 = hypot(ex4, ez4) + 1e-30
            vx4, vz4 = ex4/n4, ez4/n4

            vx = (vx1 + 2vx2 + 2vx3 + vx4)/6
            vz = (vz1 + 2vz2 + 2vz3 + vz4)/6

            x = x + h*vx
            if x < -Lx/2; x += Lx; end
            if x >  Lx/2; x -= Lx; end
            z = z + h*vz
        end
        return hit_disk, xs, zs
    end

    # Seed uniformly in the hole, just above plane 2
    function capture_fraction(nseeds::Int; want_paths::Bool=false, nplot::Int=200)
        zseed = zc_plane2 + 1.5dz
        hits = 0
        seeds = 0
        paths = Tuple{Vector{Float64},Vector{Float64}}[]
        while seeds < nseeds
            u, v = rand(), rand()
            r = r_hole*sqrt(u)
            θ = 2π*v
            xs = r*cos(θ)
            if abs(xs) ≤ Lx/2
                seeds += 1
                hit, xp, zp = trace_line(xs, zseed)
                hits += hit ? 1 : 0
                if want_paths && (length(paths) < nplot)
                    push!(paths, (xp,zp))
                end
            end
        end
        frac = hits / nseeds
        return frac, paths
    end

    # Compute capture fraction and collect a sample of streamlines
    println("\nTracing $(nseeds) field lines...")
    frac, paths = capture_fraction(nseeds; want_paths=true, nplot=nplot)
    println("✓ Capture fraction to disk: ", round(frac*100; digits=2), " %")

    # Plot potential + streamlines
    println("\nGenerating visualization...")
    plt = heatmap(
        -Lx/2:dx:Lx/2, 0:dz:Lz, V'; aspect_ratio = Lz/Lx,
        xlabel = "x (m)", ylabel = "z (m)",
        cbar = true, title = "Potential and field lines",
        color = :viridis, right_margin = 8mm, bottom_margin = 5mm
    )

    # Draw streamlines
    for (xp,zp) in paths
        plot!(plt, xp, zp, color=:white, lw=0.6, alpha=0.8, label=false)
    end

    # Outline disk (top view in x–z; draw its lateral projection as vertical bar)
    xs_disk = [-r_disk, -r_disk,  r_disk,  r_disk, -r_disk]
    zs_disk = [zc_plane3, zc_plane3 + t, zc_plane3 + t, zc_plane3, zc_plane3]
    plot!(plt, xs_disk, zs_disk, color=:red, lw=1.2, label=false)

    # Draw plane lines (thin) + hole edges on plane 2
    plot!(plt, [-Lx/2, -r_hole], [zc_plane2, zc_plane2], color=:gray, lw=1.0, label=false)
    plot!(plt, [ r_hole,  Lx/2], [zc_plane2, zc_plane2], color=:gray, lw=1.0, label=false)
    plot!(plt, [-Lx/2, Lx/2], [0.0, 0.0], color=:gray, lw=1.0, label=false)
    plot!(plt, [-Lx/2, Lx/2], [zc_plane3, zc_plane3], color=:gray, lw=1.0, label=false)

    savefig(plt, outfile)
    println("✓ Saved: $outfile")
    println("\n" * "="^60)
    println("SIMULATION COMPLETE")
    println("="^60)
end

main()