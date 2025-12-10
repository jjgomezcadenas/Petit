#!/usr/bin/env julia

"""
Check MC hits to understand true extremes.
For electrons: first hit of particle_id=1 is the start,
Bragg peak is where the primary electron ends.
"""

const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using DataFrames
using Statistics

include(joinpath(pdir, "src", "Petit.jl"))
import .Petit

# Load MC file
mc_file = ARGS[1]
println("Loading MC file: $mc_file")

dfs = Petit.get_dataset_dfs(mc_file)
hits = dfs["hits"]

event_ids = unique(hits.event_id)
println("Number of events: $(length(event_ids))")

# Analyze several events
println("\n" * "="^80)
println("ANALYSIS: Finding true MC extremes")
println("="^80)

for eid in event_ids[1:min(10, length(event_ids))]
    evt = Petit.get_event(hits, eid)

    # Get primary particle hits (particle_id = 1)
    primary = filter(row -> row.particle_id == 1, evt)

    if nrow(primary) == 0
        println("\nEvent $eid: No primary particle hits!")
        continue
    end

    # Sort primary by time
    sort!(primary, :time)

    # First and last hit of primary particle
    first_hit = first(primary)
    last_hit = last(primary)

    # Also check: last hit by time across ALL particles
    evt_sorted = sort(evt, :time)
    last_hit_overall = last(evt_sorted)

    println("\n--- Event $eid ---")
    println("Total hits: $(nrow(evt)), Primary hits: $(nrow(primary))")
    println("Primary particle time range: $(first_hit.time) to $(last_hit.time)")

    println("\nTRUE EXTREME 1 (start, first hit of primary):")
    println("  Position: ($(round(first_hit.x, digits=2)), $(round(first_hit.y, digits=2)), $(round(first_hit.z, digits=2)))")
    println("  Time: $(first_hit.time), Energy: $(round(first_hit.energy*1e6, digits=1)) eV")

    println("\nTRUE EXTREME 2 (Bragg peak, last hit of primary):")
    println("  Position: ($(round(last_hit.x, digits=2)), $(round(last_hit.y, digits=2)), $(round(last_hit.z, digits=2)))")
    println("  Time: $(last_hit.time), Energy: $(round(last_hit.energy*1e6, digits=1)) eV")

    # Track length (distance between extremes)
    d = sqrt((first_hit.x - last_hit.x)^2 + (first_hit.y - last_hit.y)^2 + (first_hit.z - last_hit.z)^2)
    println("\nTrack length (MC extremes distance): $(round(d, digits=1)) mm")

    # Check energy near the two extremes
    # Energy in sphere of radius 5mm around each extreme
    r = 5.0
    energy_near_start = sum(evt.energy[sqrt.((evt.x .- first_hit.x).^2 .+
                                              (evt.y .- first_hit.y).^2 .+
                                              (evt.z .- first_hit.z).^2) .< r])
    energy_near_end = sum(evt.energy[sqrt.((evt.x .- last_hit.x).^2 .+
                                            (evt.y .- last_hit.y).^2 .+
                                            (evt.z .- last_hit.z).^2) .< r])

    println("\nEnergy in $(r)mm sphere:")
    println("  Near start: $(round(energy_near_start * 1e3, digits=1)) keV")
    println("  Near end:   $(round(energy_near_end * 1e3, digits=1)) keV")
    println("  Ratio (end/start): $(round(energy_near_end/energy_near_start, digits=2))")
end
