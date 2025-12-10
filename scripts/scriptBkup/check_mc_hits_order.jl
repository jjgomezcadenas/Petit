#!/usr/bin/env julia

"""
Check if MC hits are ordered by time/sequence.
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

println("\nHits DataFrame columns: $(names(hits))")
println("Number of rows: $(nrow(hits))")

# Check first event
event_ids = unique(hits.event_id)
println("\nNumber of events: $(length(event_ids))")

# Look at first few events
for eid in event_ids[1:min(3, length(event_ids))]
    evt = Petit.get_event(hits, eid)
    println("\n" * "="^60)
    println("Event $eid: $(nrow(evt)) hits")
    println("="^60)

    # Show first and last few rows
    println("\nFirst 5 hits:")
    println(first(evt, 5))

    println("\nLast 5 hits:")
    println(last(evt, 5))

    # Check if time column exists and is ordered
    if :time in propertynames(evt)
        times = evt.time
        is_sorted = issorted(times)
        println("\nTime column exists: YES")
        println("Time range: $(minimum(times)) to $(maximum(times))")
        println("Is sorted by time: $is_sorted")

        # Check hit_id if exists
        if :hit_id in propertynames(evt)
            hit_ids = evt.hit_id
            println("Hit IDs range: $(minimum(hit_ids)) to $(maximum(hit_ids))")
            println("Is sorted by hit_id: $(issorted(hit_ids))")
        end
    else
        println("\nTime column exists: NO")
    end

    # Show positions of first and last hits
    println("\nFirst hit position: ($(evt.x[1]), $(evt.y[1]), $(evt.z[1]))")
    println("Last hit position:  ($(evt.x[end]), $(evt.y[end]), $(evt.z[end]))")

    # Distance between first and last
    d = sqrt((evt.x[1] - evt.x[end])^2 + (evt.y[1] - evt.y[end])^2 + (evt.z[1] - evt.z[end])^2)
    println("Distance first-last: $(round(d, digits=2)) mm")
end
