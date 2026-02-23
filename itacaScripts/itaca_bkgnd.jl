#!/usr/bin/env julia

"""
ITACA Background Budget Calculator.

Computes the background in the ROI (2400-2500 keV) for Bi-214 and Tl-208 from
detector components. Applies three efficiency factors:

1. Geometry cut efficiency (ε_geom): Fraction of events in fiducial volume
   - Source: datainfo/fiducial_efficiency.csv (from itaca_fiducial.jl)

2. Energy cut efficiency (ε_ecut): Fraction of fiducial events in ROI
   - Source: datainfo/fiducial_efficiency.csv (from itaca_fiducial.jl)

3. Single-track efficiency (ε_1trk): Fraction of events reconstructed as single track
   - Source: tables/efficiency_table.csv (from itaca_tables.jl)
   - Computed as: n_single / n_events from statistics.json files

The total fiducial efficiency is ε_fid = ε_geom × ε_ecut.
Final background = raw_events × ε_fid × ε_1trk

Usage:
    julia itaca_bkgnd.jl [options]

Options:
    --fiducial=X        Path to fiducial efficiency CSV (default: datainfo/fiducial_efficiency.csv)
    --efficiency-table=X  Path to single-track efficiency CSV (default: tables/efficiency_table.csv)
    --output-dir=X      Output directory for tables (default: tables/)
    --output-prefix=X   Prefix for output files (default: background_tables)
    --help, -h          Show this help message

Output:
    For each sigma configuration (σ=1mm, σ=10mm):
    - Console output with efficiency tables and event counts
    - CSV file with full data
    - Markdown file with formatted tables
    - LaTeX file with publication-ready tables (separate tables for Bi-214 and Tl-208)
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using Printf
using CSV
using DataFrames
using Dates

# =============================================================================
# Constants
# =============================================================================

# Seconds per year
const SECONDS_PER_YEAR = 365.25 * 24 * 3600  # 3.156e7 s

# Default path to fiducial efficiency CSV (relative to script directory)
const DEFAULT_FIDUCIAL_CSV = joinpath(dirname(@__FILE__), "datainfo", "fiducial_efficiency.csv")

# Default path to single-track efficiency table CSV (relative to script directory)
const DEFAULT_EFFICIENCY_TABLE = joinpath(dirname(@__FILE__), "tables", "efficiency_table.csv")

# Default output directory for tables (relative to script directory)
const DEFAULT_OUTPUT_DIR = joinpath(dirname(@__FILE__), "tables")

# Default output prefix for background tables
const DEFAULT_OUTPUT_PREFIX = "background_tables"

# Default path to 0νββ signal info file (relative to script directory)
const DEFAULT_SIGNAL_INFO = joinpath(dirname(@__FILE__), "datainfo", "0nubb_filtered_ecut_2400_2500_info.csv")

# Default 0νββ geometry efficiency (fraction of events in fiducial volume)
const DEFAULT_SIGNAL_GEOM_EFF = 0.95

# =============================================================================
# Data Structures
# =============================================================================

"""
Component in the ITACA detector with radioactivity data.
"""
struct DetectorComponent
    name::String
    material::String
    mass_kg::Float64
    activity_bi214_mBq::Float64
    activity_tl208_mBq::Float64
    zone::Symbol  # :copper_shell, :ptfe_barrel, or :ptfe_endcap
end

"""
Fiducial efficiency for a simulation zone (from itaca_fiducial.jl results).

Contains three efficiencies per isotope:
- eff_geom: geometry cut efficiency (events in fiducial region / total events)
- eff_ecut: energy cut efficiency (events in ROI / events in fiducial region)
- eff_fid: total fiducial efficiency (product of geometry and energy cut)
"""
struct FiducialEfficiency
    zone::Symbol
    eff_geom_bi214::Float64
    eff_geom_tl208::Float64
    eff_ecut_bi214::Float64
    eff_ecut_tl208::Float64
    eff_fid_bi214::Float64
    eff_fid_tl208::Float64
end

"""
Single-track efficiency for a simulation zone (from itaca_tables.jl results).
"""
struct SingleTrackEfficiency
    zone::Symbol
    eff_bi214::Float64
    eff_tl208::Float64
end

# =============================================================================
# Fiducial Efficiency Loading
# =============================================================================

"""
    load_fiducial_efficiencies(csv_path::String)

Load fiducial efficiencies from a CSV file produced by itaca_fiducial.jl.
The CSV must have columns:
- source: source name (e.g., bi214_copper_shell)
- fiducial_geom_efficiency: geometry cut efficiency
- fiducial_cut: energy cut efficiency
- fiducial_tot: total fiducial efficiency

Source names are parsed to extract isotope and zone:
- bi214_copper_shell → :copper_shell, Bi-214
- bi214_ptfe_barrel* → :ptfe_barrel, Bi-214
- bi214_ptfe_endcap* → :ptfe_endcap, Bi-214
- tl208_copper_shell → :copper_shell, Tl-208
- tl208_ptfe_barrel* → :ptfe_barrel, Tl-208
- tl208_ptfe_endcap* → :ptfe_endcap, Tl-208

Returns a Dict{Symbol, FiducialEfficiency}.
"""
function load_fiducial_efficiencies(csv_path::String)
    if !isfile(csv_path)
        error("Fiducial efficiency CSV not found: $csv_path")
    end

    df = CSV.read(csv_path, DataFrame; comment="#")

    # Initialize storage for efficiencies by zone
    # Each zone stores Dict of isotope → (eff_geom, eff_ecut, eff_fid)
    zone_effs = Dict{Symbol, Dict{Symbol, NamedTuple{(:geom, :ecut, :fid), Tuple{Float64, Float64, Float64}}}}()
    for zone in [:copper_shell, :ptfe_barrel, :ptfe_endcap]
        zone_effs[zone] = Dict{Symbol, NamedTuple{(:geom, :ecut, :fid), Tuple{Float64, Float64, Float64}}}()
    end

    # Parse each row
    for row in eachrow(df)
        source = row.source
        eff_geom = row.fiducial_geom_efficiency
        eff_ecut = row.fiducial_cut
        eff_fid = row.fiducial_tot

        # Determine isotope
        isotope = nothing
        if startswith(source, "bi214")
            isotope = :bi214
        elseif startswith(source, "tl208")
            isotope = :tl208
        else
            continue  # Skip unknown isotopes
        end

        # Determine zone (look for zone name in source string)
        zone = nothing
        if occursin("copper_shell", source)
            zone = :copper_shell
        elseif occursin("ptfe_barrel", source)
            zone = :ptfe_barrel
        elseif occursin("ptfe_endcap", source)
            zone = :ptfe_endcap
        else
            continue  # Skip unknown zones
        end

        zone_effs[zone][isotope] = (geom=eff_geom, ecut=eff_ecut, fid=eff_fid)
    end

    # Build FiducialEfficiency structs
    result = Dict{Symbol, FiducialEfficiency}()
    for zone in [:copper_shell, :ptfe_barrel, :ptfe_endcap]
        effs_bi214 = get(zone_effs[zone], :bi214, (geom=0.0, ecut=0.0, fid=0.0))
        effs_tl208 = get(zone_effs[zone], :tl208, (geom=0.0, ecut=0.0, fid=0.0))
        result[zone] = FiducialEfficiency(
            zone,
            effs_bi214.geom, effs_tl208.geom,
            effs_bi214.ecut, effs_tl208.ecut,
            effs_bi214.fid, effs_tl208.fid
        )
    end

    return result
end

"""
    load_single_track_efficiencies(csv_path::String, sigma::Symbol)

Load single-track efficiencies from the efficiency table CSV.
The CSV must have columns: source, sigma_1mm, sigma_10mm, ratio

Source names are parsed to extract isotope and zone:
- "bi214 (copper_shell)" → :bi214, :copper_shell
- "tl208 (ptfe_barrel)" → :tl208, :ptfe_barrel
- "0nubb" → signal (skipped for background)

Arguments:
- csv_path: Path to efficiency_table.csv
- sigma: Either :sigma_1mm or :sigma_10mm

Returns a Dict{Symbol, SingleTrackEfficiency}.
"""
function load_single_track_efficiencies(csv_path::String, sigma::Symbol)
    if !isfile(csv_path)
        error("Single-track efficiency CSV not found: $csv_path")
    end

    df = CSV.read(csv_path, DataFrame; comment="#")

    # Determine which column to use
    eff_column = sigma == :sigma_1mm ? :sigma_1mm : :sigma_10mm

    # Initialize storage for efficiencies by zone
    zone_effs = Dict{Symbol, Dict{Symbol, Float64}}()
    for zone in [:copper_shell, :ptfe_barrel, :ptfe_endcap]
        zone_effs[zone] = Dict{Symbol, Float64}()
    end

    # Parse each row
    for row in eachrow(df)
        source = row.source
        eff = row[eff_column]

        # Skip signal (0nubb) and handle missing values
        if source == "0nubb" || ismissing(eff)
            continue
        end

        # Parse source: "isotope (zone)" format
        m = match(r"^(\w+)\s*\((\w+)\)$", source)
        if m === nothing
            continue  # Skip unparseable sources
        end

        isotope_str = m.captures[1]
        zone_str = m.captures[2]

        # Determine isotope
        isotope = nothing
        if isotope_str == "bi214"
            isotope = :bi214
        elseif isotope_str == "tl208"
            isotope = :tl208
        else
            continue  # Skip unknown isotopes
        end

        # Determine zone
        zone = nothing
        if zone_str == "copper_shell"
            zone = :copper_shell
        elseif zone_str == "ptfe_barrel"
            zone = :ptfe_barrel
        elseif zone_str == "ptfe_endcap"
            zone = :ptfe_endcap
        else
            continue  # Skip unknown zones
        end

        zone_effs[zone][isotope] = eff
    end

    # Build SingleTrackEfficiency structs
    result = Dict{Symbol, SingleTrackEfficiency}()
    for zone in [:copper_shell, :ptfe_barrel, :ptfe_endcap]
        eff_bi214 = get(zone_effs[zone], :bi214, 0.0)
        eff_tl208 = get(zone_effs[zone], :tl208, 0.0)
        result[zone] = SingleTrackEfficiency(zone, eff_bi214, eff_tl208)
    end

    return result
end

"""
    load_signal_1trk_efficiency(csv_path::String, sigma::Symbol)

Load 0νββ signal single-track efficiency from the efficiency table CSV.
Returns the efficiency value for the specified sigma.
"""
function load_signal_1trk_efficiency(csv_path::String, sigma::Symbol)
    if !isfile(csv_path)
        error("Efficiency table CSV not found: $csv_path")
    end

    df = CSV.read(csv_path, DataFrame; comment="#")

    # Determine which column to use
    eff_column = sigma == :sigma_1mm ? :sigma_1mm : :sigma_10mm

    # Find the 0nubb row
    for row in eachrow(df)
        if row.source == "0nubb"
            return row[eff_column]
        end
    end

    return 0.0  # Return 0 if not found
end

"""
    load_signal_ecut_efficiency(info_path::String)

Load 0νββ signal energy cut efficiency from the info file.
Parses the "Filter efficiency:" line from the 0nubb_filtered_ecut_2400_2500_info.csv file.
Returns the efficiency as a fraction (e.g., 0.9878 for 98.78%).
"""
function load_signal_ecut_efficiency(info_path::String)
    if !isfile(info_path)
        @warn "Signal info file not found: $info_path, using default efficiency 1.0"
        return 1.0
    end

    # Parse the info file to find "Filter efficiency:"
    for line in eachline(info_path)
        if occursin("Filter efficiency:", line)
            # Extract percentage, e.g., "  Filter efficiency:     98.78%"
            m = match(r"Filter efficiency:\s*([\d.]+)%", line)
            if m !== nothing
                return parse(Float64, m.captures[1]) / 100.0
            end
        end
    end

    @warn "Could not parse filter efficiency from $info_path, using default 1.0"
    return 1.0
end

"""
    SignalEfficiency

Container for 0νββ signal efficiency components.
"""
struct SignalEfficiency
    eff_geom::Float64    # Geometry efficiency (fraction in fiducial volume)
    eff_ecut::Float64    # Energy cut efficiency (fraction in ROI)
    eff_fid::Float64     # Total fiducial efficiency (eff_geom × eff_ecut)
    eff_1trk::Float64    # Single-track efficiency
    eff_total::Float64   # Total signal efficiency (eff_fid × eff_1trk)
end

"""
    load_signal_efficiencies(efficiency_csv::String, info_path::String, sigma::Symbol; geom_eff::Float64=0.95)

Load all 0νββ signal efficiencies and compute derived quantities.

Arguments:
- efficiency_csv: Path to efficiency_table.csv (contains ε_1trk)
- info_path: Path to 0nubb_filtered_ecut_2400_2500_info.csv (contains ε_ecut)
- sigma: Either :sigma_1mm or :sigma_10mm
- geom_eff: Geometry efficiency (default 0.95, fraction of events in fiducial volume)

Returns a SignalEfficiency struct with all efficiency components.
"""
function load_signal_efficiencies(efficiency_csv::String, info_path::String, sigma::Symbol; geom_eff::Float64=0.95)
    # Load single-track efficiency
    eff_1trk = load_signal_1trk_efficiency(efficiency_csv, sigma)

    # Load energy cut efficiency
    eff_ecut = load_signal_ecut_efficiency(info_path)

    # Compute derived efficiencies
    eff_fid = geom_eff * eff_ecut
    eff_total = eff_fid * eff_1trk

    return SignalEfficiency(geom_eff, eff_ecut, eff_fid, eff_1trk, eff_total)
end

# Zone mapping:
# - copper_shell: ICS_15cm, FC_rings (copper components)
# - ptfe_barrel: BFD (barrel field degrader)
# - ptfe_endcap: SiPM_boards, Light_Guides, FAT_GEM, Cathode, MARS_* (endcap components)

const ITACA_COMPONENTS = [
    DetectorComponent("BFD", "PTFE", 171.2, 0.0171, 0.2191, :ptfe_barrel),
    DetectorComponent("ICS_15cm", "Cu", 13699.0, 16.44, 19.18, :copper_shell),
    DetectorComponent("FC_rings", "Cu", 107.7, 0.1292, 0.1508, :copper_shell),
    DetectorComponent("Light_Guides", "HDPE", 268.6, 1.67, 2.15, :ptfe_endcap),
    DetectorComponent("FAT_GEM", "HDPE", 26.00, 0.1612, 0.2080, :ptfe_endcap),
    DetectorComponent("SiPM_boards", "Kapton", 2.28, 0.1827, 0.2512, :ptfe_endcap),
    DetectorComponent("Cathode", "Fe316Ti", 0.7954, 1.51, 0.3182, :ptfe_endcap),
    DetectorComponent("MARS_Ti", "Ti", 0.3006, 0.2796, 0.0661, :ptfe_endcap),
    DetectorComponent("MARS_Kapton", "Kapton", 0.0511, 4.09e-03, 5.62e-03, :ptfe_endcap),
    DetectorComponent("MARS_HDPE", "HDPE", 0.3458, 2.14e-03, 2.77e-03, :ptfe_endcap),
]

# =============================================================================
# Calculation Functions
# =============================================================================

"""
    activity_to_events_per_year(activity_mBq::Float64)

Convert activity in mBq to events per year.
1 mBq = 1e-3 decays/second
Events/year = activity_mBq * 1e-3 * seconds_per_year
"""
function activity_to_events_per_year(activity_mBq::Float64)
    return activity_mBq * 1e-3 * SECONDS_PER_YEAR
end

"""
    get_fiducial_efficiency(fid_effs::Dict, zone::Symbol, isotope::Symbol)

Get total fiducial efficiency for a given zone and isotope.
"""
function get_fiducial_efficiency(fid_effs::Dict, zone::Symbol, isotope::Symbol)
    fid = fid_effs[zone]
    return isotope == :bi214 ? fid.eff_fid_bi214 : fid.eff_fid_tl208
end

"""
    get_fiducial_efficiency_geom(fid_effs::Dict, zone::Symbol, isotope::Symbol)

Get geometry fiducial efficiency for a given zone and isotope.
"""
function get_fiducial_efficiency_geom(fid_effs::Dict, zone::Symbol, isotope::Symbol)
    fid = fid_effs[zone]
    return isotope == :bi214 ? fid.eff_geom_bi214 : fid.eff_geom_tl208
end

"""
    get_fiducial_efficiency_ecut(fid_effs::Dict, zone::Symbol, isotope::Symbol)

Get energy cut efficiency for a given zone and isotope.
"""
function get_fiducial_efficiency_ecut(fid_effs::Dict, zone::Symbol, isotope::Symbol)
    fid = fid_effs[zone]
    return isotope == :bi214 ? fid.eff_ecut_bi214 : fid.eff_ecut_tl208
end

"""
    get_single_track_efficiency(st_effs::Dict, zone::Symbol, isotope::Symbol)

Get single-track efficiency for a given zone and isotope.
"""
function get_single_track_efficiency(st_effs::Dict, zone::Symbol, isotope::Symbol)
    st = st_effs[zone]
    return isotope == :bi214 ? st.eff_bi214 : st.eff_tl208
end

"""
    compute_events_with_fiducial(components::Vector{DetectorComponent}, fid_effs::Dict)

Compute events per year for each component, both raw and after fiducial cut.
Returns a vector of named tuples.
"""
function compute_events_with_fiducial(components::Vector{DetectorComponent}, fid_effs::Dict)
    results = NamedTuple[]

    for c in components
        # Raw events
        events_bi214_raw = activity_to_events_per_year(c.activity_bi214_mBq)
        events_tl208_raw = activity_to_events_per_year(c.activity_tl208_mBq)

        # Fiducial efficiencies
        eff_bi214 = get_fiducial_efficiency(fid_effs, c.zone, :bi214)
        eff_tl208 = get_fiducial_efficiency(fid_effs, c.zone, :tl208)

        # Events after fiducial
        events_bi214_fid = events_bi214_raw * eff_bi214
        events_tl208_fid = events_tl208_raw * eff_tl208

        push!(results, (
            name = c.name,
            material = c.material,
            mass_kg = c.mass_kg,
            zone = c.zone,
            activity_bi214_mBq = c.activity_bi214_mBq,
            activity_tl208_mBq = c.activity_tl208_mBq,
            eff_bi214 = eff_bi214,
            eff_tl208 = eff_tl208,
            events_bi214_raw = events_bi214_raw,
            events_tl208_raw = events_tl208_raw,
            events_bi214_fid = events_bi214_fid,
            events_tl208_fid = events_tl208_fid,
            events_total_raw = events_bi214_raw + events_tl208_raw,
            events_total_fid = events_bi214_fid + events_tl208_fid
        ))
    end

    return results
end

"""
    compute_events_with_single_track(components::Vector{DetectorComponent}, fid_effs::Dict, st_effs::Dict)

Compute events per year for each component, including single-track efficiency.
Returns a vector of named tuples with raw, fiducial, and final (fiducial × single-track) events.
Includes all three fiducial efficiencies: geometry, energy cut, and total.
"""
function compute_events_with_single_track(components::Vector{DetectorComponent}, fid_effs::Dict, st_effs::Dict)
    results = NamedTuple[]

    for c in components
        # Raw events
        events_bi214_raw = activity_to_events_per_year(c.activity_bi214_mBq)
        events_tl208_raw = activity_to_events_per_year(c.activity_tl208_mBq)

        # Fiducial efficiencies (all three components)
        eff_geom_bi214 = get_fiducial_efficiency_geom(fid_effs, c.zone, :bi214)
        eff_geom_tl208 = get_fiducial_efficiency_geom(fid_effs, c.zone, :tl208)
        eff_ecut_bi214 = get_fiducial_efficiency_ecut(fid_effs, c.zone, :bi214)
        eff_ecut_tl208 = get_fiducial_efficiency_ecut(fid_effs, c.zone, :tl208)
        eff_fid_bi214 = get_fiducial_efficiency(fid_effs, c.zone, :bi214)
        eff_fid_tl208 = get_fiducial_efficiency(fid_effs, c.zone, :tl208)

        # Single-track efficiencies
        eff_1trk_bi214 = get_single_track_efficiency(st_effs, c.zone, :bi214)
        eff_1trk_tl208 = get_single_track_efficiency(st_effs, c.zone, :tl208)

        # Events after fiducial (using total fiducial efficiency)
        events_bi214_fid = events_bi214_raw * eff_fid_bi214
        events_tl208_fid = events_tl208_raw * eff_fid_tl208

        # Events after fiducial AND single-track selection
        events_bi214_final = events_bi214_fid * eff_1trk_bi214
        events_tl208_final = events_tl208_fid * eff_1trk_tl208

        push!(results, (
            name = c.name,
            material = c.material,
            mass_kg = c.mass_kg,
            zone = c.zone,
            activity_bi214_mBq = c.activity_bi214_mBq,
            activity_tl208_mBq = c.activity_tl208_mBq,
            eff_geom_bi214 = eff_geom_bi214,
            eff_geom_tl208 = eff_geom_tl208,
            eff_ecut_bi214 = eff_ecut_bi214,
            eff_ecut_tl208 = eff_ecut_tl208,
            eff_fid_bi214 = eff_fid_bi214,
            eff_fid_tl208 = eff_fid_tl208,
            eff_1trk_bi214 = eff_1trk_bi214,
            eff_1trk_tl208 = eff_1trk_tl208,
            events_bi214_raw = events_bi214_raw,
            events_tl208_raw = events_tl208_raw,
            events_bi214_fid = events_bi214_fid,
            events_tl208_fid = events_tl208_fid,
            events_bi214_final = events_bi214_final,
            events_tl208_final = events_tl208_final,
            events_total_raw = events_bi214_raw + events_tl208_raw,
            events_total_fid = events_bi214_fid + events_tl208_fid,
            events_total_final = events_bi214_final + events_tl208_final
        ))
    end

    return results
end

"""
    format_scientific(x::Float64; digits::Int=2)

Format a number in scientific notation with specified digits.
"""
function format_scientific(x::Float64; digits::Int=2)
    if x == 0.0
        return "0.00"
    end
    exp = floor(Int, log10(abs(x)))
    mantissa = x / 10.0^exp
    return @sprintf("%.2fe%+d", mantissa, exp)
end

"""
    zone_to_string(zone::Symbol)

Convert zone symbol to display string.
"""
function zone_to_string(zone::Symbol)
    if zone == :copper_shell
        return "Cu shell"
    elseif zone == :ptfe_barrel
        return "PTFE barrel"
    elseif zone == :ptfe_endcap
        return "PTFE endcap"
    else
        return string(zone)
    end
end

# =============================================================================
# Output Functions
# =============================================================================

"""
    print_activity_fiducial_table(results)

Print the activity table with fiducial efficiencies.
"""
function print_activity_fiducial_table(results)
    println()
    println("=" ^ 120)
    println("ITACA DETECTOR: ACTIVITIES AND FIDUCIAL EFFICIENCIES")
    println("=" ^ 120)
    println()

    # Header
    @printf("%-15s  %-10s  %-12s  %10s  %10s  %12s  %12s\n",
            "Component", "Material", "Zone", "Bi-214(mBq)", "Tl-208(mBq)", "ε_fid(Bi)", "ε_fid(Tl)")
    println("-" ^ 120)

    # Data rows
    total_bi214 = 0.0
    total_tl208 = 0.0

    for r in results
        @printf("%-15s  %-10s  %-12s  %10.4f  %10.4f  %12.2e  %12.2e\n",
                r.name, r.material, zone_to_string(r.zone),
                r.activity_bi214_mBq, r.activity_tl208_mBq,
                r.eff_bi214, r.eff_tl208)
        total_bi214 += r.activity_bi214_mBq
        total_tl208 += r.activity_tl208_mBq
    end

    # Totals
    println("-" ^ 120)
    @printf("%-15s  %-10s  %-12s  %10.4f  %10.4f\n",
            "TOTAL", "---", "---", total_bi214, total_tl208)
    println("=" ^ 120)

    return (total_bi214=total_bi214, total_tl208=total_tl208)
end

"""
    print_events_table_raw(results)

Print raw events per year (before fiducial cut).
"""
function print_events_table_raw(results)
    println()
    println("=" ^ 90)
    println("EVENTS PER YEAR (BEFORE fiducial cut)")
    println("=" ^ 90)
    println()

    @printf("%-15s  %14s  %14s  %14s\n",
            "Component", "Events Bi-214", "Events Tl-208", "Events Total")
    println("-" ^ 90)

    total_bi214 = 0.0
    total_tl208 = 0.0

    for r in results
        @printf("%-15s  %14.2e  %14.2e  %14.2e\n",
                r.name, r.events_bi214_raw, r.events_tl208_raw, r.events_total_raw)
        total_bi214 += r.events_bi214_raw
        total_tl208 += r.events_tl208_raw
    end

    println("-" ^ 90)
    @printf("%-15s  %14.2e  %14.2e  %14.2e\n",
            "TOTAL", total_bi214, total_tl208, total_bi214 + total_tl208)
    println("=" ^ 90)

    return (total_bi214=total_bi214, total_tl208=total_tl208, total=total_bi214+total_tl208)
end

"""
    print_events_table_fiducial(results)

Print events per year after fiducial cut.
"""
function print_events_table_fiducial(results)
    println()
    println("=" ^ 90)
    println("EVENTS PER YEAR (AFTER fiducial cut)")
    println("=" ^ 90)
    println()

    @printf("%-15s  %14s  %14s  %14s\n",
            "Component", "Events Bi-214", "Events Tl-208", "Events Total")
    println("-" ^ 90)

    total_bi214 = 0.0
    total_tl208 = 0.0

    for r in results
        @printf("%-15s  %14.2e  %14.2e  %14.2e\n",
                r.name, r.events_bi214_fid, r.events_tl208_fid, r.events_total_fid)
        total_bi214 += r.events_bi214_fid
        total_tl208 += r.events_tl208_fid
    end

    println("-" ^ 90)
    @printf("%-15s  %14.2e  %14.2e  %14.2e\n",
            "TOTAL", total_bi214, total_tl208, total_bi214 + total_tl208)
    println("=" ^ 90)

    return (total_bi214=total_bi214, total_tl208=total_tl208, total=total_bi214+total_tl208)
end

"""
    print_summary(activity_totals, raw_totals, fid_totals)

Print comprehensive summary.
"""
function print_summary(activity_totals, raw_totals, fid_totals)
    println()
    println("=" ^ 70)
    println("SUMMARY")
    println("=" ^ 70)
    println()

    println("Total Activities:")
    @printf("  Bi-214: %.4f mBq\n", activity_totals.total_bi214)
    @printf("  Tl-208: %.4f mBq\n", activity_totals.total_tl208)
    @printf("  Total:  %.4f mBq\n", activity_totals.total_bi214 + activity_totals.total_tl208)
    println()

    println("Events per Year (BEFORE fiducial cut):")
    @printf("  Bi-214: %.4e events/year\n", raw_totals.total_bi214)
    @printf("  Tl-208: %.4e events/year\n", raw_totals.total_tl208)
    @printf("  Total:  %.4e events/year\n", raw_totals.total)
    println()

    println("Events per Year (AFTER fiducial cut):")
    @printf("  Bi-214: %.4e events/year\n", fid_totals.total_bi214)
    @printf("  Tl-208: %.4e events/year\n", fid_totals.total_tl208)
    @printf("  Total:  %.4e events/year\n", fid_totals.total)
    println()

    # Rejection factors
    rej_bi214 = raw_totals.total_bi214 / fid_totals.total_bi214
    rej_tl208 = raw_totals.total_tl208 / fid_totals.total_tl208
    rej_total = raw_totals.total / fid_totals.total

    println("Fiducial Rejection Factors (raw/fiducial):")
    @printf("  Bi-214: %.1f\n", rej_bi214)
    @printf("  Tl-208: %.1f\n", rej_tl208)
    @printf("  Total:  %.1f\n", rej_total)
    println()

    println("=" ^ 70)
end

"""
    print_summary_with_single_track(activity_totals, raw_totals, fid_totals, final_totals, sigma_label::String)

Print comprehensive summary including single-track efficiency.
"""
function print_summary_with_single_track(activity_totals, raw_totals, fid_totals, final_totals, sigma_label::String)
    println()
    println("=" ^ 70)
    println("SUMMARY ($sigma_label)")
    println("=" ^ 70)
    println()

    println("Total Activities:")
    @printf("  Bi-214: %.4f mBq\n", activity_totals.total_bi214)
    @printf("  Tl-208: %.4f mBq\n", activity_totals.total_tl208)
    @printf("  Total:  %.4f mBq\n", activity_totals.total_bi214 + activity_totals.total_tl208)
    println()

    println("Events per Year (BEFORE fiducial cut):")
    @printf("  Bi-214: %.4e events/year\n", raw_totals.total_bi214)
    @printf("  Tl-208: %.4e events/year\n", raw_totals.total_tl208)
    @printf("  Total:  %.4e events/year\n", raw_totals.total)
    println()

    println("Events per Year (AFTER fiducial cut):")
    @printf("  Bi-214: %.4e events/year\n", fid_totals.total_bi214)
    @printf("  Tl-208: %.4e events/year\n", fid_totals.total_tl208)
    @printf("  Total:  %.4e events/year\n", fid_totals.total)
    println()

    println("Events per Year (AFTER fiducial + single-track cut):")
    @printf("  Bi-214: %.4e events/year\n", final_totals.total_bi214)
    @printf("  Tl-208: %.4e events/year\n", final_totals.total_tl208)
    @printf("  Total:  %.4e events/year\n", final_totals.total)
    println()

    # Rejection factors
    rej_fid_bi214 = raw_totals.total_bi214 / fid_totals.total_bi214
    rej_fid_tl208 = raw_totals.total_tl208 / fid_totals.total_tl208
    rej_fid_total = raw_totals.total / fid_totals.total

    println("Fiducial Rejection Factors (raw/fiducial):")
    @printf("  Bi-214: %.1f\n", rej_fid_bi214)
    @printf("  Tl-208: %.1f\n", rej_fid_tl208)
    @printf("  Total:  %.1f\n", rej_fid_total)
    println()

    # Total rejection factors (raw to final)
    rej_total_bi214 = raw_totals.total_bi214 / final_totals.total_bi214
    rej_total_tl208 = raw_totals.total_tl208 / final_totals.total_tl208
    rej_total_total = raw_totals.total / final_totals.total

    println("Total Rejection Factors (raw/final):")
    @printf("  Bi-214: %.1f\n", rej_total_bi214)
    @printf("  Tl-208: %.1f\n", rej_total_tl208)
    @printf("  Total:  %.1f\n", rej_total_total)
    println()

    println("=" ^ 70)
end

"""
    print_efficiency_table(results, sigma_label::String)

Print table with all efficiencies (fiducial components and single-track).
"""
function print_efficiency_table(results, sigma_label::String)
    println()
    println("=" ^ 180)
    println("ITACA DETECTOR: EFFICIENCIES ($sigma_label)")
    println("=" ^ 180)
    println()

    # Header
    @printf("%-15s  %-12s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n",
            "Component", "Zone",
            "ε_geom(Bi)", "ε_geom(Tl)", "ε_ecut(Bi)", "ε_ecut(Tl)",
            "ε_fid(Bi)", "ε_fid(Tl)", "ε_1trk(Bi)", "ε_1trk(Tl)")
    println("-" ^ 180)

    for r in results
        @printf("%-15s  %-12s  %10.2e  %10.2e  %10.4f  %10.4f  %10.2e  %10.2e  %10.4f  %10.4f\n",
                r.name, zone_to_string(r.zone),
                r.eff_geom_bi214, r.eff_geom_tl208,
                r.eff_ecut_bi214, r.eff_ecut_tl208,
                r.eff_fid_bi214, r.eff_fid_tl208,
                r.eff_1trk_bi214, r.eff_1trk_tl208)
    end

    println("=" ^ 180)
end

"""
    print_events_table_final(results, sigma_label::String)

Print events per year after fiducial AND single-track cut.
"""
function print_events_table_final(results, sigma_label::String)
    println()
    println("=" ^ 90)
    println("EVENTS PER YEAR (AFTER fiducial + single-track cut) - $sigma_label")
    println("=" ^ 90)
    println()

    @printf("%-15s  %14s  %14s  %14s\n",
            "Component", "Events Bi-214", "Events Tl-208", "Events Total")
    println("-" ^ 90)

    total_bi214 = 0.0
    total_tl208 = 0.0

    for r in results
        @printf("%-15s  %14.2e  %14.2e  %14.2e\n",
                r.name, r.events_bi214_final, r.events_tl208_final, r.events_total_final)
        total_bi214 += r.events_bi214_final
        total_tl208 += r.events_tl208_final
    end

    println("-" ^ 90)
    @printf("%-15s  %14.2e  %14.2e  %14.2e\n",
            "TOTAL", total_bi214, total_tl208, total_bi214 + total_tl208)
    println("=" ^ 90)

    return (total_bi214=total_bi214, total_tl208=total_tl208, total=total_bi214+total_tl208)
end

"""
    print_single_track_efficiency_info(st_effs::Dict, csv_path::String, sigma_label::String)

Print information about single-track efficiencies used.
"""
function print_single_track_efficiency_info(st_effs::Dict, csv_path::String, sigma_label::String)
    println()
    println("=" ^ 70)
    println("SINGLE-TRACK EFFICIENCIES ($sigma_label)")
    println("=" ^ 70)
    println()
    println("Source: $csv_path")
    println()
    @printf("%-15s  %15s  %15s\n", "Zone", "ε_1trk (Bi-214)", "ε_1trk (Tl-208)")
    println("-" ^ 70)
    for (zone, st) in sort(collect(st_effs), by=x->string(x[1]))
        @printf("%-15s  %15.6f  %15.6f\n",
                zone_to_string(zone), st.eff_bi214, st.eff_tl208)
    end
    println("=" ^ 70)
end

"""
    print_fiducial_efficiency_info(fid_effs::Dict, csv_path::String)

Print information about fiducial efficiencies used.
Shows all three efficiency components: geometry, energy cut, and total.
"""
function print_fiducial_efficiency_info(fid_effs::Dict, csv_path::String)
    println()
    println("=" ^ 130)
    println("FIDUCIAL EFFICIENCIES (from CSV)")
    println("=" ^ 130)
    println()
    println("Source: $csv_path")
    println()
    @printf("%-15s  %12s  %12s  %12s  %12s  %12s  %12s\n",
            "Zone", "ε_geom(Bi)", "ε_geom(Tl)", "ε_ecut(Bi)", "ε_ecut(Tl)", "ε_fid(Bi)", "ε_fid(Tl)")
    println("-" ^ 130)
    for (zone, fid) in sort(collect(fid_effs), by=x->string(x[1]))
        @printf("%-15s  %12.2e  %12.2e  %12.4f  %12.4f  %12.2e  %12.2e\n",
                zone_to_string(zone),
                fid.eff_geom_bi214, fid.eff_geom_tl208,
                fid.eff_ecut_bi214, fid.eff_ecut_tl208,
                fid.eff_fid_bi214, fid.eff_fid_tl208)
    end
    println("=" ^ 130)
    println()
    println("Legend:")
    println("  ε_geom: Geometry cut efficiency (events in fiducial region / total events)")
    println("  ε_ecut: Energy cut efficiency (events in ROI / events in fiducial region)")
    println("  ε_fid:  Total fiducial efficiency (ε_geom × ε_ecut)")
    println()
    println("Zone → Component mapping:")
    println("  copper_shell → ICS_15cm, FC_rings")
    println("  ptfe_barrel  → BFD")
    println("  ptfe_endcap  → Light_Guides, FAT_GEM, SiPM_boards, Cathode, MARS_*")
end

# =============================================================================
# File Output Functions
# =============================================================================

"""
    write_background_table_csv(results, output_path::String, sigma_label::String)

Write background table to CSV file with all three fiducial efficiency components.
"""
function write_background_table_csv(results, output_path::String, sigma_label::String)
    # Prepare data
    rows = []
    for r in results
        push!(rows, (
            component = r.name,
            material = r.material,
            zone = zone_to_string(r.zone),
            activity_bi214_mBq = r.activity_bi214_mBq,
            activity_tl208_mBq = r.activity_tl208_mBq,
            eff_geom_bi214 = r.eff_geom_bi214,
            eff_geom_tl208 = r.eff_geom_tl208,
            eff_ecut_bi214 = r.eff_ecut_bi214,
            eff_ecut_tl208 = r.eff_ecut_tl208,
            eff_fid_bi214 = r.eff_fid_bi214,
            eff_fid_tl208 = r.eff_fid_tl208,
            eff_1trk_bi214 = r.eff_1trk_bi214,
            eff_1trk_tl208 = r.eff_1trk_tl208,
            events_bi214_raw = r.events_bi214_raw,
            events_tl208_raw = r.events_tl208_raw,
            events_bi214_fid = r.events_bi214_fid,
            events_tl208_fid = r.events_tl208_fid,
            events_bi214_final = r.events_bi214_final,
            events_tl208_final = r.events_tl208_final,
            events_total_final = r.events_total_final
        ))
    end

    # Add totals row
    total_bi214_raw = sum(r.events_bi214_raw for r in results)
    total_tl208_raw = sum(r.events_tl208_raw for r in results)
    total_bi214_fid = sum(r.events_bi214_fid for r in results)
    total_tl208_fid = sum(r.events_tl208_fid for r in results)
    total_bi214_final = sum(r.events_bi214_final for r in results)
    total_tl208_final = sum(r.events_tl208_final for r in results)
    total_final = total_bi214_final + total_tl208_final

    push!(rows, (
        component = "TOTAL",
        material = "-",
        zone = "-",
        activity_bi214_mBq = sum(r.activity_bi214_mBq for r in results),
        activity_tl208_mBq = sum(r.activity_tl208_mBq for r in results),
        eff_geom_bi214 = NaN,
        eff_geom_tl208 = NaN,
        eff_ecut_bi214 = NaN,
        eff_ecut_tl208 = NaN,
        eff_fid_bi214 = NaN,
        eff_fid_tl208 = NaN,
        eff_1trk_bi214 = NaN,
        eff_1trk_tl208 = NaN,
        events_bi214_raw = total_bi214_raw,
        events_tl208_raw = total_tl208_raw,
        events_bi214_fid = total_bi214_fid,
        events_tl208_fid = total_tl208_fid,
        events_bi214_final = total_bi214_final,
        events_tl208_final = total_tl208_final,
        events_total_final = total_final
    ))

    df = DataFrame(rows)
    CSV.write(output_path, df)
    println("  CSV: $output_path")
end

"""
    write_background_table_md(results, output_path::String, sigma_label::String)

Write background table to Markdown file with all three fiducial efficiency components.
"""
function write_background_table_md(results, output_path::String, sigma_label::String)
    open(output_path, "w") do io
        println(io, "# ITACA Background Budget ($sigma_label)")
        println(io)
        println(io, "## Final Background Events per Year")
        println(io)
        println(io, "| Component | Zone | Bi-214 (events/yr) | Tl-208 (events/yr) | Total (events/yr) |")
        println(io, "|-----------|------|-------------------:|-------------------:|------------------:|")

        total_bi214 = 0.0
        total_tl208 = 0.0

        for r in results
            @printf(io, "| %s | %s | %.2e | %.2e | %.2e |\n",
                    r.name, zone_to_string(r.zone),
                    r.events_bi214_final, r.events_tl208_final, r.events_total_final)
            total_bi214 += r.events_bi214_final
            total_tl208 += r.events_tl208_final
        end

        total = total_bi214 + total_tl208
        @printf(io, "| **TOTAL** | - | **%.2e** | **%.2e** | **%.2e** |\n",
                total_bi214, total_tl208, total)

        println(io)
        println(io, "## Efficiencies Used")
        println(io)
        println(io, "| Component | Zone | ε_geom (Bi) | ε_geom (Tl) | ε_ecut (Bi) | ε_ecut (Tl) | ε_fid (Bi) | ε_fid (Tl) | ε_1trk (Bi) | ε_1trk (Tl) |")
        println(io, "|-----------|------|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|")

        for r in results
            @printf(io, "| %s | %s | %.2e | %.2e | %.4f | %.4f | %.2e | %.2e | %.4f | %.4f |\n",
                    r.name, zone_to_string(r.zone),
                    r.eff_geom_bi214, r.eff_geom_tl208,
                    r.eff_ecut_bi214, r.eff_ecut_tl208,
                    r.eff_fid_bi214, r.eff_fid_tl208,
                    r.eff_1trk_bi214, r.eff_1trk_tl208)
        end

        println(io)
        println(io, "**Legend:**")
        println(io, "- ε_geom: Geometry cut efficiency (events in fiducial region / total events)")
        println(io, "- ε_ecut: Energy cut efficiency (events in ROI / events in fiducial region)")
        println(io, "- ε_fid: Total fiducial efficiency (ε_geom × ε_ecut)")
        println(io, "- ε_1trk: Single-track efficiency")

        println(io)
        println(io, "## Activities")
        println(io)
        println(io, "| Component | Material | Mass (kg) | Bi-214 (mBq) | Tl-208 (mBq) |")
        println(io, "|-----------|----------|----------:|-------------:|-------------:|")

        total_bi214_act = 0.0
        total_tl208_act = 0.0

        for r in results
            @printf(io, "| %s | %s | %.2f | %.4f | %.4f |\n",
                    r.name, r.material, r.mass_kg,
                    r.activity_bi214_mBq, r.activity_tl208_mBq)
            total_bi214_act += r.activity_bi214_mBq
            total_tl208_act += r.activity_tl208_mBq
        end

        @printf(io, "| **TOTAL** | - | - | **%.4f** | **%.4f** |\n",
                total_bi214_act, total_tl208_act)
    end

    println("  Markdown: $output_path")
end

"""
    sigma_label_to_latex(sigma_label::String)

Convert sigma label to LaTeX-safe format.
"""
function sigma_label_to_latex(sigma_label::String)
    # Replace σ with \\sigma in math mode
    latex_label = replace(sigma_label, "σ" => "\$\\sigma\$")
    return latex_label
end

"""
    sigma_label_to_label_id(sigma_label::String)

Convert sigma label to a clean label ID for LaTeX references.
"""
function sigma_label_to_label_id(sigma_label::String)
    # Create clean label: "σ=1mm, pitch=2mm" -> "sigma1mm_pitch2mm"
    label_id = replace(sigma_label, "σ=" => "sigma")
    label_id = replace(label_id, ", " => "_")
    label_id = replace(label_id, "=" => "")
    return label_id
end

"""
    write_background_table_tex(results, output_path::String, sigma_label::String; signal_eff::Union{SignalEfficiency, Nothing}=nothing)

Write background table to LaTeX file with all three fiducial efficiency components.
Optionally includes 0νββ signal efficiency table if signal_eff is provided.
"""
function write_background_table_tex(results, output_path::String, sigma_label::String; signal_eff::Union{SignalEfficiency, Nothing}=nothing)
    latex_label = sigma_label_to_latex(sigma_label)
    label_id = sigma_label_to_label_id(sigma_label)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")

    open(output_path, "w") do io
        # LaTeX preamble
        println(io, "% Generated: $timestamp")
        println(io, "\\documentclass{article}")
        println(io, "\\usepackage{booktabs}")
        println(io, "\\usepackage{siunitx}")
        println(io)
        println(io, "% ITACA Background Budget ($sigma_label)")
        println(io, "% Generated by itaca_bkgnd.jl")
        println(io)
        println(io, "\\begin{document}")
        println(io)

        # Final background table
        println(io, "% Final Background Events per Year")
        println(io, "\\begin{table}[htbp]")
        println(io, "\\centering")
        println(io, "\\caption{ITACA Background Budget ($latex_label)}")
        println(io, "\\label{tab:background_$label_id}")
        println(io, "\\begin{tabular}{lrrr}")
        println(io, "\\toprule")
        println(io, "Component & Bi-214 (ev/yr) & Tl-208 (ev/yr) & Total (ev/yr) \\\\")
        println(io, "\\midrule")

        total_bi214 = 0.0
        total_tl208 = 0.0

        for r in results
            @printf(io, "%s & %.2e & %.2e & %.2e \\\\\n",
                    format_component_name_latex(r.name),
                    r.events_bi214_final, r.events_tl208_final, r.events_total_final)
            total_bi214 += r.events_bi214_final
            total_tl208 += r.events_tl208_final
        end

        total = total_bi214 + total_tl208
        println(io, "\\midrule")
        @printf(io, "\\textbf{TOTAL} & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.2e} \\\\\n",
                total_bi214, total_tl208, total)

        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")

        println(io)

        # Efficiencies table for Bi-214
        println(io, "% Efficiencies Used - Bi-214")
        println(io, "\\begin{table}[htbp]")
        println(io, "\\centering")
        println(io, "\\caption{ITACA Efficiencies for Bi-214 ($latex_label)}")
        println(io, "\\label{tab:efficiency_bi214_$label_id}")
        println(io, "\\begin{tabular}{lrrrr}")
        println(io, "\\toprule")
        println(io, "Component & \$\\varepsilon_{\\mathrm{geom}}\$ & \$\\varepsilon_{\\mathrm{ecut}}\$ & \$\\varepsilon_{\\mathrm{fid}}\$ & \$\\varepsilon_{\\mathrm{1trk}}\$ \\\\")
        println(io, "\\midrule")

        for r in results
            @printf(io, "%s & %.2e & %.4f & %.2e & %.4f \\\\\n",
                    format_component_name_latex(r.name),
                    r.eff_geom_bi214, r.eff_ecut_bi214, r.eff_fid_bi214, r.eff_1trk_bi214)
        end

        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")

        println(io)

        # Efficiencies table for Tl-208
        println(io, "% Efficiencies Used - Tl-208")
        println(io, "\\begin{table}[htbp]")
        println(io, "\\centering")
        println(io, "\\caption{ITACA Efficiencies for Tl-208 ($latex_label)}")
        println(io, "\\label{tab:efficiency_tl208_$label_id}")
        println(io, "\\begin{tabular}{lrrrr}")
        println(io, "\\toprule")
        println(io, "Component & \$\\varepsilon_{\\mathrm{geom}}\$ & \$\\varepsilon_{\\mathrm{ecut}}\$ & \$\\varepsilon_{\\mathrm{fid}}\$ & \$\\varepsilon_{\\mathrm{1trk}}\$ \\\\")
        println(io, "\\midrule")

        for r in results
            @printf(io, "%s & %.2e & %.4f & %.2e & %.4f \\\\\n",
                    format_component_name_latex(r.name),
                    r.eff_geom_tl208, r.eff_ecut_tl208, r.eff_fid_tl208, r.eff_1trk_tl208)
        end

        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")

        println(io)

        # Signal efficiency table (0νββ)
        if signal_eff !== nothing
            println(io, "% Signal Efficiency - 0νββ")
            println(io, "\\begin{table}[htbp]")
            println(io, "\\centering")
            println(io, "\\caption{ITACA Signal Efficiency for \$0\\nu\\beta\\beta\$ ($latex_label)}")
            println(io, "\\label{tab:efficiency_0nubb_$label_id}")
            println(io, "\\begin{tabular}{lrrrr}")
            println(io, "\\toprule")
            println(io, " & \$\\varepsilon_{\\mathrm{geom}}\$ & \$\\varepsilon_{\\mathrm{ecut}}\$ & \$\\varepsilon_{\\mathrm{fid}}\$ & \$\\varepsilon_{\\mathrm{1trk}}\$ \\\\")
            println(io, "\\midrule")
            @printf(io, "\$0\\nu\\beta\\beta\$ & %.4f & %.4f & %.4f & %.4f \\\\\n",
                    signal_eff.eff_geom, signal_eff.eff_ecut, signal_eff.eff_fid, signal_eff.eff_1trk)
            println(io, "\\midrule")
            @printf(io, "\\textbf{Total signal efficiency} & \\multicolumn{4}{r}{\\textbf{%.4f}} \\\\\n", signal_eff.eff_total)
            println(io, "\\bottomrule")
            println(io, "\\end{tabular}")
            println(io, "\\end{table}")
            println(io)
        end

        println(io, "% Legend:")
        println(io, "% ε_geom: Geometry cut efficiency (events in fiducial region / total events)")
        println(io, "% ε_ecut: Energy cut efficiency (events in ROI / events in fiducial region)")
        println(io, "% ε_fid: Total fiducial efficiency (ε_geom × ε_ecut)")
        println(io, "% ε_1trk: Single-track efficiency")

        println(io)
        println(io, "\\end{document}")
    end

    println("  LaTeX: $output_path")
end

"""
    escape_latex(s::String)

Escape special LaTeX characters.
"""
function escape_latex(s::String)
    s = replace(s, "_" => "\\_")
    s = replace(s, "&" => "\\&")
    s = replace(s, "%" => "\\%")
    s = replace(s, "#" => "\\#")
    return s
end

"""
    format_component_name_latex(name::String)

Format component names for LaTeX tables.
Converts underscore-separated names to readable format:
- ICS_15cm → ICS (15 cm)
- FC_rings → FC rings
- Light_Guides → Light Guides
- FAT_GEM → FAT GEM
- SiPM_boards → SiPM boards
- MARS_Ti → MARS Ti
- MARS_Kapton → MARS Kapton
- MARS_HDPE → MARS HDPE
"""
function format_component_name_latex(name::String)
    # Special case for ICS_15cm
    if name == "ICS_15cm"
        return "ICS (15 cm)"
    end
    # Replace underscores with spaces for other names
    return replace(name, "_" => " ")
end

"""
    write_comparison_table_tex(results_1mm, results_10mm, output_path::String)

Write comparison table (σ=1mm vs σ=10mm) to LaTeX file.
"""
function write_comparison_table_tex(results_1mm, results_10mm, output_path::String)
    # Compute ratios
    ratio_bi214 = results_10mm.final_totals.total_bi214 / results_1mm.final_totals.total_bi214
    ratio_tl208 = results_10mm.final_totals.total_tl208 / results_1mm.final_totals.total_tl208
    ratio_total = results_10mm.final_totals.total / results_1mm.final_totals.total
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")

    open(output_path, "w") do io
        println(io, "% Generated: $timestamp")
        println(io, "\\documentclass{article}")
        println(io, "\\usepackage{booktabs}")
        println(io, "\\usepackage{siunitx}")
        println(io)
        println(io, "% ITACA Background Comparison: σ=1mm vs σ=10mm")
        println(io, "% Generated by itaca_bkgnd.jl")
        println(io)
        println(io, "\\begin{document}")
        println(io)

        println(io, "\\begin{table}[htbp]")
        println(io, "\\centering")
        println(io, "\\caption{ITACA Background Comparison: \$\\sigma\$=1mm vs \$\\sigma\$=10mm}")
        println(io, "\\label{tab:background_comparison}")
        println(io, "\\begin{tabular}{lrrr}")
        println(io, "\\toprule")
        println(io, "Isotope & \$\\sigma\$=1mm (ev/yr) & \$\\sigma\$=10mm (ev/yr) & Ratio \\\\")
        println(io, "\\midrule")

        @printf(io, "Bi-214 & %.4e & %.4e & %.2f \\\\\n",
                results_1mm.final_totals.total_bi214, results_10mm.final_totals.total_bi214, ratio_bi214)
        @printf(io, "Tl-208 & %.4e & %.4e & %.2f \\\\\n",
                results_1mm.final_totals.total_tl208, results_10mm.final_totals.total_tl208, ratio_tl208)
        println(io, "\\midrule")
        @printf(io, "\\textbf{Total} & \\textbf{%.4e} & \\textbf{%.4e} & \\textbf{%.2f} \\\\\n",
                results_1mm.final_totals.total, results_10mm.final_totals.total, ratio_total)

        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")
        println(io)
        println(io, "\\end{document}")
    end

    println("  LaTeX comparison: $output_path")
end

"""
    write_comparison_table_md(results_1mm, results_10mm, output_path::String)

Write comparison table (σ=1mm vs σ=10mm) to Markdown file.
"""
function write_comparison_table_md(results_1mm, results_10mm, output_path::String)
    # Compute ratios
    ratio_bi214 = results_10mm.final_totals.total_bi214 / results_1mm.final_totals.total_bi214
    ratio_tl208 = results_10mm.final_totals.total_tl208 / results_1mm.final_totals.total_tl208
    ratio_total = results_10mm.final_totals.total / results_1mm.final_totals.total

    open(output_path, "w") do io
        println(io, "# ITACA Background Comparison: σ=1mm vs σ=10mm")
        println(io)
        println(io, "Final background (after fiducial + single-track cuts):")
        println(io)
        println(io, "| Isotope | σ=1mm (ev/yr) | σ=10mm (ev/yr) | Ratio |")
        println(io, "|---------|-------------:|---------------:|------:|")
        @printf(io, "| Bi-214 | %.4e | %.4e | %.2f |\n",
                results_1mm.final_totals.total_bi214, results_10mm.final_totals.total_bi214, ratio_bi214)
        @printf(io, "| Tl-208 | %.4e | %.4e | %.2f |\n",
                results_1mm.final_totals.total_tl208, results_10mm.final_totals.total_tl208, ratio_tl208)
        @printf(io, "| **Total** | **%.4e** | **%.4e** | **%.2f** |\n",
                results_1mm.final_totals.total, results_10mm.final_totals.total, ratio_total)
    end

    println("  Markdown comparison: $output_path")
end

"""
    write_combined_tables_tex(results_ion, results_ele, signal_eff_ion::SignalEfficiency, signal_eff_ele::SignalEfficiency, output_path::String)

Write combined efficiency and background tables for ion (σ=1mm) and electron (σ=10mm) tracks.
Since fiducial efficiencies are the same for both, we show them once with two columns for ε_1trk.

Tables generated:
1. Combined Efficiencies for Bi-214 (ε_geom, ε_ecut, ε_fid, ε_1trk^ion, ε_1trk^ele)
2. Combined Efficiencies for Tl-208 (same structure)
3. Combined Efficiencies for 0νββ signal
4. Bi-214 Background Events (ion vs electron columns)
5. Tl-208 Background Events (ion vs electron columns)
6. Total Background Summary (ion vs electron with ratio)
"""
function write_combined_tables_tex(results_ion, results_ele, signal_eff_ion::SignalEfficiency, signal_eff_ele::SignalEfficiency, output_path::String)
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    open(output_path, "w") do io
        # LaTeX preamble
        println(io, "% Generated: $timestamp")
        println(io, "\\documentclass{article}")
        println(io, "\\usepackage{booktabs}")
        println(io, "\\usepackage{siunitx}")
        println(io)
        println(io, "% ITACA Combined Efficiency Tables (Ion and Electron Tracks)")
        println(io, "% Ion track: σ=1mm (small diffusion)")
        println(io, "% Electron track: σ=10mm (large diffusion)")
        println(io, "% Generated by itaca_bkgnd.jl")
        println(io)
        println(io, "\\begin{document}")
        println(io)

        # =================================================================
        # Table 1: Combined Efficiencies for Bi-214
        # =================================================================
        println(io, "% Combined Efficiencies - Bi-214")
        println(io, "\\begin{table}[htbp]")
        println(io, "\\centering")
        println(io, "\\caption{ITACA Efficiencies for Bi-214 (Ion and Electron Tracks)}")
        println(io, "\\label{tab:efficiency_bi214_combined}")
        println(io, "\\begin{tabular}{lrrrrr}")
        println(io, "\\toprule")
        println(io, "Component & \$\\varepsilon_{\\mathrm{geom}}\$ & \$\\varepsilon_{\\mathrm{ecut}}\$ & \$\\varepsilon_{\\mathrm{fid}}\$ & \$\\varepsilon_{\\mathrm{1trk}}^{\\mathrm{ion}}\$ & \$\\varepsilon_{\\mathrm{1trk}}^{\\mathrm{ele}}\$ \\\\")
        println(io, "\\midrule")

        for (r_ion, r_ele) in zip(results_ion, results_ele)
            @printf(io, "%s & %.2e & %.4f & %.2e & %.4f & %.4f \\\\\n",
                    format_component_name_latex(r_ion.name),
                    r_ion.eff_geom_bi214, r_ion.eff_ecut_bi214, r_ion.eff_fid_bi214,
                    r_ion.eff_1trk_bi214, r_ele.eff_1trk_bi214)
        end

        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")
        println(io)

        # =================================================================
        # Table 2: Combined Efficiencies for Tl-208
        # =================================================================
        println(io, "% Combined Efficiencies - Tl-208")
        println(io, "\\begin{table}[htbp]")
        println(io, "\\centering")
        println(io, "\\caption{ITACA Efficiencies for Tl-208 (Ion and Electron Tracks)}")
        println(io, "\\label{tab:efficiency_tl208_combined}")
        println(io, "\\begin{tabular}{lrrrrr}")
        println(io, "\\toprule")
        println(io, "Component & \$\\varepsilon_{\\mathrm{geom}}\$ & \$\\varepsilon_{\\mathrm{ecut}}\$ & \$\\varepsilon_{\\mathrm{fid}}\$ & \$\\varepsilon_{\\mathrm{1trk}}^{\\mathrm{ion}}\$ & \$\\varepsilon_{\\mathrm{1trk}}^{\\mathrm{ele}}\$ \\\\")
        println(io, "\\midrule")

        for (r_ion, r_ele) in zip(results_ion, results_ele)
            @printf(io, "%s & %.2e & %.4f & %.2e & %.4f & %.4f \\\\\n",
                    format_component_name_latex(r_ion.name),
                    r_ion.eff_geom_tl208, r_ion.eff_ecut_tl208, r_ion.eff_fid_tl208,
                    r_ion.eff_1trk_tl208, r_ele.eff_1trk_tl208)
        end

        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")
        println(io)

        # =================================================================
        # Table 3: Combined Efficiencies for 0νββ Signal
        # =================================================================
        println(io, "% Combined Efficiencies - 0νββ Signal")
        println(io, "\\begin{table}[htbp]")
        println(io, "\\centering")
        println(io, "\\caption{ITACA Signal Efficiency for \$0\\nu\\beta\\beta\$ (Ion and Electron Tracks)}")
        println(io, "\\label{tab:efficiency_0nubb_combined}")
        println(io, "\\begin{tabular}{lrrrrr}")
        println(io, "\\toprule")
        println(io, " & \$\\varepsilon_{\\mathrm{geom}}\$ & \$\\varepsilon_{\\mathrm{ecut}}\$ & \$\\varepsilon_{\\mathrm{fid}}\$ & \$\\varepsilon_{\\mathrm{1trk}}^{\\mathrm{ion}}\$ & \$\\varepsilon_{\\mathrm{1trk}}^{\\mathrm{ele}}\$ \\\\")
        println(io, "\\midrule")
        @printf(io, "\$0\\nu\\beta\\beta\$ & %.4f & %.4f & %.4f & %.4f & %.4f \\\\\n",
                signal_eff_ion.eff_geom, signal_eff_ion.eff_ecut, signal_eff_ion.eff_fid,
                signal_eff_ion.eff_1trk, signal_eff_ele.eff_1trk)
        println(io, "\\midrule")
        @printf(io, "\\textbf{Total signal efficiency} & & & & \\textbf{%.4f} & \\textbf{%.4f} \\\\\n",
                signal_eff_ion.eff_total, signal_eff_ele.eff_total)
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")
        println(io)

        # =================================================================
        # Table 4: Bi-214 Background Events
        # =================================================================
        println(io, "% Bi-214 Background Events")
        println(io, "\\begin{table}[htbp]")
        println(io, "\\centering")
        println(io, "\\caption{ITACA Bi-214 Background (Ion and Electron Tracks)}")
        println(io, "\\label{tab:background_bi214_combined}")
        println(io, "\\begin{tabular}{lrr}")
        println(io, "\\toprule")
        println(io, "Component & Ion (ev/yr) & Electron (ev/yr) \\\\")
        println(io, "\\midrule")

        total_ion_bi214 = 0.0
        total_ele_bi214 = 0.0
        for (r_ion, r_ele) in zip(results_ion, results_ele)
            @printf(io, "%s & %.2e & %.2e \\\\\n",
                    format_component_name_latex(r_ion.name),
                    r_ion.events_bi214_final, r_ele.events_bi214_final)
            total_ion_bi214 += r_ion.events_bi214_final
            total_ele_bi214 += r_ele.events_bi214_final
        end

        println(io, "\\midrule")
        @printf(io, "\\textbf{TOTAL} & \\textbf{%.2e} & \\textbf{%.2e} \\\\\n",
                total_ion_bi214, total_ele_bi214)
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")
        println(io)

        # =================================================================
        # Table 5: Tl-208 Background Events
        # =================================================================
        println(io, "% Tl-208 Background Events")
        println(io, "\\begin{table}[htbp]")
        println(io, "\\centering")
        println(io, "\\caption{ITACA Tl-208 Background (Ion and Electron Tracks)}")
        println(io, "\\label{tab:background_tl208_combined}")
        println(io, "\\begin{tabular}{lrr}")
        println(io, "\\toprule")
        println(io, "Component & Ion (ev/yr) & Electron (ev/yr) \\\\")
        println(io, "\\midrule")

        total_ion_tl208 = 0.0
        total_ele_tl208 = 0.0
        for (r_ion, r_ele) in zip(results_ion, results_ele)
            @printf(io, "%s & %.2e & %.2e \\\\\n",
                    format_component_name_latex(r_ion.name),
                    r_ion.events_tl208_final, r_ele.events_tl208_final)
            total_ion_tl208 += r_ion.events_tl208_final
            total_ele_tl208 += r_ele.events_tl208_final
        end

        println(io, "\\midrule")
        @printf(io, "\\textbf{TOTAL} & \\textbf{%.2e} & \\textbf{%.2e} \\\\\n",
                total_ion_tl208, total_ele_tl208)
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")
        println(io)

        # =================================================================
        # Table 6: Total Background Summary
        # =================================================================
        println(io, "% Total Background Summary")
        println(io, "\\begin{table}[htbp]")
        println(io, "\\centering")
        println(io, "\\caption{ITACA Total Background Summary (Ion and Electron Tracks)}")
        println(io, "\\label{tab:background_total_combined}")
        println(io, "\\begin{tabular}{lrrr}")
        println(io, "\\toprule")
        println(io, "Isotope & Ion (ev/yr) & Electron (ev/yr) & Ratio \\\\")
        println(io, "\\midrule")

        ratio_bi214 = total_ele_bi214 / total_ion_bi214
        ratio_tl208 = total_ele_tl208 / total_ion_tl208
        total_ion = total_ion_bi214 + total_ion_tl208
        total_ele = total_ele_bi214 + total_ele_tl208
        ratio_total = total_ele / total_ion

        @printf(io, "Bi-214 & %.2e & %.2e & %.2f \\\\\n", total_ion_bi214, total_ele_bi214, ratio_bi214)
        @printf(io, "Tl-208 & %.2e & %.2e & %.2f \\\\\n", total_ion_tl208, total_ele_tl208, ratio_tl208)
        println(io, "\\midrule")
        @printf(io, "\\textbf{Total} & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.2f} \\\\\n",
                total_ion, total_ele, ratio_total)
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
        println(io, "\\end{table}")
        println(io)

        # Legend
        println(io, "% Legend:")
        println(io, "% ε_geom: Geometry cut efficiency (events in fiducial region / total events)")
        println(io, "% ε_ecut: Energy cut efficiency (events in ROI / events in fiducial region)")
        println(io, "% ε_fid: Total fiducial efficiency (ε_geom × ε_ecut)")
        println(io, "% ε_1trk^ion: Single-track efficiency for ion track (σ=1mm)")
        println(io, "% ε_1trk^ele: Single-track efficiency for electron track (σ=10mm)")
        println(io)
        println(io, "\\end{document}")
    end

    println("  LaTeX combined: $output_path")
end

"""
    write_background_tables(results, output_dir::String, prefix::String, sigma_suffix::String, sigma_label::String; signal_eff::Union{SignalEfficiency, Nothing}=nothing)

Write background tables to all formats (CSV, Markdown, LaTeX).
Optionally includes 0νββ signal efficiency in LaTeX output if signal_eff is provided.
"""
function write_background_tables(results, output_dir::String, prefix::String, sigma_suffix::String, sigma_label::String; signal_eff::Union{SignalEfficiency, Nothing}=nothing)
    # Ensure output directory exists
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    println()
    println("Writing background tables for $sigma_label:")

    # Write CSV
    csv_path = joinpath(output_dir, "$(prefix)_$(sigma_suffix).csv")
    write_background_table_csv(results, csv_path, sigma_label)

    # Write Markdown
    md_path = joinpath(output_dir, "$(prefix)_$(sigma_suffix).md")
    write_background_table_md(results, md_path, sigma_label)

    # Write LaTeX (with signal efficiency if provided)
    tex_path = joinpath(output_dir, "$(prefix)_$(sigma_suffix).tex")
    write_background_table_tex(results, tex_path, sigma_label; signal_eff=signal_eff)
end

# =============================================================================
# Main
# =============================================================================

"""
    run_analysis_for_sigma(fid_effs::Dict, st_effs::Dict, sigma::Symbol, sigma_label::String,
                           fiducial_csv::String, efficiency_csv::String;
                           output_dir::String="", output_prefix::String="",
                           signal_info::String=DEFAULT_SIGNAL_INFO, signal_geom_eff::Float64=DEFAULT_SIGNAL_GEOM_EFF)

Run the full background analysis for a specific sigma configuration.
If output_dir is provided, writes tables to files in all formats.

Arguments:
- signal_info: Path to 0νββ info file containing energy cut efficiency
- signal_geom_eff: 0νββ geometry efficiency (default 0.95)
"""
function run_analysis_for_sigma(fid_effs::Dict, st_effs::Dict, sigma::Symbol, sigma_label::String,
                                fiducial_csv::String, efficiency_csv::String;
                                output_dir::String="", output_prefix::String="",
                                signal_info::String=DEFAULT_SIGNAL_INFO, signal_geom_eff::Float64=DEFAULT_SIGNAL_GEOM_EFF)
    println()
    println("#" ^ 80)
    println("#" * " " ^ 30 * sigma_label * " " ^ (48 - length(sigma_label)) * "#")
    println("#" ^ 80)

    # Print efficiency info
    print_fiducial_efficiency_info(fid_effs, fiducial_csv)
    print_single_track_efficiency_info(st_effs, efficiency_csv, sigma_label)

    # Compute events with both efficiencies
    results = compute_events_with_single_track(ITACA_COMPONENTS, fid_effs, st_effs)

    # Print efficiency table
    print_efficiency_table(results, sigma_label)

    # Print raw events table
    raw_totals = print_events_table_raw(results)

    # Print fiducial events table
    fid_totals = print_events_table_fiducial(results)

    # Print final events table (after single-track cut)
    final_totals = print_events_table_final(results, sigma_label)

    # Compute activity totals
    total_bi214 = sum(r.activity_bi214_mBq for r in results)
    total_tl208 = sum(r.activity_tl208_mBq for r in results)
    activity_totals = (total_bi214=total_bi214, total_tl208=total_tl208)

    # Print summary
    print_summary_with_single_track(activity_totals, raw_totals, fid_totals, final_totals, sigma_label)

    # Load signal efficiencies for 0νββ
    signal_eff = load_signal_efficiencies(efficiency_csv, signal_info, sigma; geom_eff=signal_geom_eff)

    # Write tables to files if output directory is specified
    if !isempty(output_dir)
        sigma_suffix = sigma == :sigma_1mm ? "sigma1mm" : "sigma10mm"
        write_background_tables(results, output_dir, output_prefix, sigma_suffix, sigma_label; signal_eff=signal_eff)
    end

    return (results=results, signal_eff=signal_eff, raw_totals=raw_totals, fid_totals=fid_totals, final_totals=final_totals)
end

function main()
    # Default values
    fiducial_csv = DEFAULT_FIDUCIAL_CSV
    efficiency_csv = DEFAULT_EFFICIENCY_TABLE
    output_dir = DEFAULT_OUTPUT_DIR
    output_prefix = DEFAULT_OUTPUT_PREFIX

    # Parse arguments
    for arg in ARGS
        if arg == "--help" || arg == "-h"
            println("""
Usage: julia itaca_bkgnd.jl [options]

Computes the ITACA detector background budget from radioactivity data.
Generates separate tables for σ=1mm and σ=10mm diffusion scenarios.

Options:
    --fiducial=X        Path to fiducial efficiency CSV (default: datainfo/fiducial_efficiency.csv)
    --efficiency-table=X  Path to single-track efficiency CSV (default: tables/efficiency_table.csv)
    --output-dir=X      Output directory for tables (default: tables/)
    --output-prefix=X   Prefix for output files (default: background_tables)
    --help, -h          Show this help message

Output:
    - Efficiency tables (fiducial and single-track) for each sigma
    - Events per year: raw, after fiducial, after fiducial+single-track
    - Summary with rejection factors for each sigma scenario
    - Tables in CSV, Markdown, and LaTeX formats
""")
            return
        elseif startswith(arg, "--fiducial=")
            fiducial_csv = String(split(arg, '=')[2])
        elseif startswith(arg, "--efficiency-table=")
            efficiency_csv = String(split(arg, '=')[2])
        elseif startswith(arg, "--output-dir=")
            output_dir = String(split(arg, '=')[2])
        elseif startswith(arg, "--output-prefix=")
            output_prefix = String(split(arg, '=')[2])
        else
            println("Warning: Unknown argument: $arg")
        end
    end

    println()
    println("*" ^ 80)
    println("*" * " " ^ 24 * "ITACA BACKGROUND CALCULATOR" * " " ^ 27 * "*")
    println("*" ^ 80)
    println()
    println("Fiducial efficiency CSV:     $fiducial_csv")
    println("Single-track efficiency CSV: $efficiency_csv")
    println("Output directory:            $output_dir")
    println("Output prefix:               $output_prefix")

    # Load fiducial efficiencies (same for both sigma cases)
    fid_effs = load_fiducial_efficiencies(fiducial_csv)

    # Load single-track efficiencies for both sigma cases
    st_effs_1mm = load_single_track_efficiencies(efficiency_csv, :sigma_1mm)
    st_effs_10mm = load_single_track_efficiencies(efficiency_csv, :sigma_10mm)

    # Run analysis for σ=1mm
    results_1mm = run_analysis_for_sigma(fid_effs, st_effs_1mm, :sigma_1mm, "σ=1mm, pitch=2mm",
                                          fiducial_csv, efficiency_csv;
                                          output_dir=output_dir, output_prefix=output_prefix)

    # Run analysis for σ=10mm
    results_10mm = run_analysis_for_sigma(fid_effs, st_effs_10mm, :sigma_10mm, "σ=10mm, pitch=20mm",
                                           fiducial_csv, efficiency_csv;
                                           output_dir=output_dir, output_prefix=output_prefix)

    # Print comparison
    println()
    println("=" ^ 80)
    println("COMPARISON: σ=1mm vs σ=10mm")
    println("=" ^ 80)
    println()
    println("Final background (after fiducial + single-track):")
    println()
    @printf("  %-20s  %15s  %15s  %10s\n", "Isotope", "σ=1mm", "σ=10mm", "Ratio")
    println("-" ^ 80)
    ratio_bi214 = results_10mm.final_totals.total_bi214 / results_1mm.final_totals.total_bi214
    ratio_tl208 = results_10mm.final_totals.total_tl208 / results_1mm.final_totals.total_tl208
    ratio_total = results_10mm.final_totals.total / results_1mm.final_totals.total
    @printf("  %-20s  %15.4e  %15.4e  %10.2f\n", "Bi-214",
            results_1mm.final_totals.total_bi214, results_10mm.final_totals.total_bi214, ratio_bi214)
    @printf("  %-20s  %15.4e  %15.4e  %10.2f\n", "Tl-208",
            results_1mm.final_totals.total_tl208, results_10mm.final_totals.total_tl208, ratio_tl208)
    @printf("  %-20s  %15.4e  %15.4e  %10.2f\n", "Total",
            results_1mm.final_totals.total, results_10mm.final_totals.total, ratio_total)
    println("=" ^ 80)

    # Write comparison tables to files
    println()
    println("Writing comparison tables:")
    comparison_tex = joinpath(output_dir, "$(output_prefix)_comparison.tex")
    write_comparison_table_tex(results_1mm, results_10mm, comparison_tex)
    comparison_md = joinpath(output_dir, "$(output_prefix)_comparison.md")
    write_comparison_table_md(results_1mm, results_10mm, comparison_md)

    # Write combined ion/electron efficiency tables
    println()
    println("Writing combined ion/electron tables:")
    combined_tex = joinpath(output_dir, "efficiencies_combined.tex")
    write_combined_tables_tex(results_1mm.results, results_10mm.results,
                               results_1mm.signal_eff, results_10mm.signal_eff,
                               combined_tex)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
