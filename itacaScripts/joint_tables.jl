#!/usr/bin/env julia

"""
Generate a joint LaTeX document with consolidated efficiency tables.

Reads individual .tex tables from the tables/ directory and creates a single
LaTeX article document with:
1. Table of fiducial efficiencies (ε_fid)
2. Table of single-track efficiencies (ε_1trk) with "Ion track" and "Electron track" columns

Usage:
    julia joint_tables.jl [--output=<file>]

Options:
    --output=X    Output file (default: tables/joint_tables.tex)
"""

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using Printf

# =============================================================================
# Configuration
# =============================================================================

const TABLES_DIR = joinpath(dirname(@__FILE__), "tables")

# Files to read (in order)
const TABLE_FILES = [
    (isotope="0nubb", zone=nothing, file="efficiency_0nubb.tex"),
    (isotope="bi214", zone="copper_shell", file="efficiency_bi214_copper_shell.tex"),
    (isotope="bi214", zone="copper_endcaps", file="efficiency_bi214_copper_endcaps.tex"),
    (isotope="bi214", zone="weighted_mean", file="efficiency_bi214_weighted_mean.tex"),
    (isotope="tl208", zone="copper_shell", file="efficiency_tl208_copper_shell.tex"),
    (isotope="tl208", zone="copper_endcaps", file="efficiency_tl208_copper_endcaps.tex"),
    (isotope="tl208", zone="weighted_mean", file="efficiency_tl208_weighted_mean.tex"),
]

# Display names for rows
const ROW_NAMES = Dict(
    ("0nubb", nothing) => "\$0\\nu\\beta\\beta\$",
    ("bi214", "copper_shell") => "\$^{214}\$Bi copper shell",
    ("bi214", "copper_endcaps") => "\$^{214}\$Bi copper endcaps",
    ("bi214", "weighted_mean") => "\$^{214}\$Bi average",
    ("tl208", "copper_shell") => "\$^{208}\$Tl copper shell",
    ("tl208", "copper_endcaps") => "\$^{208}\$Tl copper endcaps",
    ("tl208", "weighted_mean") => "\$^{208}\$Tl average",
)

# =============================================================================
# Parsing Functions
# =============================================================================

"""
    parse_scientific_notation(s::String)

Parse LaTeX scientific notation like "\$2.93 \\times 10^{-4}\$" to Float64.
"""
function parse_scientific_notation(s::AbstractString)
    # Extract mantissa and exponent from pattern like "$2.93 \times 10^{-4}$"
    # Pattern: $X.XX \times 10^{Y}$
    # Use raw string to avoid escaping issues
    pattern = r"\$([0-9.]+)\s*\\times\s*10\^\{([+-]?[0-9]+)\}\$"
    m = match(pattern, s)
    if m === nothing
        error("Could not parse scientific notation: $s")
    end
    mantissa = parse(Float64, m.captures[1])
    exponent = parse(Int, m.captures[2])
    return mantissa * 10.0^exponent
end

"""
    parse_tex_table(filepath::String)

Parse a .tex table file and extract efficiency values.
Returns (ε_fid_2mm, ε_fid_20mm, ε_1trk_2mm, ε_1trk_20mm).
"""
function parse_tex_table(filepath::String)
    if !isfile(filepath)
        error("File not found: $filepath")
    end

    content = read(filepath, String)

    # Find ε_fid row: $\epsilon_{\text{fid}}$ & value1 & value2 \\
    # Capture groups: value between & markers, ending at the line break \\
    # Use .+? for non-greedy matching to capture the full LaTeX expression
    fid_match = match(r"\$\\epsilon_\{\\text\{fid\}\}\$\s*&\s*(\$[^$]+\$)\s*&\s*(\$[^$]+\$)", content)
    if fid_match === nothing
        error("Could not find ε_fid row in: $filepath")
    end

    # Find ε_1trk row: $\epsilon_{1\text{trk}}$ & value1 & value2 \\
    trk_match = match(r"\$\\epsilon_\{1\\text\{trk\}\}\$\s*&\s*(\$[^$]+\$)\s*&\s*(\$[^$]+\$)", content)
    if trk_match === nothing
        error("Could not find ε_1trk row in: $filepath")
    end

    ε_fid_2mm = parse_scientific_notation(strip(fid_match.captures[1]))
    ε_fid_20mm = parse_scientific_notation(strip(fid_match.captures[2]))
    ε_1trk_2mm = parse_scientific_notation(strip(trk_match.captures[1]))
    ε_1trk_20mm = parse_scientific_notation(strip(trk_match.captures[2]))

    return (ε_fid_2mm=ε_fid_2mm, ε_fid_20mm=ε_fid_20mm,
            ε_1trk_2mm=ε_1trk_2mm, ε_1trk_20mm=ε_1trk_20mm)
end

# =============================================================================
# LaTeX Formatting
# =============================================================================

"""
    format_latex_scientific(x::Float64)

Format number in LaTeX scientific notation.
"""
function format_latex_scientific(x::Float64)
    if x == 0.0
        return "0"
    end
    exp = floor(Int, log10(abs(x)))
    mantissa = x / 10.0^exp
    return @sprintf("\$%.2f \\times 10^{%d}\$", mantissa, exp)
end

# =============================================================================
# Document Generation
# =============================================================================

"""
    write_joint_document(data::Vector, output_path::String)

Write the joint LaTeX document with two tables.
"""
function write_joint_document(data::Vector, output_path::String)
    open(output_path, "w") do f
        # Document preamble
        println(f, "\\documentclass[11pt]{article}")
        println(f, "\\usepackage[margin=1in]{geometry}")
        println(f, "\\usepackage{booktabs}")
        println(f, "\\usepackage{amsmath}")
        println(f, "\\usepackage{caption}")
        println(f, "")
        println(f, "\\begin{document}")
        println(f, "")
        println(f, "\\section{Efficiency Tables}")
        println(f, "")

        # Table 1: Fiducial Efficiency
        println(f, "\\begin{table}[htbp]")
        println(f, "\\centering")
        println(f, "\\caption{Fiducial efficiency (\$\\epsilon_{\\text{fid}}\$) for signal and background samples.}")
        println(f, "\\label{tab:fiducial_efficiency}")
        println(f, "\\begin{tabular}{lc}")
        println(f, "\\toprule")
        println(f, "Sample & \$\\epsilon_{\\text{fid}}\$ \\\\")
        println(f, "\\midrule")

        for d in data
            row_name = ROW_NAMES[(d.isotope, d.zone)]
            # ε_fid is same for both maxd, use 2mm value
            ε_fid_str = format_latex_scientific(d.values.ε_fid_2mm)
            println(f, "$row_name & $ε_fid_str \\\\")

            # Add small spacing after average rows
            if d.zone == "weighted_mean"
                println(f, "\\addlinespace")
            end
        end

        println(f, "\\bottomrule")
        println(f, "\\end{tabular}")
        println(f, "\\end{table}")
        println(f, "")

        # Table 2: Single-Track Efficiency
        println(f, "\\begin{table}[htbp]")
        println(f, "\\centering")
        println(f, "\\caption{Single-track efficiency (\$\\epsilon_{1\\text{trk}}\$) for ion track (\$\\sigma_t = 1\$ mm) and electron track (\$\\sigma_t = 10\$ mm) conditions.}")
        println(f, "\\label{tab:single_track_efficiency}")
        println(f, "\\begin{tabular}{lcc}")
        println(f, "\\toprule")
        println(f, "Sample & Ion track & Electron track \\\\")
        println(f, "\\midrule")

        for d in data
            row_name = ROW_NAMES[(d.isotope, d.zone)]
            ε_ion_str = format_latex_scientific(d.values.ε_1trk_2mm)
            ε_elec_str = format_latex_scientific(d.values.ε_1trk_20mm)
            println(f, "$row_name & $ε_ion_str & $ε_elec_str \\\\")

            # Add small spacing after average rows
            if d.zone == "weighted_mean"
                println(f, "\\addlinespace")
            end
        end

        println(f, "\\bottomrule")
        println(f, "\\end{tabular}")
        println(f, "\\end{table}")
        println(f, "")

        # End document
        println(f, "\\end{document}")
    end

    println("Saved: $output_path")
end

# =============================================================================
# Main
# =============================================================================

function main()
    # Parse arguments
    output_path = joinpath(TABLES_DIR, "joint_tables.tex")

    for arg in ARGS
        if startswith(arg, "--")
            parts = split(arg[3:end], '=')
            if length(parts) == 2
                key, value = parts
                if key == "output"
                    output_path = String(value)
                else
                    println("Warning: Unknown argument: --$key")
                end
            end
        end
    end

    println("=" ^ 60)
    println("JOINT EFFICIENCY TABLES")
    println("=" ^ 60)
    println()

    # Read all tables
    data = []

    for tf in TABLE_FILES
        filepath = joinpath(TABLES_DIR, tf.file)
        println("Reading: $(tf.file)")

        try
            values = parse_tex_table(filepath)
            push!(data, (isotope=tf.isotope, zone=tf.zone, values=values))
            println("  ε_fid = $(format_latex_scientific(values.ε_fid_2mm))")
            println("  ε_1trk (ion) = $(format_latex_scientific(values.ε_1trk_2mm))")
            println("  ε_1trk (elec) = $(format_latex_scientific(values.ε_1trk_20mm))")
        catch e
            println("  ERROR: $e")
        end
    end

    println()

    if isempty(data)
        println("ERROR: No tables read successfully!")
        return
    end

    # Write joint document
    write_joint_document(data, output_path)

    println()
    println("=" ^ 60)
    println("Done! Output: $output_path")
    println("=" ^ 60)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
