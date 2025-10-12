#!/usr/bin/env julia

# Add project path and activate BEFORE using packages
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

# Now load packages
using ArgParse
using CSV
using DataFrames
using Statistics
using Printf
using PrettyTables

"""
    parse_commandline()

Parse command line arguments for summary tables script.
"""
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input", "-i"
            help = "Input statistics CSV file(s)"
            arg_type = String
            nargs = '+'
            required = true
        "--output", "-o"
            help = "Output file for formatted tables (supports .txt, .tex, .md)"
            arg_type = String
            default = ""
        "--format", "-f"
            help = "Output format: text, latex, markdown, html"
            arg_type = String
            default = "text"
        "--table-type", "-t"
            help = "Table type: selection, energy, length, confidence, all"
            arg_type = String
            default = "all"
        "--sort-by"
            help = "Column to sort by (e.g., frac_selected, mean_etrk_all)"
            arg_type = String
            default = "file"
        "--ascending"
            help = "Sort in ascending order (default is descending)"
            action = :store_true
        "--precision"
            help = "Number of decimal places for numerical values"
            arg_type = Int
            default = 2
    end

    return parse_args(s)
end

"""
    read_and_combine_stats(input_files::Vector{String})

Read multiple statistics CSV files and combine them into a single DataFrame.
"""
function read_and_combine_stats(input_files::Vector{String})
    all_dfs = DataFrame[]

    for file in input_files
        if !isfile(file)
            println("WARNING: File not found: $file")
            continue
        end

        df = CSV.read(file, DataFrame)

        # Add source file column if processing multiple files
        if length(input_files) > 1
            df.source_stats_file = fill(basename(file), nrow(df))
        end

        push!(all_dfs, df)
    end

    if isempty(all_dfs)
        error("No valid input files found")
    end

    # Combine all dataframes
    combined_df = vcat(all_dfs...)

    return combined_df
end

"""
    create_selection_summary_table(df::DataFrame; precision=2)

Create a summary table of selection statistics.
"""
function create_selection_summary_table(df::DataFrame; precision=2)
    # Extract relevant columns
    summary_df = select(df,
        :file => "File",
        :n_total => "Total Events",
        :n_selected => "Selected",
        :n_not_selected => "Not Selected",
        :frac_selected => (x -> round.(x * 100, digits=precision)) => "Selected %",
        :frac_not_selected => (x -> round.(x * 100, digits=precision)) => "Not Selected %"
    )

    return summary_df
end

"""
    create_energy_comparison_table(df::DataFrame; precision=2)

Create a table comparing energy statistics between selected and not selected events.
"""
function create_energy_comparison_table(df::DataFrame; precision=2)
    # Check which energy columns are available
    available_cols = names(df)

    # Start with file column
    summary_df = DataFrame(File = df.file)

    # Track Energy
    if "mean_etrk_all" in available_cols
        summary_df[!, "Track E (All)"] = round.(df.mean_etrk_all, digits=precision)
    end
    if "mean_etrk_selected" in available_cols
        summary_df[!, "Track E (Sel)"] = round.(df.mean_etrk_selected, digits=precision)
    end
    if "mean_etrk_not_selected" in available_cols
        summary_df[!, "Track E (Not)"] = round.(df.mean_etrk_not_selected, digits=precision)
    end

    # Blob 1 Energy
    if "mean_eblob1_all" in available_cols
        summary_df[!, "Blob1 E (All)"] = round.(df.mean_eblob1_all, digits=precision)
    end
    if "mean_eblob1_selected" in available_cols
        summary_df[!, "Blob1 E (Sel)"] = round.(df.mean_eblob1_selected, digits=precision)
    end
    if "mean_eblob1_not_selected" in available_cols
        summary_df[!, "Blob1 E (Not)"] = round.(df.mean_eblob1_not_selected, digits=precision)
    end

    # Blob 2 Energy
    if "mean_eblob2_all" in available_cols
        summary_df[!, "Blob2 E (All)"] = round.(df.mean_eblob2_all, digits=precision)
    end
    if "mean_eblob2_selected" in available_cols
        summary_df[!, "Blob2 E (Sel)"] = round.(df.mean_eblob2_selected, digits=precision)
    end
    if "mean_eblob2_not_selected" in available_cols
        summary_df[!, "Blob2 E (Not)"] = round.(df.mean_eblob2_not_selected, digits=precision)
    end

    return summary_df
end

"""
    create_track_length_table(df::DataFrame; precision=2)

Create a table for track length statistics.
"""
function create_track_length_table(df::DataFrame; precision=2)
    available_cols = names(df)

    # Start with file column
    summary_df = DataFrame(File = df.file)

    if "mean_trkl_all" in available_cols
        summary_df[!, "Track L (All)"] = round.(df.mean_trkl_all, digits=precision)
    end
    if "mean_trkl_selected" in available_cols
        summary_df[!, "Track L (Sel)"] = round.(df.mean_trkl_selected, digits=precision)
    end
    if "mean_trkl_not_selected" in available_cols
        summary_df[!, "Track L (Not)"] = round.(df.mean_trkl_not_selected, digits=precision)
    end

    return summary_df
end

"""
    create_confidence_table(df::DataFrame; precision=3)

Create a table for confidence statistics.
"""
function create_confidence_table(df::DataFrame; precision=3)
    available_cols = names(df)

    # Start with file column
    summary_df = DataFrame(File = df.file)

    if "mean_confidence_all" in available_cols
        summary_df[!, "Conf (All)"] = round.(df.mean_confidence_all, digits=precision)
    end
    if "mean_confidence_selected" in available_cols
        summary_df[!, "Conf (Sel)"] = round.(df.mean_confidence_selected, digits=precision)
    end
    if "mean_confidence_not_selected" in available_cols
        summary_df[!, "Conf (Not)"] = round.(df.mean_confidence_not_selected, digits=precision)
    end

    return summary_df
end

"""
    format_table(df::DataFrame, format::String, title::String="")

Format a DataFrame as a pretty table in the specified format.
"""
function format_table(df::DataFrame, format::String, title::String="")
    io = IOBuffer()

    if format == "latex"
        # LaTeX format
        pretty_table(io, df,
                    backend = Val(:latex),
                    title = title,
                    alignment = :c)
    elseif format == "markdown"
        # Markdown format
        pretty_table(io, df,
                    backend = Val(:text),
                    tf = tf_markdown,
                    title = title,
                    alignment = :c)
    elseif format == "html"
        # HTML format
        pretty_table(io, df,
                    backend = Val(:html),
                    title = title,
                    alignment = :c)
    else
        # Default text format with nice borders
        pretty_table(io, df,
                    backend = Val(:text),
                    title = title,
                    alignment = :c,
                    crop = :none,
                    tf = tf_unicode_rounded)
    end

    return String(take!(io))
end

"""
    create_comprehensive_summary(df::DataFrame; precision=2)

Create a comprehensive summary table with all key metrics.
"""
function create_comprehensive_summary(df::DataFrame; precision=2)
    available_cols = names(df)

    # Start with file name
    summary_df = DataFrame(File = df.file)

    # Add selection statistics
    if "n_total" in available_cols
        summary_df.Events = df.n_total
    end

    if "frac_selected" in available_cols
        summary_df[!, "Sel%"] = round.(df.frac_selected * 100, digits=1)
    end

    # Add mean track energy
    if "mean_etrk_selected" in available_cols && "mean_etrk_not_selected" in available_cols
        summary_df[!, "Track E(S)"] = round.(df.mean_etrk_selected, digits=0)
        summary_df[!, "Track E(N)"] = round.(df.mean_etrk_not_selected, digits=0)
        summary_df[!, "ΔE_trk"] = round.(df.mean_etrk_selected - df.mean_etrk_not_selected, digits=0)
    end

    # Add mean blob2 energy
    if "mean_eblob2_selected" in available_cols && "mean_eblob2_not_selected" in available_cols
        summary_df[!, "Blob2(S)"] = round.(df.mean_eblob2_selected, digits=0)
        summary_df[!, "Blob2(N)"] = round.(df.mean_eblob2_not_selected, digits=0)
        summary_df[!, "ΔBlob2"] = round.(df.mean_eblob2_selected - df.mean_eblob2_not_selected, digits=0)
    end

    # Add track length
    if "mean_trkl_selected" in available_cols && "mean_trkl_not_selected" in available_cols
        summary_df[!, "L(S)"] = round.(df.mean_trkl_selected, digits=0)
        summary_df[!, "L(N)"] = round.(df.mean_trkl_not_selected, digits=0)
    end

    return summary_df
end

"""
    main()

Main function for creating summary tables.
"""
function main()
    # Parse command line arguments
    args = parse_commandline()

    println("="^70)
    println("Summary Tables Generator")
    println("="^70)
    println("Input files: $(join(args["input"], ", "))")
    println("Output format: $(args["format"])")
    println("Table type: $(args["table-type"])")
    println("="^70)

    # Read and combine input files
    df = read_and_combine_stats(args["input"])
    println("\nLoaded $(nrow(df)) rows from $(length(args["input"])) file(s)")

    # Sort if requested
    if args["sort-by"] != "file" && args["sort-by"] in names(df)
        sort!(df, args["sort-by"], rev=!args["ascending"])
        println("Sorted by: $(args["sort-by"]) ($(args["ascending"] ? "ascending" : "descending"))")
    end

    # Create tables based on requested type
    tables = Dict{String, DataFrame}()

    if args["table-type"] == "all" || args["table-type"] == "selection"
        tables["Selection Summary"] = create_selection_summary_table(df, precision=args["precision"])
    end

    if args["table-type"] == "all" || args["table-type"] == "energy"
        tables["Energy Comparison (keV)"] = create_energy_comparison_table(df, precision=args["precision"])
    end

    if args["table-type"] == "all" || args["table-type"] == "length"
        tables["Track Length (mm)"] = create_track_length_table(df, precision=args["precision"])
    end

    if args["table-type"] == "all" || args["table-type"] == "confidence"
        tables["Confidence Values"] = create_confidence_table(df, precision=3)
    end

    if args["table-type"] == "all"
        tables["Comprehensive Summary"] = create_comprehensive_summary(df, precision=args["precision"])
    end

    # Format and display/save tables
    all_output = String[]

    # Sort tables by title (key) only
    sorted_tables = sort(collect(tables), by = x -> x[1])

    for (title, table_df) in sorted_tables
        println("\n" * "="^70)
        println(title)
        println("="^70)

        formatted_table = format_table(table_df, args["format"], title)
        println(formatted_table)

        push!(all_output, formatted_table)
    end

    # Save to file if requested
    if !isempty(args["output"])
        output_content = join(all_output, "\n\n")

        # Add document wrapper for LaTeX
        if args["format"] == "latex"
            output_content = """
            \\documentclass{article}
            \\usepackage{booktabs}
            \\usepackage{array}
            \\usepackage{float}
            \\begin{document}

            $output_content

            \\end{document}
            """
        elseif args["format"] == "html"
            output_content = """
            <!DOCTYPE html>
            <html>
            <head>
                <title>Summary Tables</title>
                <style>
                    table {
                        border-collapse: collapse;
                        margin: 20px;
                    }
                    th, td {
                        border: 1px solid #ddd;
                        padding: 8px;
                        text-align: center;
                    }
                    th {
                        background-color: #f2f2f2;
                    }
                    h2 {
                        margin-top: 30px;
                    }
                </style>
            </head>
            <body>
            $output_content
            </body>
            </html>
            """
        end

        open(args["output"], "w") do io
            write(io, output_content)
        end

        println("\n✓ Tables saved to: $(args["output"])")
    end

    # Print summary statistics
    println("\n" * "="^70)
    println("Summary Statistics Across All Files:")
    println("="^70)

    if "n_total" in names(df)
        total_events = sum(df.n_total)
        println("Total events: $total_events")
    end

    if "n_selected" in names(df) && "n_total" in names(df)
        total_selected = sum(df.n_selected)
        total_events = sum(df.n_total)
        overall_rate = total_selected / total_events * 100
        println("Overall selection rate: $(round(overall_rate, digits=1))%")
    end

    if "mean_etrk_all" in names(df) && "n_total" in names(df)
        # Weighted average of track energy
        weighted_etrk = sum(df.mean_etrk_all .* df.n_total) / sum(df.n_total)
        println("Weighted avg track energy: $(round(weighted_etrk, digits=1)) keV")
    end

    println("="^70)
    println("\n✓ Analysis complete!")
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end