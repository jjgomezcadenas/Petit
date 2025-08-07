module histos
#using Revise
using DataFrames
using Plots
using Statistics
using StatsBase
import Base.length
# histograms

export hist1d, hist2d, p1df, step_hist, in_range, get_histo1d, Histo1d
export scan_level, get_vals_from_sparse
export save_histo1d, load_histo1d
export save_histos, load_histos

include("util.jl")

"""Define an histogram"""
struct Histo1d
	edges::Vector{Float64}
	weights::Vector{Float64}
	centers::Vector{Float64}
end



"""Take an object of type Histogram and return a Histo1d"""
get_histo1d(h::Histogram) = Histo1d(edges(h), h.weights, centers(h))


"""Returns a Histo1d struct""" 
function hist1d(x::AbstractVector{T}, bins::AbstractVector{T}, norm=false) where {T<:Real} 
    h = fit(Histogram, x, bins)
    if norm
        h = StatsBase.normalize(h, mode=:pdf)
    end
    get_histo1d(h)
end


"""Return a Histo1d struct"""
function hist1d(x::AbstractVector{T}; 
                nbins::Integer=50, 
                xlim::Union{Nothing, Tuple{Float64, Float64}} = nothing,
                norm=false) where {T<:Real}
    
    edges = isnothing(xlim) ? nothing :
            collect(range(xlim[1], xlim[2]; length = nbins+1))
    h = isnothing(edges) ?
        fit(Histogram, x; nbins=nbins) :
        fit(Histogram, x, edges)
    #h = fit(Histogram, x; nbins=nbins)
    if norm
        h = StatsBase.normalize(h, mode=:pdf)
    end
    get_histo1d(h)
end



"""
    edges(h::Histogram)

edges of the histogram
"""
function edges(h::Histogram, index::Int64=1)
    collect(h.edges[index])
end


"""
    centers(h::Histogram)

return the centres of the histogram bins.
"""
centers(h::Histogram) = centers(edges(h))


"""
    centers(edges::Vector{<:Real})

Calculate the bin centers from a vector of bin edges
"""
function centers(edges::Vector{<:Real})
    edges[1:end-1] + .-(@view(edges[2:end]), @view(edges[1:end-1])) / 2
end


"""
    hist_weights

return histogram weights of a set of data.
"""
function hist_weights(edges::Vector{<:Real})
  function get_weights(y::SubArray{<:Real})
    histo = fit(Histogram, y, edges)
    return histo.weights
  end
  return get_weights
end


"""
return 1D histogram and plot axes for data x and number of bins nbins.
"""
function hist1d(x::AbstractVector{T}, nbins::Integer, xs::String;
                norm::Bool=false, 
                datap::Bool=true,
                markersize::Int64=3, 
                fwl::Bool=false,
                label::String="", 
                legend::Bool=false) where {T<:Real} 

    h = hist1d(x, nbins, norm)
    h, plot_histogram(h, datap, markersize, fwl, label, legend)
    
end


"""
return 1D histogram and plot axes for data x and number of bins nbins.
"""
function hist1d(x::AbstractVector{T},  bins::AbstractVector{T}, xs::String;
                norm::Bool=false, 
                datap::Bool=true,
                markersize::Int64=3, 
                fwl::Bool=false,
                label::String="", 
                legend::Bool=false) where {T<:Real} 

    h = hist1d(x, bins, norm)
    h, plot_histogram(h, datap, markersize, fwl, label, xs, legend)
    
end


function plot_histogram(h::Histo1d, datap::Bool,
                        markersize::Int64, 
                        fwl::Bool,
                        label::String,
                        xs::String, 
                        legend::Bool)
    xg = h.centers
    yg = eltype(xg).(h.weights)
    
    if datap    
        p  = scatter(xg, yg, yerr=sqrt.(yg),
                     markersize=markersize, label=label, legend=legend)
        if fwl
            p = plot!(p, xg, yg, yerr=sqrt.(yg), fmt=:png,
                     linewidth=1, label=label, legend=legend)
        end
    else
        p = plot(h, xlabel=xs, fmt=:png, yl="frequency")
    end
    xlabel!(xs)
    ylabel!("frequency")

    return p
end


"""
    digitize(x, bins)

    Return the indices of the bins to which each value in input array belongs.
"""
function digitize(x::Vector{<:Real}, bins::LinRange{<:Real})
    return searchsortedlast.(Ref(bins), x)
end
# Notice that in this case the bins are expected to be sorted
function digitize(x::Vector{<:Real}, bins::Vector{<:Real})
    return searchsortedlast.(Ref(bins), x)
end


"""
    hist2d(x::Vector{T}, y::Vector{T}, nbins::Integer,
           xl::String, yl::String,
           xmin::T=typemin(T), xmax::T=typemax(T),
           ymin::T=typemin(T), ymax::T=typemax(T)) where T <: Real

ggiven data x, y, labels xl and yl and limits for each axis.
"""
function hist2d(x::Vector{T},y::Vector{T}, nbins::Integer,
                xl::String, yl::String,
                xmin::T=typemin(T), xmax::T=typemax(T),
                ymin::T=typemin(T), ymax::T=typemax(T); title="") where T <: Real
    function select_data(xval::T, yval::T)
        range_bound(xmin, xmax, OpenBound)(xval) .&
        range_bound(ymin, ymax, OpenBound)(yval)
    end
    mask = select_data.(x, y)
    data = (y[mask], x[mask])
    h    = fit(Histogram, data, nbins=nbins)
    xe   = h.edges[1]
    ye   = h.edges[2]
    hm   = heatmap(xe, ye, h.weights)
    xlabel!(xl)
    ylabel!(yl)
    if title != ""
        title!(title)
    end
    return h, hm
end


"""
    fmean_std

returns the weighted mean and std of an array given a minimum
bin weight for consideration.
"""
function fmean_std(x::Vector{<:Real}, min_prop::Float64=0.1)
    mask_func = weights -> weights .> min_prop * maximum(weights)
    function filt_mean(w::SubArray{<:Real})
        mask = mask_func(w)
        return [mean_and_std(x[mask], FrequencyWeights(w[mask]), corrected=true)]
    end
    return filt_mean
end


"""
    p1df(x, y, nbins; ybin_width, ymin, ymax, min_proportion)

return a profile DataFrame. This is a DF in which the variable y is
histogrammed as a function of the average of variable x in each bin.
The keyword arguments allow for a filtering in the y variable by
range and min_proportion of the maximum in the bin.
TODO: Protect against fine binning that results in zeros.
"""
function p1df(x::Vector{<:Real}, y::Vector{<:Real}, nbins::Integer;
              ybin_width::Real=0.1, ymin::Real=minimum(y), ymax::Real=maximum(y),
              min_proportion::Real=0.0)
    df           = DataFrame(:y => y)
    x_upper      = nextfloat(maximum(x))
    bin_edges    = LinRange(minimum(x), x_upper, nbins + 1)
    df[!, "bin"] = digitize(x, bin_edges)
    bin_centers  = centers(collect(bin_edges))
    bin_width    = .-(@view(bin_edges[2:end]), @view(bin_edges[1:end-1]))
    ## Bin in y so filtering can be done on bin content
    y_upper      = nextfloat(ymax)
    y_bins       = collect(ymin:ybin_width:y_upper)
    y_centers    = centers(y_bins)
    y_weights    = combine(groupby(df, :bin),
                           :y => hist_weights(y_bins) => :yweights,
                           ungroup = false)
    bin_stats    = fmean_std(y_centers, min_proportion)
    ndf = combine(y_weights,
                  :yweights => bin_stats => [:y_mean, :y_std])
    ndf = ndf[:, [:y_mean, :y_std]]
    ndf[!, :x_mean] = bin_centers
    ndf[!, :x_std]  = bin_width / 2

    p1 = plot(ndf.x_mean,ndf.y_mean, yerror=ndf.y_std, fmt = :png,
              shape = :circle, color = :black, legend=false)
    return ndf, p1
end


"""
    step_hist(data::Vector{<:Real};
              nbins::Int=50,
              xlim::Union{Nothing,Tuple{Float64,Float64}}=nothing,
              logy::Bool=false,
              xlabel::String="",
              ylabel::String="")

Plot an unfilled step‑style histogram, adding a small offset to avoid zeros on a log scale.

# Keywords
- `nbins`: Number of bins.
- `xlim`: Optional `(min, max)` bin range.
- `logy`: If `true`, use log10 on the y‑axis.
- `xlabel`, `ylabel`: Axis labels.
"""
function step_hist(data::Vector{<:Real};
                   nbins::Int = 50,
                   xlim::Union{Nothing, Tuple{Float64, Float64}} = nothing,
                   logy::Bool = false,
                   xlabel::String = "",
                   ylabel::String = "",
                   title::String="")

    # 1. Compute histogram
    edges = isnothing(xlim) ? nothing :
            collect(range(xlim[1], xlim[2]; length = nbins+1))
    h = isnothing(edges) ?
        fit(Histogram, data; nbins=nbins) :
        fit(Histogram, data, edges)

    bin_edges = h.edges[1]   # nbins+1 edges
    counts    = h.weights    # nbins counts

    # 2. Build staircase coords
    x_step = repeat(bin_edges, inner=2)
    # pad counts with zeros at ends, then add +1 to every y‑value
    y_pad  = vcat(0.0, counts, 0.0) .+ 1.0
    y_step = repeat(y_pad, inner=2)[2:end-1]

    # 3. Plot the step histogram
    plt = plot(x_step, y_step;
               seriestype = :step,
               linecolor  = :black,
               lw         = 1,
               xlabel     = xlabel,
               ylabel     = ylabel,
               title     = title,
               legend     = false)

    # 4. Optional log scale
    if logy
        yaxis!(plt, :log10)
    end

    # 5. Draw vertical boundary lines
    baseline = logy ? minimum(y_step) : 0.0

    # Left boundary
    e1 = bin_edges[1]
    h1 = y_pad[2]   # first bin’s height after +1
    plot!(plt, [e1, e1], [baseline, h1]; linecolor=:black, lw=1)

    # Right boundary
    en = bin_edges[end]
    hn = y_pad[end-1]  # last bin’s height after +1
    plot!(plt, [en, en], [baseline, hn]; linecolor=:black, lw=1)

    return get_histo1d(h), plt
end


"""
    save_histo1d(h::Histo1d, filename::String)

Save a Histo1d histogram to a text file.

# Arguments
- `h::Histo1d`: The histogram to save
- `filename::String`: Path to the output file

# File Format
The file contains:
- Line 1: Number of bins
- Line 2: Edges (space-separated)
- Line 3: Weights (space-separated)
- Line 4: Centers (space-separated)

# Example
```julia
h = hist1d([1.0, 2.0, 3.0, 4.0, 5.0]; nbins=5)
save_histo1d(h, "my_histogram.txt")
```
"""
function save_histo1d(h::Histo1d, filename::String)
    open(filename, "w") do io
        # Write number of bins
        println(io, length(h.weights))
        
        # Write edges (space-separated)
        println(io, join(h.edges, " "))
        
        # Write weights (space-separated)
        println(io, join(h.weights, " "))
        
        # Write centers (space-separated)
        println(io, join(h.centers, " "))
    end
end


"""
    load_histo1d(filename::String) -> Histo1d

Load a Histo1d histogram from a text file created by save_histo1d.

# Arguments
- `filename::String`: Path to the input file

# Returns
- `Histo1d`: The loaded histogram

# Example
```julia
h = load_histo1d("my_histogram.txt")
```
"""
function load_histo1d(filename::String)
    lines = readlines(filename)
    
    if length(lines) != 4
        error("Invalid histogram file format. Expected 4 lines, got $(length(lines))")
    end
    
    # Parse number of bins (not strictly needed but good for validation)
    n_bins = parse(Int, strip(lines[1]))
    
    # Parse edges
    edges = parse.(Float64, split(strip(lines[2])))
    
    # Parse weights
    weights = parse.(Float64, split(strip(lines[3])))
    
    # Parse centers
    centers = parse.(Float64, split(strip(lines[4])))
    
    # Validate dimensions
    if length(weights) != n_bins
        error("Mismatch: expected $n_bins weights, got $(length(weights))")
    end
    
    if length(edges) != n_bins + 1
        error("Mismatch: expected $(n_bins + 1) edges, got $(length(edges))")
    end
    
    if length(centers) != n_bins
        error("Mismatch: expected $n_bins centers, got $(length(centers))")
    end
    
    return Histo1d(edges, weights, centers)
end


"""
    save_histos(histos::Dict{String, Histo1d}, filename::String)

Save a collection of named histograms to a single file.

# Arguments
- `histos::Dict{String, Histo1d}`: Dictionary mapping histogram names to Histo1d objects
- `filename::String`: Path to the output file

# File Format
The file contains:
- Line 1: Number of histograms
- For each histogram:
  - Line: Histogram name
  - Line: Number of bins
  - Line: Edges (space-separated)
  - Line: Weights (space-separated)
  - Line: Centers (space-separated)

# Example
```julia
h1 = hist1d([1.0, 2.0, 3.0]; nbins=3)
h2 = hist1d([4.0, 5.0, 6.0]; nbins=2)
histos = Dict("energy" => h1, "position" => h2)
save_histos(histos, "my_histograms.txt")
```
"""
function save_histos(histos::Dict{String, Histo1d}, filename::String)
    open(filename, "w") do io
        # Write number of histograms
        println(io, length(histos))
        
        # Write each histogram
        for (name, h) in histos
            # Write histogram name
            println(io, name)
            
            # Write number of bins
            println(io, length(h.weights))
            
            # Write edges (space-separated)
            println(io, join(h.edges, " "))
            
            # Write weights (space-separated)
            println(io, join(h.weights, " "))
            
            # Write centers (space-separated)  
            println(io, join(h.centers, " "))
        end
    end
end


"""
    save_histos(histos::Vector{Histo1d}, names::Vector{String}, filename::String)

Save a collection of histograms with given names to a single file.

# Arguments
- `histos::Vector{Histo1d}`: Vector of Histo1d objects
- `names::Vector{String}`: Vector of names for the histograms
- `filename::String`: Path to the output file

# Example
```julia
h1 = hist1d([1.0, 2.0, 3.0]; nbins=3)
h2 = hist1d([4.0, 5.0, 6.0]; nbins=2)
save_histos([h1, h2], ["energy", "position"], "my_histograms.txt")
```
"""
function save_histos(histos::Vector{Histo1d}, names::Vector{String}, filename::String)
    if length(histos) != length(names)
        error("Number of histograms ($(length(histos))) must match number of names ($(length(names)))")
    end
    
    histo_dict = Dict(name => hist for (name, hist) in zip(names, histos))
    save_histos(histo_dict, filename)
end


"""
    load_histos(filename::String) -> Dict{String, Histo1d}

Load a collection of named histograms from a file created by save_histos.

# Arguments
- `filename::String`: Path to the input file

# Returns
- `Dict{String, Histo1d}`: Dictionary mapping histogram names to Histo1d objects

# Example
```julia
histos = load_histos("my_histograms.txt")
energy_hist = histos["energy"]
position_hist = histos["position"]
```
"""
function load_histos(filename::String)
    lines = readlines(filename)
    
    if length(lines) < 1
        error("Invalid histogram collection file: file is empty")
    end
    
    # Parse number of histograms
    n_histos = parse(Int, strip(lines[1]))
    
    if n_histos == 0
        return Dict{String, Histo1d}()
    end
    
    # Expected total lines: 1 (count) + n_histos * 5 (name + 4 data lines per histogram)
    expected_lines = 1 + n_histos * 5
    if length(lines) != expected_lines
        error("Invalid histogram collection file format. Expected $expected_lines lines, got $(length(lines))")
    end
    
    histos = Dict{String, Histo1d}()
    line_idx = 2  # Start after the count line
    
    for i in 1:n_histos
        # Parse histogram name
        name = strip(lines[line_idx])
        line_idx += 1
        
        # Parse number of bins
        n_bins = parse(Int, strip(lines[line_idx]))
        line_idx += 1
        
        # Parse edges
        edges = parse.(Float64, split(strip(lines[line_idx])))
        line_idx += 1
        
        # Parse weights
        weights = parse.(Float64, split(strip(lines[line_idx])))
        line_idx += 1
        
        # Parse centers
        centers = parse.(Float64, split(strip(lines[line_idx])))
        line_idx += 1
        
        # Validate dimensions
        if length(weights) != n_bins
            error("Histogram '$name': expected $n_bins weights, got $(length(weights))")
        end
        
        if length(edges) != n_bins + 1
            error("Histogram '$name': expected $(n_bins + 1) edges, got $(length(edges))")
        end
        
        if length(centers) != n_bins
            error("Histogram '$name': expected $n_bins centers, got $(length(centers))")
        end
        
        # Create histogram and add to collection
        histos[name] = Histo1d(edges, weights, centers)
    end
    
    return histos
end

end # module 