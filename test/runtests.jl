using Petit
using Test
using DataFrames
using Graphs
using StatsBase

@testset "Petit.jl" begin
    # Include all test files
    #include("test_basic.jl")
    #include("test_voxels_tracks.jl")
    #include("test_histograms.jl")
    #include("test_analysis.jl")
    include("test_io.jl")
    #include("test_plots.jl")
end