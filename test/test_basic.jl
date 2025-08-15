# Test basic functionality
using Petit
using Test
using DataFrames
using Graphs

@testset "Basic Functionality" begin
    @test isdefined(Petit, :Tracks)
    @test isdefined(Petit, :AnalysisResults)
    @test isdefined(Petit, :euclidean_distance)
    @test isdefined(Petit, :hist1d)
    @test isdefined(Petit, :voxelize_hits)
    @test isdefined(Petit, :build_tracks)
    @test isdefined(Petit, :select_events)
    @test isdefined(Petit, :analysis_loop)
    @test isdefined(Petit, :voxel_distances)
    @test isdefined(Petit, :voxel_closest_distance)
    @test isdefined(Petit, :voxel_energy)
    @test isdefined(Petit, :hits_per_event)
    @test isdefined(Petit, :hits_per_all_events)
end

@testset "Euclidean Distance" begin
    # Test distance calculation
    @test euclidean_distance(0.0, 0.0, 0.0, 1.0, 1.0, 1.0) ≈ sqrt(3)
    @test euclidean_distance(0.0, 0.0, 0.0, 0.0, 0.0, 0.0) ≈ 0.0
    @test euclidean_distance(1.0, 0.0, 0.0, 4.0, 0.0, 0.0) ≈ 3.0
    @test euclidean_distance(0.0, 0.0, 0.0, 3.0, 4.0, 0.0) ≈ 5.0  # 3-4-5 triangle
end

@testset "Tracks Structure" begin
    # Test Tracks constructor
    test_df = DataFrame(x=[1.0, 2.0], y=[1.0, 2.0], z=[1.0, 2.0], energy=[0.5, 0.5])
    test_graph = SimpleGraph(2)
    test_components = [[1, 2]]
    
    tracks = Tracks(test_df, test_graph, test_components)
    @test tracks isa Tracks
    @test nrow(tracks.voxels) == 2
    @test nv(tracks.graph) == 2
    @test length(tracks.components) == 1
end