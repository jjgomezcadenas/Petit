# Test plotting functionality
using Petit
using Test
using DataFrames
using Graphs

@testset "Plot Functions" begin
    # Create test data for plotting
    test_df = DataFrame(
        event_id = [1, 1, 1, 2, 2, 2],
        x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        y = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        z = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        energy = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    )
    
    # Create test Tracks object
    test_graph = SimpleGraph(3)
    add_edge!(test_graph, 1, 2)
    add_edge!(test_graph, 2, 3)
    test_components = [[1, 2, 3]]
    test_voxels = DataFrame(x=[0.0, 1.0, 2.0], y=[0.0, 1.0, 2.0], 
                           z=[0.0, 1.0, 2.0], energy=[0.1, 0.2, 0.3])
    test_track = Tracks(test_voxels, test_graph, test_components)
    
    # Test that plot functions are defined and callable
    @test isdefined(Petit, :plot_hits_trk)
    @test isdefined(Petit, :plot_hits_evt)
    @test isdefined(Petit, :plot_hits)
    
    # We could test that they return plot objects if needed:
    # p1 = plot_hits_trk(test_track; nbins=50)
    # @test p1 isa Plots.Plot
    
    # p2 = plot_hits_evt(test_df, 1; nbins=50)
    # @test p2 isa Plots.Plot
    
    # p3 = plot_hits(test_voxels; nbins=50)
    # @test p3 isa Plots.Plot
    
    @testset "Make Tracks" begin
        # Create test hit data with multiple events
        test_hits = DataFrame(
            event_id = [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3],
            x = [0.0, 0.5, 5.0, 10.0, 15.0, 0.0, 1.0, 20.0, 0.0, 0.1, 0.2, 0.3],
            y = [0.0, 0.5, 5.0, 10.0, 15.0, 0.0, 1.0, 20.0, 0.0, 0.1, 0.2, 0.3],
            z = [0.0, 0.5, 5.0, 10.0, 15.0, 0.0, 1.0, 20.0, 0.0, 0.1, 0.2, 0.3],
            energy = [0.5, 0.6, 0.005, 0.008, 0.7, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4]  # MeV units
        )
        
        # Test with first 3 events
        results = make_tracks(test_hits, 1:3)
        @test results isa DataFrame
        @test "event_id" in names(results)
        @test "ntracks" in names(results)
        @test "n1trk" in names(results)
        @test "n2trk" in names(results)
        @test "n3trk" in names(results)
        
        # Should have results for 3 events
        @test nrow(results) == 3
        
        # Check event IDs
        @test results.event_id == [1, 2, 3]
        
        # Check that ntracks values are reasonable
        @test all(results.ntracks .>= 0)
        @test all(results.n1trk .>= 0)
        @test all(results.n2trk .>= 0)
        @test all(results.n3trk .>= 0)
        
        # Test with subset of events
        results_subset = make_tracks(test_hits, [1, 3])
        @test nrow(results_subset) == 2
        @test results_subset.event_id == [1, 3]
        
        # Test with single event
        results_single = make_tracks(test_hits, [2])
        @test nrow(results_single) == 1
        @test results_single.event_id == [2]
        
        # Test with custom parameters
        results_custom = make_tracks(test_hits, 1:2; 
                                    voxel_size_mm=5.0, 
                                    max_distance_mm=10.0, 
                                    energy_threshold_kev=10.0)
        @test nrow(results_custom) == 2
        
        # Test with empty event list
        results_empty = make_tracks(test_hits, Int[])
        @test nrow(results_empty) == 0
        @test "event_id" in names(results_empty)
        @test "ntracks" in names(results_empty)
    end
end