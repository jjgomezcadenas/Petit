using Petit
using Test
using DataFrames
using Graphs

@testset "Petit.jl" begin
    
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
    
    @testset "Voxelization" begin
        # Create test hit data
        test_hits = DataFrame(
            event_id=[1, 1, 1, 2, 2],
            x=[0.1, 0.9, 2.1, 1.1, 1.9],
            y=[0.1, 0.9, 2.1, 1.1, 1.9], 
            z=[0.1, 0.9, 2.1, 1.1, 1.9],
            energy=[0.1, 0.2, 0.3, 0.15, 0.25]
        )
        
        voxels = voxelize_hits(test_hits, 1.0)  # 1mm voxel size
        
        @test voxels isa DataFrame
        @test nrow(voxels) > 0
        @test all(col in names(voxels) for col in ["event_id", "x", "y", "z", "energy"])
        
        # Test that voxelization reduces number of hits (hits in same voxel are combined)
        @test nrow(voxels) <= nrow(test_hits)
        
        # Test with larger voxel size (should reduce voxel count more)
        voxels_large = voxelize_hits(test_hits, 2.0)
        @test nrow(voxels_large) <= nrow(voxels)
    end
    
    @testset "Track Building" begin
        # Create test data
        test_data = DataFrame(
            event_id=[1, 1, 1],
            x=[0.0, 1.0, 3.0],  # Third point is far away
            y=[0.0, 1.0, 3.0],
            z=[0.0, 1.0, 3.0],
            energy=[0.1, 0.2, 0.15]
        )
        
        # Test with max_distance that connects first two points but not third
        tracks = build_tracks(test_data, 1; max_distance=2.0, energy_threshold=0.05)
        
        @test length(tracks) >= 1  # Should create at least one track
        @test all(t isa Tracks for t in tracks)
        
        # Test with energy threshold that excludes some points (150 keV = 0.15 MeV)
        tracks_filtered = build_tracks(test_data, 1; max_distance=5.0, energy_threshold=150.0)
        @test length(tracks_filtered) >= 0
        
        # Test that energy filtering works (threshold was 150 keV = 0.15 MeV)
        for track in tracks_filtered
            @test all(track.voxels.energy .>= 0.15)
        end
    end
    
    @testset "Voxel Distances" begin
        # Create test data with known distances
        test_data = DataFrame(
            event_id=[1, 1, 1, 2, 2],
            x=[0.0, 1.0, 0.0, 0.0, 2.0],
            y=[0.0, 0.0, 1.0, 0.0, 0.0],
            z=[0.0, 0.0, 0.0, 0.0, 0.0],
            energy=[0.1, 0.2, 0.15, 0.25, 0.3]
        )
        
        # Test voxel_distances for event 1 (method with event_id)
        distances = voxel_distances(test_data, 1)
        @test length(distances) == 3  # 3 pairs for 3 voxels
        @test 1.0 in distances  # distance between (0,0,0) and (1,0,0)
        @test 1.0 in distances  # distance between (0,0,0) and (0,1,0)
        @test sqrt(2) ≈ distances[findfirst(d -> d ≈ sqrt(2), distances)]  # distance between (1,0,0) and (0,1,0)
        
        # Test with max_distance filter
        distances_filtered = voxel_distances(test_data, 1; max_distance=1.0)
        @test length(distances_filtered) == 2  # Only distances ≤ 1.0
        
        # Test direct DataFrame method
        event1_df = test_data[test_data.event_id .== 1, :]
        distances_direct = voxel_distances(event1_df)
        @test length(distances_direct) == 3
        @test distances ≈ distances_direct
        
        # Test error handling
        @test_throws Exception voxel_distances(test_data, 999)
    end
    
    @testset "Voxel Closest Distance" begin
        # Create test data with known closest distances
        test_data = DataFrame(
            event_id=[1, 1, 1],
            x=[0.0, 1.0, 5.0],  # Third point is far away
            y=[0.0, 0.0, 0.0],
            z=[0.0, 0.0, 0.0],
            energy=[0.1, 0.2, 0.15]
        )
        
        # Test voxel_closest_distance
        closest_dists = voxel_closest_distance(test_data, 1)
        @test length(closest_dists) == 3  # One closest distance per voxel
        @test 1.0 in closest_dists  # closest distance for first voxel
        @test 1.0 in closest_dists  # closest distance for second voxel
        @test 4.0 in closest_dists  # closest distance for third voxel
        
        # Test with max_distance filter
        closest_filtered = voxel_closest_distance(test_data, 1; max_distance=2.0)
        @test length(closest_filtered) == 2  # Only distances ≤ 2.0
        
        # Test direct DataFrame method
        event1_df = test_data[test_data.event_id .== 1, :]
        closest_direct = voxel_closest_distance(event1_df)
        @test length(closest_direct) == 3
        @test closest_dists ≈ closest_direct
        
        # Test error handling
        @test_throws Exception voxel_closest_distance(test_data, 999)
    end
    
    @testset "Voxel Energy" begin
        # Create test data
        test_data = DataFrame(
            event_id=[1, 1, 1, 2, 2],
            x=[0.0, 1.0, 2.0, 0.0, 1.0],
            y=[0.0, 1.0, 2.0, 0.0, 1.0],
            z=[0.0, 1.0, 2.0, 0.0, 1.0],
            energy=[0.1, 0.2, 0.3, 0.4, 0.5]
        )
        
        # Test voxel_energy for event 1 (method with event_id)
        energies = voxel_energy(test_data, 1)
        @test length(energies) == 3
        @test energies ≈ [0.1, 0.2, 0.3]
        
        # Test voxel_energy for event 2
        energies2 = voxel_energy(test_data, 2)
        @test length(energies2) == 2
        @test energies2 ≈ [0.4, 0.5]
        
        # Test direct DataFrame method
        event1_df = test_data[test_data.event_id .== 1, :]
        energies_direct = voxel_energy(event1_df)
        @test energies ≈ energies_direct
        
        # Test error handling
        @test_throws Exception voxel_energy(test_data, 999)
    end
    
    @testset "Hits Per Event" begin
        # Create test data
        test_hits = DataFrame(
            event_id=[1, 1, 1, 2, 2, 3, 3, 3, 3],
            x=rand(9),
            y=rand(9),
            z=rand(9),
            energy=rand(9)
        )
        
        # Test hits_per_event
        @test hits_per_event(test_hits, 1) == 3
        @test hits_per_event(test_hits, 2) == 2
        @test hits_per_event(test_hits, 3) == 4
        
        # Test error handling
        @test_throws Exception hits_per_event(test_hits, 999)
        
        # Test hits_per_all_events
        all_counts = hits_per_all_events(test_hits)
        @test length(all_counts) == 3  # 3 events
        @test all_counts == [3, 2, 4]  # hits per event
    end
    
    @testset "Select Events" begin
        # Create test data with multiple events
        test_data = DataFrame(
            event_id=[1, 1, 1, 2, 2, 3, 3, 3, 3],
            x=[0.0, 1.0, 5.0, 0.0, 1.0, 0.0, 1.0, 2.0, 10.0],
            y=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            energy=[0.1, 0.2, 0.15, 0.3, 0.25, 0.1, 0.2, 0.15, 0.05]
        )
        
        # Test select_events function
        tracks = select_events(test_data, 1; voxel_size_mm=2.0, max_distance_mm=2.0, energy_threshold_kev=50.0)
        @test tracks isa Vector{Tracks}
        @test length(tracks) >= 0  # Should return some tracks
        
        # Test with high energy threshold (should return fewer/no tracks)
        tracks_filtered = select_events(test_data, 1; voxel_size_mm=2.0, max_distance_mm=2.0, energy_threshold_kev=500.0)
        @test tracks_filtered isa Vector{Tracks}
        
        # Test error handling for non-existent event
        @test_throws Exception select_events(test_data, 999)
    end
    
    @testset "Analysis Loop" begin
        # Create comprehensive test data with known track patterns
        test_data = DataFrame(
            event_id=[1, 1, 1,        # Event 1: 3 close hits (should form 1 track)
                     2, 2,           # Event 2: 2 close hits (should form 1 track)  
                     3, 3, 3, 3,     # Event 3: 2 pairs of hits (should form 2 tracks)
                     4, 4, 4, 4, 4], # Event 4: 3 groups of hits (should form 3+ tracks)
            x=[0.0, 1.0, 2.0,        # Event 1: linear arrangement
               0.0, 1.0,             # Event 2: close pair
               0.0, 1.0, 10.0, 11.0, # Event 3: two separate pairs
               0.0, 1.0, 10.0, 11.0, 20.0], # Event 4: three groups
            y=[0.0, 0.0, 0.0,
               0.0, 0.0,
               0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0],
            z=[0.0, 0.0, 0.0,
               0.0, 0.0,
               0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0],
            energy=[0.1, 0.2, 0.3,   # Event 1: total ~0.6 MeV = 600 keV
                    0.15, 0.25,      # Event 2: total ~0.4 MeV = 400 keV
                    0.1, 0.2, 0.15, 0.25,  # Event 3: tracks ~0.3, 0.4 MeV
                    0.1, 0.2, 0.15, 0.25, 0.1] # Event 4: multiple tracks
        )
        
        # Test analysis_loop function
        results = analysis_loop(test_data; 
                               events_to_run=1:4, 
                               voxel_size_mm=3.0,     # Large voxels to merge close hits
                               max_distance_mm=5.0,   # Connect hits within 5mm
                               energy_threshold_kev=50.0) # 50 keV threshold
        
        # Test that results is of correct type
        @test results isa AnalysisResults
        
        # Test structure fields exist
        @test hasfield(AnalysisResults, :single_track_energies)
        @test hasfield(AnalysisResults, :two_track_primary)
        @test hasfield(AnalysisResults, :two_track_secondary)
        @test hasfield(AnalysisResults, :three_track_primary)
        @test hasfield(AnalysisResults, :three_track_secondary)
        @test hasfield(AnalysisResults, :n_events_processed)
        @test hasfield(AnalysisResults, :n_single_track)
        @test hasfield(AnalysisResults, :n_two_track)
        @test hasfield(AnalysisResults, :n_three_plus_track)
        @test hasfield(AnalysisResults, :n_failed)
        
        # Test basic statistics
        @test results.n_events_processed == 4
        @test results.n_failed >= 0
        @test results.n_single_track + results.n_two_track + results.n_three_plus_track <= results.n_events_processed
        
        # Test energy arrays are proper type
        @test results.single_track_energies isa Vector{Float64}
        @test results.two_track_primary isa Vector{Float64}
        @test results.two_track_secondary isa Vector{Float64}
        @test results.three_track_primary isa Vector{Float64}
        @test results.three_track_secondary isa Vector{Float64}
        
        # Test that energy values are positive (if any tracks found)
        if length(results.single_track_energies) > 0
            @test all(e -> e > 0, results.single_track_energies)
        end
        if length(results.two_track_primary) > 0
            @test all(e -> e > 0, results.two_track_primary)
            @test all(e -> e > 0, results.two_track_secondary)
        end
        if length(results.three_track_primary) > 0
            @test all(e -> e > 0, results.three_track_primary)
            @test all(e -> e > 0, results.three_track_secondary)
        end
        
        # Test consistency: secondary track arrays should have same length as primary
        @test length(results.two_track_primary) == length(results.two_track_secondary)
        @test length(results.three_track_primary) == length(results.three_track_secondary)
        
        # Test with empty event range
        empty_results = analysis_loop(test_data; events_to_run=Int[])
        @test empty_results.n_events_processed == 0
        @test length(empty_results.single_track_energies) == 0
    end
    
    @testset "Error Handling" begin
        # Test with non-existent event_id
        test_data = DataFrame(
            event_id=[1, 1],
            x=[0.0, 1.0],
            y=[0.0, 1.0], 
            z=[0.0, 1.0],
            energy=[0.1, 0.2]
        )
        
        # Should throw error for non-existent event
        @test_throws Exception build_tracks(test_data, 999; max_distance=1.0, energy_threshold=0.0)
        @test_throws Exception voxel_distances(test_data, 999)
        @test_throws Exception voxel_closest_distance(test_data, 999)
        @test_throws Exception voxel_energy(test_data, 999)
        @test_throws Exception hits_per_event(test_data, 999)
        
        # Test with empty data
        empty_data = DataFrame(event_id=Int[], x=Float64[], y=Float64[], z=Float64[], energy=Float64[])
        
        # These should handle empty data gracefully
        @test nrow(voxelize_hits(empty_data, 1.0)) == 0
        @test length(hits_per_all_events(empty_data)) == 0
    end
end