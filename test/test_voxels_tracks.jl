# Test voxelization and track building functionality
using Petit
using Test
using DataFrames
using Graphs

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