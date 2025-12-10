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

@testset "Voxelize Event" begin
    # Create test hit data with multiple events
    test_hits = DataFrame(
        event_id=[1, 1, 1, 1, 2, 2, 2],
        x=[0.1, 0.9, 2.1, 2.3, 1.1, 1.9, 3.5],
        y=[0.1, 0.9, 2.1, 2.4, 1.1, 1.9, 3.2],
        z=[0.1, 0.9, 2.1, 2.2, 1.1, 1.9, 3.8],
        energy=[0.1, 0.2, 0.3, 0.15, 0.25, 0.35, 0.4]
    )

    # Test basic functionality - voxelize event 1
    voxels_evt1 = voxelize_event(test_hits, 1, 1.0)

    @test voxels_evt1 isa DataFrame
    @test nrow(voxels_evt1) > 0
    @test all(col in names(voxels_evt1) for col in ["event_id", "x", "y", "z", "energy", "electrons"])

    # All voxels should belong to event 1
    @test all(voxels_evt1.event_id .== 1)

    # Should only include hits from event 1 (4 hits)
    @test nrow(voxels_evt1) <= 4

    # Test voxelize event 2
    voxels_evt2 = voxelize_event(test_hits, 2, 1.0)

    @test voxels_evt2 isa DataFrame
    @test nrow(voxels_evt2) > 0
    @test all(voxels_evt2.event_id .== 2)
    @test nrow(voxels_evt2) <= 3  # Event 2 has 3 hits

    # Test energy aggregation - hits in same voxel should have combined energy
    test_same_voxel = DataFrame(
        event_id=[1, 1, 1],
        x=[0.1, 0.2, 0.3],  # All in same voxel (voxel size = 1.0)
        y=[0.1, 0.2, 0.3],
        z=[0.1, 0.2, 0.3],
        energy=[0.1, 0.2, 0.3]
    )

    voxels_combined = voxelize_event(test_same_voxel, 1, 1.0)
    @test nrow(voxels_combined) == 1  # All hits in one voxel
    @test voxels_combined.energy[1] ≈ 0.6  # Sum of all energies
    @test voxels_combined.electrons[1] == 60000  # energy * 1e5 = 0.6 * 1e5 = 60000

    # Test voxel centers are computed correctly
    # Voxel index i corresponds to center (i + 0.5) * voxel_size
    # For hit at x=0.1 with voxel_size=1.0: floor(0.1/1.0) = 0, center = 0.5
    @test voxels_combined.x[1] ≈ 0.5
    @test voxels_combined.y[1] ≈ 0.5
    @test voxels_combined.z[1] ≈ 0.5

    # Test with larger voxel size
    voxels_large = voxelize_event(test_hits, 1, 2.0)
    @test nrow(voxels_large) <= nrow(voxels_evt1)  # Larger voxels → fewer voxels

    # Test voxel center calculation with different position
    test_negative = DataFrame(
        event_id=[1, 1],
        x=[-1.5, -1.2],  # Both in voxel -2
        y=[0.0, 0.0],
        z=[0.0, 0.0],
        energy=[0.1, 0.2]
    )

    voxels_neg = voxelize_event(test_negative, 1, 1.0)
    @test nrow(voxels_neg) == 1
    # floor(-1.5/1.0) = -2, center = (-2 + 0.5) * 1.0 = -1.5
    @test voxels_neg.x[1] ≈ -1.5
    @test voxels_neg.energy[1] ≈ 0.3

    # Test with hits in different voxels
    test_different_voxels = DataFrame(
        event_id=[1, 1, 1],
        x=[0.5, 1.5, 2.5],  # Three different voxels
        y=[0.5, 1.5, 2.5],
        z=[0.5, 1.5, 2.5],
        energy=[0.1, 0.2, 0.3]
    )

    voxels_diff = voxelize_event(test_different_voxels, 1, 1.0)
    @test nrow(voxels_diff) == 3  # Should stay 3 voxels
    @test sum(voxels_diff.energy) ≈ 0.6  # Total energy conserved

    # Test consistency with voxelize_hits
    # voxelize_event(df, evt_id, size) should give same result as
    # filtering voxelize_hits(df, size) for that event
    voxels_from_hits = voxelize_hits(test_hits, 1.0)
    voxels_from_event = voxelize_event(test_hits, 1, 1.0)

    voxels_hits_filtered = voxels_from_hits[voxels_from_hits.event_id .== 1, :]

    # Should have same number of voxels
    @test nrow(voxels_from_event) == nrow(voxels_hits_filtered)

    # Should have same total energy (allowing for floating point errors)
    @test sum(voxels_from_event.energy) ≈ sum(voxels_hits_filtered.energy)

    # Test error handling - nonexistent event should throw error
    @test_throws Exception voxelize_event(test_hits, 999, 1.0)
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

@testset "Make Tracks" begin
    # Create test data with a linear track and a separate isolated voxel
    test_data = DataFrame(
        event_id=[1, 1, 1, 1, 1],
        x=[0.0, 1.0, 2.0, 3.0, 10.0],  # Last point is isolated
        y=[0.0, 0.0, 0.0, 0.0, 10.0],
        z=[0.0, 0.0, 0.0, 0.0, 10.0],
        energy=[0.5, 0.3, 0.4, 0.2, 0.6]  # Isolated point has highest individual energy
    )

    # Test make_tracks with parameters that connect first 4 points but not the 5th
    tracks = make_tracks(test_data, 1; max_dist=1.5, energy_thr=0.1)

    @test length(tracks) == 2  # Should create 2 tracks
    @test all(t isa Tracks for t in tracks)

    # Tracks should be sorted by total energy (first track has more total energy)
    track_energies = [sum(t.voxels.energy) for t in tracks]
    @test track_energies[1] >= track_energies[2]
    @test track_energies[1] ≈ 1.4  # Sum of first 4 points
    @test track_energies[2] ≈ 0.6  # Isolated point

    # Test with higher energy threshold that excludes low energy voxels
    # energy_thr is in keV, our data is in MeV, so 250 keV = 0.25 MeV
    tracks_filtered = make_tracks(test_data, 1; max_dist=1.5, energy_thr=250.0)
    total_voxels = sum(nrow(t.voxels) for t in tracks_filtered)
    @test total_voxels == 4  # Only 4 voxels have energy >= 0.25 MeV

    # Test with no valid tracks (high threshold and small max_dist)
    tracks_empty = make_tracks(test_data, 1; max_dist=0.1, energy_thr=10000.0)
    @test length(tracks_empty) == 0
end

@testset "Find Track Extremes" begin
    # Test Case 1: Linear track with clear endpoints
    linear_data = DataFrame(
        x=[0.0, 1.0, 2.0, 3.0],
        y=[0.0, 0.0, 0.0, 0.0],
        z=[0.0, 0.0, 0.0, 0.0],
        energy=[0.1, 0.2, 0.3, 0.4]
    )

    # Create a linear graph
    linear_graph = SimpleGraph(4)
    add_edge!(linear_graph, 1, 2)
    add_edge!(linear_graph, 2, 3)
    add_edge!(linear_graph, 3, 4)

    linear_track = Tracks(linear_data, linear_graph, [[1, 2, 3, 4]])

    extreme1, extreme2, path, confidence = find_track_extremes(linear_track)

    @test extreme1 == 1 || extreme1 == 4  # Should be one endpoint
    @test extreme2 == 1 || extreme2 == 4  # Should be other endpoint
    @test extreme1 != extreme2
    @test length(path) == 4  # Path should include all 4 vertices
    @test path == [1, 2, 3, 4] || path == [4, 3, 2, 1]

    # Test Case 2: Branched track (Y-shape)
    branch_data = DataFrame(
        x=[0.0, 1.0, 2.0, 2.0],
        y=[0.0, 0.0, 1.0, -1.0],
        z=[0.0, 0.0, 0.0, 0.0],
        energy=[0.1, 0.2, 0.3, 0.4]
    )

    branch_graph = SimpleGraph(4)
    add_edge!(branch_graph, 1, 2)
    add_edge!(branch_graph, 2, 3)
    add_edge!(branch_graph, 2, 4)

    branch_track = Tracks(branch_data, branch_graph, [[1, 2, 3, 4]])

    extreme1_b, extreme2_b, path_b, confidence_b = find_track_extremes(branch_track)

    # Should find endpoints (vertices with degree 1)
    @test extreme1_b in [1, 3, 4]
    @test extreme2_b in [1, 3, 4]
    @test extreme1_b != extreme2_b

    # Test Case 3: Circular track (no clear endpoints)
    circle_data = DataFrame(
        x=[1.0, 0.0, -1.0, 0.0],
        y=[0.0, 1.0, 0.0, -1.0],
        z=[0.0, 0.0, 0.0, 0.0],
        energy=[0.1, 0.2, 0.3, 0.4]
    )

    circle_graph = SimpleGraph(4)
    add_edge!(circle_graph, 1, 2)
    add_edge!(circle_graph, 2, 3)
    add_edge!(circle_graph, 3, 4)
    add_edge!(circle_graph, 4, 1)

    circle_track = Tracks(circle_data, circle_graph, [[1, 2, 3, 4]])

    extreme1_c, extreme2_c, path_c, confidence_c = find_track_extremes(circle_track)

    # Should find furthest apart vertices (1 and 3 are furthest)
    @test (extreme1_c == 1 && extreme2_c == 3) || (extreme1_c == 3 && extreme2_c == 1)
    @test length(path_c) >= 2  # Should have a valid path

    # Test Case 4: Single vertex
    single_data = DataFrame(x=[0.0], y=[0.0], z=[0.0], energy=[0.1])
    single_graph = SimpleGraph(1)
    single_track = Tracks(single_data, single_graph, [[1]])

    extreme1_s, extreme2_s, path_s, confidence_s = find_track_extremes(single_track)

    @test extreme1_s == 1
    @test extreme2_s == 1
    @test path_s == [1]

    # Test Case 5: Empty track
    empty_data = DataFrame(x=Float64[], y=Float64[], z=Float64[], energy=Float64[])
    empty_graph = SimpleGraph(0)
    empty_track = Tracks(empty_data, empty_graph, Vector{Int}[])

    extreme1_e, extreme2_e, path_e, confidence_e = find_track_extremes(empty_track)

    @test isnothing(extreme1_e)
    @test isnothing(extreme2_e)
    @test isempty(path_e)

    # Test confidence scores
    @test confidence > 0.0  # Should have some confidence
    @test confidence_b > 0.0
    @test confidence_c >= 0.0  # Circular tracks might have low confidence
    @test confidence_s >= 0.0
    @test confidence_e == 0.0  # Empty track should have zero confidence

    # Test multiple calls consistency
    for _ in 1:5
        extreme1_test, extreme2_test, path_test, conf_test = find_track_extremes(linear_track)
        @test (extreme1_test, extreme2_test) == (extreme1, extreme2) || (extreme1_test, extreme2_test) == (extreme2, extreme1)
        @test conf_test == confidence
    end
end

@testset "Improved Track Extremes Algorithms" begin
    # Test track with sharp turn (L-shape)
    sharp_turn_data = DataFrame(
        x=[0.0, 1.0, 2.0, 3.0, 3.0, 3.0],  # Horizontal then vertical
        y=[0.0, 0.0, 0.0, 0.0, 1.0, 2.0],  # Sharp 90° turn
        z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        energy=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    )

    sharp_turn_graph = SimpleGraph(6)
    for i in 1:5
        add_edge!(sharp_turn_graph, i, i+1)
    end
    sharp_turn_track = Tracks(sharp_turn_data, sharp_turn_graph, [[1, 2, 3, 4, 5, 6]])

    # Test different methods
    topo_result = find_track_extremes_improved(sharp_turn_track; method=:topology)
    curv_result = find_track_extremes_improved(sharp_turn_track; method=:curvature)
    comb_result = find_track_extremes_improved(sharp_turn_track; method=:combined)

    @test length(topo_result) == 4  # (extreme1, extreme2, path, confidence)
    @test length(curv_result) == 4
    @test length(comb_result) == 4

    # All methods should find the true endpoints (vertices 1 and 6)
    @test (topo_result[1] == 1 && topo_result[2] == 6) || (topo_result[1] == 6 && topo_result[2] == 1)
    @test topo_result[4] > 0.8  # High confidence for clear endpoints

    # Test curvature calculation
    curvatures = calculate_vertex_curvatures(sharp_turn_track)
    @test length(curvatures) == 6
    @test curvatures[1] == 0.0  # Endpoint should have zero curvature
    @test curvatures[6] == 0.0  # Endpoint should have zero curvature
    @test curvatures[4] > 0.5   # Sharp turn vertex should have high curvature

    # Test path length calculation
    path_length = calculate_path_length(sharp_turn_track, [1, 2, 3, 4, 5, 6])
    expected_length = 5.0  # 5 unit segments
    @test abs(path_length - expected_length) < 0.1

    # Test vertex proximity check
    @test are_vertices_too_close(sharp_turn_track, 1, 2, 1.5)   # Adjacent vertices are close
    @test !are_vertices_too_close(sharp_turn_track, 1, 6, 1.5)  # Endpoints are far apart
end

@testset "Walk Track From Extremes" begin
    # Create a test linear track
    linear_data = DataFrame(
        x=[0.0, 1.0, 2.0, 3.0, 4.0],
        y=[0.0, 0.0, 0.0, 0.0, 0.0],
        z=[0.0, 0.0, 0.0, 0.0, 0.0],
        energy=[0.1, 0.2, 0.3, 0.4, 0.5]
    )

    linear_graph = SimpleGraph(5)
    add_edge!(linear_graph, 1, 2)
    add_edge!(linear_graph, 2, 3)
    add_edge!(linear_graph, 3, 4)
    add_edge!(linear_graph, 4, 5)

    linear_track = Tracks(linear_data, linear_graph, [[1, 2, 3, 4, 5]])

    result = walk_track_from_extremes(linear_track)

    # Test extremes
    start_voxel, end_voxel = result.extremes
    @test !isnothing(start_voxel)
    @test !isnothing(end_voxel)
    @test start_voxel.x == 0.0 || start_voxel.x == 4.0
    @test end_voxel.x == 0.0 || end_voxel.x == 4.0
    @test start_voxel.x != end_voxel.x

    # Test path
    @test length(result.path_indices) == 5
    @test result.path_indices == [1, 2, 3, 4, 5] || result.path_indices == [5, 4, 3, 2, 1]

    # Test path voxels
    @test nrow(result.path_voxels) == 5
    @test result.path_voxels.energy == linear_data.energy[result.path_indices]

    # Test total length
    @test result.total_length ≈ 4.0  # 4 segments of length 1.0

    # Test confidence score
    @test result.confidence > 0.0  # Should have positive confidence

    # Test with empty track
    empty_data = DataFrame(x=Float64[], y=Float64[], z=Float64[], energy=Float64[])
    empty_graph = SimpleGraph(0)
    empty_track = Tracks(empty_data, empty_graph, Vector{Int}[])

    empty_result = walk_track_from_extremes(empty_track)

    @test isnothing(empty_result.extremes[1])
    @test isnothing(empty_result.extremes[2])
    @test isempty(empty_result.path_indices)
    @test nrow(empty_result.path_voxels) == 0
    @test empty_result.total_length == 0.0
    @test empty_result.confidence == 0.0  # Empty track should have zero confidence

    # Test with single vertex
    single_data = DataFrame(x=[1.5], y=[2.5], z=[3.5], energy=[0.75])
    single_graph = SimpleGraph(1)
    single_track = Tracks(single_data, single_graph, [[1]])

    single_result = walk_track_from_extremes(single_track)

    @test single_result.extremes[1] == single_result.extremes[2]
    @test single_result.path_indices == [1]
    @test nrow(single_result.path_voxels) == 1
    @test single_result.total_length == 0.0
    @test single_result.confidence == 1.0  # Single vertex should have perfect confidence
end

@testset "Find Path BFS" begin
    # Create a simple graph for testing
    g = SimpleGraph(5)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 4)
    add_edge!(g, 2, 5)  # Branch at vertex 2

    # Test direct path
    path = find_path_bfs(g, 1, 4)
    @test path == [1, 2, 3, 4]

    # Test path to branch
    path_branch = find_path_bfs(g, 1, 5)
    @test path_branch == [1, 2, 5]

    # Test reverse path
    path_reverse = find_path_bfs(g, 4, 1)
    @test path_reverse == [4, 3, 2, 1]

    # Test same vertex
    path_same = find_path_bfs(g, 2, 2)
    @test path_same == [2]

    # Test disconnected vertices
    g_disconnected = SimpleGraph(4)
    add_edge!(g_disconnected, 1, 2)
    add_edge!(g_disconnected, 3, 4)  # 3-4 not connected to 1-2

    path_disconnected = find_path_bfs(g_disconnected, 1, 4)
    @test isempty(path_disconnected)
end

@testset "Energy in Spheres Around Extremes" begin
    # Create a linear test track with known positions and energies
    linear_data = DataFrame(
        x=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        y=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        energy=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    )

    # Create linear graph: 1-2-3-4-5-6
    linear_graph = SimpleGraph(6)
    add_edge!(linear_graph, 1, 2)
    add_edge!(linear_graph, 2, 3)
    add_edge!(linear_graph, 3, 4)
    add_edge!(linear_graph, 4, 5)
    add_edge!(linear_graph, 5, 6)

    linear_track = Tracks(linear_data, linear_graph, [[1, 2, 3, 4, 5, 6]])

    @testset "Basic functionality" begin
        # Test with radius that includes only extreme voxels
        result = energy_in_spheres_around_extremes(linear_track, 0.5)

        @test result.start_sphere_energy ≈ 0.1  # Only voxel 1
        @test result.end_sphere_energy ≈ 0.6    # Only voxel 6
        @test result.start_voxel_count == 1
        @test result.end_voxel_count == 1
        @test result.start_center == (0.0, 0.0, 0.0)
        @test result.end_center == (5.0, 0.0, 0.0)
    end

    @testset "Larger radius including neighbors" begin
        # Test with radius that includes neighboring voxels
        result = energy_in_spheres_around_extremes(linear_track, 1.1)

        @test result.start_sphere_energy ≈ 0.3   # Voxels 1 and 2: 0.1 + 0.2
        @test result.end_sphere_energy ≈ 1.1     # Voxels 5 and 6: 0.5 + 0.6
        @test result.start_voxel_count == 2
        @test result.end_voxel_count == 2
    end

    @testset "Very large radius including all voxels" begin
        # Test with radius that includes all voxels
        result = energy_in_spheres_around_extremes(linear_track, 10.0)

        total_energy = sum(linear_data.energy)  # 2.1
        @test result.start_sphere_energy ≈ total_energy
        @test result.end_sphere_energy ≈ total_energy
        @test result.start_voxel_count == 6
        @test result.end_voxel_count == 6
    end

    @testset "With walk_result parameter" begin
        # Test the version that takes walk_result as parameter
        walk_result = walk_track_from_extremes(linear_track)
        result = energy_in_spheres_around_extremes(linear_track, walk_result, 1.1)

        @test result.start_sphere_energy ≈ 0.3   # Voxels 1 and 2
        @test result.end_sphere_energy ≈ 1.1     # Voxels 5 and 6
        @test result.start_voxel_count == 2
        @test result.end_voxel_count == 2
    end

    @testset "Edge cases" begin
        # Test with empty track
        empty_data = DataFrame(x=Float64[], y=Float64[], z=Float64[], energy=Float64[])
        empty_graph = SimpleGraph(0)
        empty_track = Tracks(empty_data, empty_graph, Vector{Int}[])

        empty_result = energy_in_spheres_around_extremes(empty_track, 1.0)
        @test empty_result.start_sphere_energy == 0.0
        @test empty_result.end_sphere_energy == 0.0
        @test empty_result.start_voxel_count == 0
        @test empty_result.end_voxel_count == 0
        @test isnan(empty_result.start_center[1])
        @test isnan(empty_result.end_center[1])

        # Test with single voxel track
        single_data = DataFrame(x=[2.5], y=[3.5], z=[1.5], energy=[0.75])
        single_graph = SimpleGraph(1)
        single_track = Tracks(single_data, single_graph, [[1]])

        single_result = energy_in_spheres_around_extremes(single_track, 0.1)
        @test single_result.start_sphere_energy ≈ 0.75
        @test single_result.end_sphere_energy ≈ 0.75  # Same voxel is both start and end
        @test single_result.start_voxel_count == 1
        @test single_result.end_voxel_count == 1
        @test single_result.start_center == (2.5, 3.5, 1.5)
        @test single_result.end_center == (2.5, 3.5, 1.5)

        # Test with zero radius
        zero_result = energy_in_spheres_around_extremes(linear_track, 0.0)
        @test zero_result.start_sphere_energy ≈ 0.1  # Only exact position
        @test zero_result.end_sphere_energy ≈ 0.6
        @test zero_result.start_voxel_count == 1
        @test zero_result.end_voxel_count == 1
    end

    @testset "3D track with varying coordinates" begin
        # Create a track with voxels in different 3D positions
        data_3d = DataFrame(
            x=[0.0, 1.0, 2.0, 3.0, 4.0],
            y=[0.0, 1.0, 2.0, 3.0, 4.0],
            z=[0.0, 1.0, 2.0, 3.0, 4.0],
            energy=[0.2, 0.3, 0.4, 0.5, 0.6]
        )

        graph_3d = SimpleGraph(5)
        add_edge!(graph_3d, 1, 2)
        add_edge!(graph_3d, 2, 3)
        add_edge!(graph_3d, 3, 4)
        add_edge!(graph_3d, 4, 5)

        track_3d = Tracks(data_3d, graph_3d, [[1, 2, 3, 4, 5]])

        # Distance from (0,0,0) to (1,1,1) is sqrt(3) ≈ 1.73
        # Distance from (4,4,4) to (3,3,3) is sqrt(3) ≈ 1.73
        result_3d = energy_in_spheres_around_extremes(track_3d, 1.8)

        @test result_3d.start_sphere_energy ≈ 0.5   # Voxels 1 and 2: 0.2 + 0.3
        @test result_3d.end_sphere_energy ≈ 1.1     # Voxels 4 and 5: 0.5 + 0.6
        @test result_3d.start_voxel_count == 2
        @test result_3d.end_voxel_count == 2
        @test result_3d.start_center == (0.0, 0.0, 0.0)
        @test result_3d.end_center == (4.0, 4.0, 4.0)
    end

    @testset "Blob1/Blob2 format (new return values)" begin
        # Test that blob1 always has higher energy than blob2
        result = energy_in_spheres_around_extremes(linear_track, 1.1)

        # In this case: start_sphere_energy = 0.3, end_sphere_energy = 1.1
        # So blob1 should be the end sphere (higher energy)
        @test result.blob1_energy ≈ 1.1      # Higher energy sphere
        @test result.blob2_energy ≈ 0.3      # Lower energy sphere
        @test result.blob1_voxel_count == 2  # End sphere voxel count
        @test result.blob2_voxel_count == 2  # Start sphere voxel count
        @test result.blob1_center == (5.0, 0.0, 0.0)  # End center
        @test result.blob2_center == (0.0, 0.0, 0.0)  # Start center

        # Test case where start sphere has higher energy
        # Create track with higher energy at start
        reversed_data = DataFrame(
            x=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            y=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            energy=[0.8, 0.7, 0.3, 0.2, 0.1, 0.05]  # Higher energy at start
        )
        reversed_track = Tracks(reversed_data, linear_graph, [[1, 2, 3, 4, 5, 6]])

        result_rev = energy_in_spheres_around_extremes(reversed_track, 1.1)

        # start_sphere_energy = 0.8 + 0.7 = 1.5, end_sphere_energy = 0.1 + 0.05 = 0.15
        # So blob1 should be the start sphere (higher energy)
        @test result_rev.blob1_energy ≈ 1.5       # Higher energy sphere (start)
        @test result_rev.blob2_energy ≈ 0.15      # Lower energy sphere (end)
        @test result_rev.blob1_voxel_count == 2   # Start sphere voxel count
        @test result_rev.blob2_voxel_count == 2   # End sphere voxel count
        @test result_rev.blob1_center == (0.0, 0.0, 0.0)  # Start center
        @test result_rev.blob2_center == (5.0, 0.0, 0.0)  # End center

        # Test equal energy case (blob1 should be start sphere by convention)
        equal_data = DataFrame(
            x=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            y=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            energy=[0.3, 0.3, 0.3, 0.3, 0.3, 0.3]  # Equal energies
        )
        equal_track = Tracks(equal_data, linear_graph, [[1, 2, 3, 4, 5, 6]])

        result_equal = energy_in_spheres_around_extremes(equal_track, 1.1)

        # Both spheres have same energy: 0.6, blob1 should be start (by >= condition)
        @test result_equal.blob1_energy ≈ 0.6     # Start sphere energy
        @test result_equal.blob2_energy ≈ 0.6     # End sphere energy
        @test result_equal.blob1_center == (0.0, 0.0, 0.0)  # Start center
        @test result_equal.blob2_center == (5.0, 0.0, 0.0)  # End center

        # Verify backward compatibility: blob values should match start/end appropriately
        @test result_equal.start_sphere_energy == result_equal.blob1_energy
        @test result_equal.end_sphere_energy == result_equal.blob2_energy
    end
end