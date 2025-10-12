# Comprehensive tests for analysis and event loop functions
using Petit
using Test
using DataFrames
using Graphs
using HDF5
using JSON

@testset "Determine Events to Process" begin
    # Test data: event IDs [1, 3, 5, 7, 9] (not sequential)
    unique_event_ids = [1, 3, 5, 7, 9]

    @testset "Process all events (nothing)" begin
        result = determine_events_to_process(nothing, unique_event_ids)
        @test result == unique_event_ids
        @test length(result) == 5
    end

    @testset "Process first N events (Integer)" begin
        # Process first 3 events
        result = determine_events_to_process(3, unique_event_ids)
        @test result == [1, 3, 5]
        @test length(result) == 3

        # Request more than available
        result = determine_events_to_process(10, unique_event_ids)
        @test result == unique_event_ids
        @test length(result) == 5

        # Process exactly available number
        result = determine_events_to_process(5, unique_event_ids)
        @test result == unique_event_ids

        # Process zero events
        result = determine_events_to_process(0, unique_event_ids)
        @test isempty(result)
    end

    @testset "Process specific events (Range)" begin
        # Range that includes some existing IDs
        result = determine_events_to_process(1:5, unique_event_ids)
        @test result == [1, 3, 5]  # Only these exist in unique_event_ids
        @test length(result) == 3

        # Range that includes all existing IDs
        result = determine_events_to_process(1:10, unique_event_ids)
        @test result == unique_event_ids

        # Range with no matching IDs
        result = determine_events_to_process(10:20, unique_event_ids)
        @test isempty(result)

        # Empty range
        result = determine_events_to_process(5:4, unique_event_ids)
        @test isempty(result)
    end

    @testset "Process specific events (Vector)" begin
        # Vector with some existing IDs
        result = determine_events_to_process([1, 5, 9, 11], unique_event_ids)
        @test result == [1, 5, 9]  # Only these exist
        @test length(result) == 3

        # Vector with all existing IDs
        result = determine_events_to_process([9, 1, 5, 3, 7], unique_event_ids)
        @test Set(result) == Set(unique_event_ids)  # Order might differ

        # Vector with no matching IDs
        result = determine_events_to_process([2, 4, 6], unique_event_ids)
        @test isempty(result)

        # Empty vector
        result = determine_events_to_process(Int[], unique_event_ids)
        @test isempty(result)
    end

    @testset "Error handling" begin
        # Invalid input type
        @test_throws ArgumentError determine_events_to_process("invalid", unique_event_ids)
        @test_throws ArgumentError determine_events_to_process(1.5, unique_event_ids)
        @test_throws ArgumentError determine_events_to_process(Dict(), unique_event_ids)
    end

    @testset "Edge cases" begin
        # Empty unique_event_ids
        empty_ids = Int[]
        result = determine_events_to_process(nothing, empty_ids)
        @test isempty(result)

        result = determine_events_to_process(5, empty_ids)
        @test isempty(result)

        result = determine_events_to_process([1, 2, 3], empty_ids)
        @test isempty(result)

        # Single event ID
        single_id = [42]
        result = determine_events_to_process(nothing, single_id)
        @test result == [42]

        result = determine_events_to_process(1, single_id)
        @test result == [42]

        result = determine_events_to_process([42], single_id)
        @test result == [42]

        result = determine_events_to_process([41, 43], single_id)
        @test isempty(result)
    end
end

@testset "Analysis Loop Comprehensive" begin
    # Create comprehensive test data with known track patterns
    test_data = DataFrame(
        event_id = vcat(
            fill(1, 5),    # Event 1: Single track (5 voxels)
            fill(2, 8),    # Event 2: Two tracks (4+4 voxels)
            fill(3, 12),   # Event 3: Three tracks (5+4+3 voxels)
            fill(4, 3),    # Event 4: Single track (3 voxels)
            fill(5, 6)     # Event 5: Two tracks (3+3 voxels)
        ),
        x = vcat(
            # Event 1: Linear track
            [0.0, 1.0, 2.0, 3.0, 4.0],
            # Event 2: Two separate clusters
            [0.0, 1.0, 2.0, 3.0, 20.0, 21.0, 22.0, 23.0],
            # Event 3: Three separate clusters
            [0.0, 1.0, 2.0, 3.0, 4.0, 50.0, 51.0, 52.0, 53.0, 100.0, 101.0, 102.0],
            # Event 4: Short linear track
            [10.0, 11.0, 12.0],
            # Event 5: Two small clusters
            [30.0, 31.0, 32.0, 70.0, 71.0, 72.0]
        ),
        y = vcat(
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),
        z = vcat(
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),
        energy = vcat(
            [0.5, 0.4, 0.6, 0.3, 0.2],      # Event 1: Total = 2.0 MeV
            [0.3, 0.3, 0.2, 0.2, 0.4, 0.4, 0.3, 0.3],  # Event 2: Track1=1.0, Track2=1.4 MeV
            [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.3, 0.2, 0.2, 0.1, 0.1, 0.1],  # Event 3: Various
            [0.3, 0.4, 0.3],                # Event 4: Total = 1.0 MeV
            [0.2, 0.2, 0.2, 0.3, 0.3, 0.3]  # Event 5: Track1=0.6, Track2=0.9 MeV
        )
    )

    @testset "Basic functionality" begin
        results = analysis_loop(test_data; events_to_run=nothing,
                               voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)

        # Check return type
        @test results isa AnalysisResults

        # Check that all events were processed
        @test results.n_events_processed == 5

        # Check track count statistics
        @test results.n_single_track == 2    # Events 1 and 4
        @test results.n_two_track == 2       # Events 2 and 5
        @test results.n_three_plus_track == 1  # Event 3
        @test results.n_failed == 0

        # Verify total adds up
        @test results.n_single_track + results.n_two_track + results.n_three_plus_track + results.n_failed == results.n_events_processed

        # Check DataFrame structures
        @test results.single_track isa DataFrame
        @test results.two_track_primary isa DataFrame
        @test results.two_track_secondary isa DataFrame
        @test results.three_track_primary isa DataFrame
        @test results.three_track_secondary isa DataFrame

        # Check DataFrame columns
        expected_cols = [:event_id, :energy, :x, :y, :z]
        for df in [results.single_track, results.two_track_primary, results.two_track_secondary,
                   results.three_track_primary, results.three_track_secondary]
            @test Set(names(df)) == Set(string.(expected_cols))
        end
    end

    @testset "Event selection parameters" begin
        # Test processing first 3 events
        results = analysis_loop(test_data; events_to_run=3,
                               voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test results.n_events_processed == 3
        @test results.n_single_track == 1  # Only event 1
        @test results.n_two_track == 1     # Only event 2
        @test results.n_three_plus_track == 1  # Only event 3

        # Test processing specific events
        results = analysis_loop(test_data; events_to_run=[1, 4, 5],
                               voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test results.n_events_processed == 3
        @test results.n_single_track == 2  # Events 1 and 4
        @test results.n_two_track == 1     # Event 5
        @test results.n_three_plus_track == 0

        # Test processing with range
        results = analysis_loop(test_data; events_to_run=2:4,
                               voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test results.n_events_processed == 3
        @test results.n_single_track == 1  # Event 4
        @test results.n_two_track == 1     # Event 2 (note: event 5 not in range)
        @test results.n_three_plus_track == 1  # Event 3
    end

    @testset "Parameter effects" begin
        # Test with higher energy threshold
        results = analysis_loop(test_data; events_to_run=nothing,
                               voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=300.0)
        # Some low-energy voxels should be filtered out
        @test results.n_events_processed == 5
        # Track counts might change due to energy filtering

        # Test with larger max_distance (should merge more voxels into fewer tracks)
        results = analysis_loop(test_data; events_to_run=nothing,
                               voxel_size_mm=0.5, max_distance_mm=25.0, energy_threshold_kev=1.0)
        @test results.n_events_processed == 5
        # Some two-track events might become single-track due to larger connection distance

        # Test with larger voxel size
        results = analysis_loop(test_data; events_to_run=nothing,
                               voxel_size_mm=2.0, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test results.n_events_processed == 5
    end

    @testset "DataFrame content validation" begin
        results = analysis_loop(test_data; events_to_run=[1, 2],
                               voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)

        # Check single track data (Event 1)
        @test nrow(results.single_track) == 5  # Event 1 has 5 voxels
        @test all(results.single_track.event_id .== 1)
        # Event 1 total energy is 0.5+0.4+0.6+0.3+0.2 = 2.0 MeV = 2000 keV
        expected_energy = 1000 * (0.5 + 0.4 + 0.6 + 0.3 + 0.2)  # Convert to keV
        @test all(results.single_track.energy .≈ expected_energy)

        # Check two track data (Event 2)
        primary_rows = nrow(results.two_track_primary)
        secondary_rows = nrow(results.two_track_secondary)
        @test primary_rows + secondary_rows == 8  # Event 2 has 8 total voxels
        @test all(results.two_track_primary.event_id .== 2)
        @test all(results.two_track_secondary.event_id .== 2)

        # Energies should be consistent within each track type
        if primary_rows > 0
            @test length(unique(results.two_track_primary.energy)) == 1
        end
        if secondary_rows > 0
            @test length(unique(results.two_track_secondary.energy)) == 1
        end
    end

    @testset "Empty and edge cases" begin
        # Test with no events to process
        results = analysis_loop(test_data; events_to_run=Int[],
                               voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test results.n_events_processed == 0
        @test nrow(results.single_track) == 0
        @test nrow(results.two_track_primary) == 0
        @test nrow(results.two_track_secondary) == 0
        @test nrow(results.three_track_primary) == 0
        @test nrow(results.three_track_secondary) == 0

        # Test with non-existent events
        results = analysis_loop(test_data; events_to_run=[99, 100],
                               voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test results.n_events_processed == 0

        # Test with empty DataFrame
        empty_data = DataFrame(event_id=Int[], x=Float64[], y=Float64[], z=Float64[], energy=Float64[])
        results = analysis_loop(empty_data; events_to_run=nothing,
                               voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test results.n_events_processed == 0
        @test results.n_single_track == 0
        @test results.n_two_track == 0
        @test results.n_three_plus_track == 0
        @test results.n_failed == 0
    end
end

@testset "Analysis Loop Single Track" begin
    # Use same test data as comprehensive analysis loop tests
    test_data = DataFrame(
        event_id = vcat(fill(1, 5), fill(2, 8), fill(3, 12), fill(4, 3), fill(5, 6)),
        x = vcat(
            [0.0, 1.0, 2.0, 3.0, 4.0],
            [0.0, 1.0, 2.0, 3.0, 20.0, 21.0, 22.0, 23.0],
            [0.0, 1.0, 2.0, 3.0, 4.0, 50.0, 51.0, 52.0, 53.0, 100.0, 101.0, 102.0],
            [10.0, 11.0, 12.0],
            [30.0, 31.0, 32.0, 70.0, 71.0, 72.0]
        ),
        y = vcat(
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),
        z = vcat(
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ),
        energy = vcat(
            [0.5, 0.4, 0.6, 0.3, 0.2],
            [0.3, 0.3, 0.2, 0.2, 0.4, 0.4, 0.3, 0.3],
            [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.3, 0.2, 0.2, 0.1, 0.1, 0.1],
            [0.3, 0.4, 0.3],
            [0.2, 0.2, 0.2, 0.3, 0.3, 0.3]
        )
    )

    @testset "Basic functionality" begin
        tracks = analysis_loop_single_track(test_data; events_to_run=nothing,
                                           voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)

        # Check return type
        @test tracks isa Vector{Tracks}

        # Should only return single-track events (events 1 and 4)
        @test length(tracks) == 2

        # Check that each element is a Tracks object
        @test all(t isa Tracks for t in tracks)

        # Check that tracks have the expected structure
        for track in tracks
            @test track.voxels isa DataFrame
            @test track.graph isa SimpleGraph
            @test track.components isa Vector{Vector{Int}}
            @test nrow(track.voxels) > 0
        end
    end

    @testset "Track properties" begin
        tracks = analysis_loop_single_track(test_data; events_to_run=[1, 4],
                                           voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)

        @test length(tracks) == 2

        # Find tracks by number of voxels (Event 1 has 5, Event 4 has 3)
        track_sizes = [nrow(t.voxels) for t in tracks]
        @test Set(track_sizes) == Set([5, 3])

        # Check that voxels have required columns
        expected_cols = ["x", "y", "z", "energy"]
        for track in tracks
            for col in expected_cols
                @test col in names(track.voxels)
            end
        end
    end

    @testset "Event filtering" begin
        # Process only first 2 events (1 and 2) - only event 1 is single track
        tracks = analysis_loop_single_track(test_data; events_to_run=2,
                                           voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test length(tracks) == 1  # Only event 1 is single track
        @test nrow(tracks[1].voxels) == 5  # Event 1 has 5 voxels

        # Process events [2, 3, 5] - none are single track
        tracks = analysis_loop_single_track(test_data; events_to_run=[2, 3, 5],
                                           voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test length(tracks) == 0  # No single-track events

        # Process specific single-track events
        tracks = analysis_loop_single_track(test_data; events_to_run=[1, 4],
                                           voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test length(tracks) == 2  # Both are single track
    end

    @testset "Parameter effects" begin
        # Test with larger max_distance - might merge some tracks
        tracks = analysis_loop_single_track(test_data; events_to_run=nothing,
                                           voxel_size_mm=0.5, max_distance_mm=25.0, energy_threshold_kev=1.0)
        # Some multi-track events might become single-track with larger connection distance
        @test length(tracks) >= 2  # At least the originally single-track events

        # Test with high energy threshold
        tracks = analysis_loop_single_track(test_data; events_to_run=nothing,
                                           voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=500.0)
        # Some tracks might be eliminated or changed due to energy filtering
        @test length(tracks) <= 2
    end

    @testset "Edge cases" begin
        # Test with no events
        tracks = analysis_loop_single_track(test_data; events_to_run=Int[],
                                           voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test length(tracks) == 0

        # Test with non-existent events
        tracks = analysis_loop_single_track(test_data; events_to_run=[99, 100],
                                           voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test length(tracks) == 0

        # Test with empty DataFrame
        empty_data = DataFrame(event_id=Int[], x=Float64[], y=Float64[], z=Float64[], energy=Float64[])
        tracks = analysis_loop_single_track(empty_data; events_to_run=nothing,
                                           voxel_size_mm=0.5, max_distance_mm=1.5, energy_threshold_kev=1.0)
        @test length(tracks) == 0
    end
end

@testset "Event Loop Helper Functions" begin
    @testset "Validate Event Loop Parameters" begin
        # Test valid parameters
        @test_nowarn validate_event_loop_parameters(100, 5.0, 10.0, 1.0)
        @test_nowarn validate_event_loop_parameters(1, 0.1, 0.1, 0.0)

        # Test invalid events_to_run
        @test_throws ArgumentError validate_event_loop_parameters(0, 5.0, 10.0, 1.0)
        @test_throws ArgumentError validate_event_loop_parameters(-5, 5.0, 10.0, 1.0)

        # Test invalid voxel_size_mm
        @test_throws ArgumentError validate_event_loop_parameters(100, 0.0, 10.0, 1.0)
        @test_throws ArgumentError validate_event_loop_parameters(100, -1.0, 10.0, 1.0)

        # Test invalid max_distance_mm
        @test_throws ArgumentError validate_event_loop_parameters(100, 5.0, 0.0, 1.0)
        @test_throws ArgumentError validate_event_loop_parameters(100, 5.0, -2.0, 1.0)

        # Test invalid energy_threshold_kev
        @test_throws ArgumentError validate_event_loop_parameters(100, 5.0, 10.0, -1.0)

        # Test edge case: energy_threshold_kev = 0 should be valid
        @test_nowarn validate_event_loop_parameters(100, 5.0, 10.0, 0.0)
    end

    @testset "Load and Prepare Data" begin
        # Create a temporary HDF5 file for testing
        temp_dir = mktempdir()
        test_file = joinpath(temp_dir, "test_data.h5")

        # Create test data
        test_hits = DataFrame(
            event_id = [1, 1, 2, 2, 3, 3],
            x = [-50.0, 50.0, -200.0, 200.0, 0.0, 0.0],  # Event 2 fails fiducial cuts
            y = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            z = [150.0, 150.0, 150.0, 150.0, 150.0, 150.0],  # All pass z cut
            energy = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        )

        # Save test data to HDF5 with the expected structure (MC group)
        h5open(test_file, "w") do file
            mc_group = create_group(file, "MC")
            hits_group = create_group(mc_group, "hits")
            for col in names(test_hits)
                hits_group[col] = test_hits[!, col]
            end
        end

        try
            # Test successful loading and preparation
            result_df = load_and_prepare_data(temp_dir, "test_data.h5", 100.0, 100.0)

            @test result_df isa DataFrame
            @test Set(names(result_df)) == Set(names(test_hits))

            # Check fiducial cuts were applied (Event 2 should be removed due to |x| > 100)
            unique_events = unique(result_df.event_id)
            @test length(unique_events) == 2  # Events 1 and 3 should pass
            @test 1 in unique_events
            @test 3 in unique_events
            @test !(2 in unique_events)  # Event 2 fails due to x coordinates

            # Test with different fiducial cuts
            result_df2 = load_and_prepare_data(temp_dir, "test_data.h5", 250.0, 100.0)
            unique_events2 = unique(result_df2.event_id)
            @test length(unique_events2) == 3  # All events should pass with larger xyc

            # Test with non-existent file
            @test_throws ArgumentError load_and_prepare_data(temp_dir, "nonexistent.h5", 100.0, 100.0)

            # Test with non-existent directory
            @test_throws ArgumentError load_and_prepare_data("/nonexistent/dir", "test_data.h5", 100.0, 100.0)

        finally
            # Clean up
            rm(temp_dir, recursive=true)
        end
    end
end

# Note: Full integration tests for event_loop and event_loop_single_track
# require actual HDF5 files with proper structure, which would be better
# tested with real data files or mock files that match the expected format.

@testset "Event Loop Integration (Mock)" begin
    @testset "Parameter Integration" begin
        # Test that parameter validation is called
        @test_throws ArgumentError event_loop(".", events_to_run=0)
        @test_throws ArgumentError event_loop(".", voxel_size_mm=0.0)
        @test_throws ArgumentError event_loop(".", max_distance_mm=-1.0)
        @test_throws ArgumentError event_loop(".", energy_threshold_kev=-5.0)

        @test_throws ArgumentError event_loop_single_track(".", events_to_run=0)
        @test_throws ArgumentError event_loop_single_track(".", voxel_size_mm=0.0)
        @test_throws ArgumentError event_loop_single_track(".", max_distance_mm=-1.0)
        @test_throws ArgumentError event_loop_single_track(".", energy_threshold_kev=-5.0)
    end

    @testset "File Handling" begin
        # Test with non-existent directory
        @test_throws ArgumentError event_loop("/nonexistent/directory")
        @test_throws ArgumentError event_loop_single_track("/nonexistent/directory")

        # Test with directory that doesn't contain the expected file
        temp_dir = mktempdir()
        try
            @test_throws ArgumentError event_loop(temp_dir, input_file="nonexistent.h5")
            @test_throws ArgumentError event_loop_single_track(temp_dir, input_file="nonexistent.h5")
        finally
            rm(temp_dir, recursive=true)
        end
    end
end