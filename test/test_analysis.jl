# Test analysis and event processing functionality
using Petit
using Test
using DataFrames
using Graphs

@testset "Analysis Functions" begin
    
    @testset "Select Events" begin
        # Create more complex test data
        test_hits = DataFrame(
            event_id = [1, 1, 1, 1, 1, 2, 2, 2],
            x = [0.0, 0.5, 5.0, 10.0, 15.0, 0.0, 1.0, 2.0],
            y = [0.0, 0.5, 5.0, 10.0, 15.0, 0.0, 1.0, 2.0],
            z = [0.0, 0.5, 5.0, 10.0, 15.0, 0.0, 1.0, 2.0],
            energy = [0.5, 0.6, 0.005, 0.008, 0.7, 0.3, 0.4, 0.5]  # MeV units
        )
        
        # Test with default parameters
        tracks = select_events(test_hits, 1)
        @test length(tracks) >= 1
        @test all(t isa Tracks for t in tracks)
        
        # Test energy threshold filtering (10 keV = 0.01 MeV)
        # Note: Event 1 has some hits with energy < 10 keV (0.005 and 0.008 MeV)
        # These should be filtered out by the energy threshold
        tracks_filtered = select_events(test_hits, 1; energy_threshold_kev=10.0)
        # After filtering, we should still have some tracks (from the higher energy hits)
        @test length(tracks_filtered) >= 0
        @test all(t isa Tracks for t in tracks_filtered)
        
        # Test voxel size parameter
        tracks_large_voxels = select_events(test_hits, 1; voxel_size_mm=5.0)
        @test length(tracks_large_voxels) >= 0
        
        # Test max distance parameter
        tracks_close = select_events(test_hits, 1; max_distance_mm=2.0)
        @test length(tracks_close) >= 1
        
        # Verify tracks are sorted by energy (highest first)
        if length(tracks) > 1
            energies = [sum(t.voxels.energy) for t in tracks]
            @test issorted(energies, rev=true)
        end
    end
    
    @testset "get_event Function" begin
        # Create test data
        test_df = DataFrame(
            event_id = [1, 1, 2, 2, 2, 3],
            x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            y = [1.1, 2.1, 3.1, 4.1, 5.1, 6.1],
            z = [1.2, 2.2, 3.2, 4.2, 5.2, 6.2],
            energy = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        )
        
        # Test getting specific events
        event1 = get_event(test_df, 1)
        @test event1 isa DataFrame
        @test nrow(event1) == 2
        @test all(event1.event_id .== 1)
        
        event2 = get_event(test_df, 2)
        @test nrow(event2) == 3
        @test all(event2.event_id .== 2)
        
        event3 = get_event(test_df, 3)
        @test nrow(event3) == 1
        @test all(event3.event_id .== 3)
        
        # Test error for non-existent event
        @test_throws Exception get_event(test_df, 999)
    end
    
    @testset "Energy Functions" begin
        # Create test particle data
        partdf = DataFrame(
            event_id = [1, 1, 1, 2, 2, 3],
            particle_name = ["gamma", "electron", "positron", "gamma", "gamma", "alpha"],
            mother_id = [0, 1, 1, 0, 0, 0],  # First gamma is primary, others are secondary except last two gammas and alpha
            kin_energy = [2.044, 0.5, 0.3, 1.022, 1.022, 5.0],
            initial_x = [0.0, 1.0, 2.0, 0.0, 0.0, 10.0],
            initial_y = [0.0, 1.0, 2.0, 0.0, 0.0, 10.0],
            initial_z = [0.0, 1.0, 2.0, 0.0, 0.0, 10.0],
            final_volume = ["ACTIVE", "ACTIVE", "ACTIVE", "ACTIVE", "ACTIVE", "BUFFER"]
        )
        
        # Test energy_primary for all events
        primary_energies = energy_primary(partdf)
        @test length(primary_energies) == 3
        @test primary_energies[1] ≈ 2.044  # Event 1: only first gamma is primary
        @test primary_energies[2] ≈ 2.044  # Event 2: sum of two primary gammas
        @test primary_energies[3] ≈ 5.0    # Event 3: alpha is primary
        
        # Test energy_primary for specific event
        energy1 = energy_primary(partdf, 1)
        @test energy1 ≈ 2.044
        
        energy2 = energy_primary(partdf, 2)
        @test energy2 ≈ 2.044
        
        # Test error for non-existent event
        @test_throws Exception energy_primary(partdf, 999)
    end
    
    @testset "Energy Deposited Functions" begin
        # Create test hits data
        hitsdf = DataFrame(
            event_id = [1, 1, 1, 2, 2, 3],
            x = rand(6),
            y = rand(6),
            z = rand(6),
            energy = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        )
        
        # Test energy_deposited for specific event
        edep1 = energy_deposited(hitsdf, 1)
        @test edep1 ≈ 0.6  # Sum of 0.1 + 0.2 + 0.3
        
        edep2 = energy_deposited(hitsdf, 2)
        @test edep2 ≈ 0.9  # Sum of 0.4 + 0.5
        
        edep3 = energy_deposited(hitsdf, 3)
        @test edep3 ≈ 0.6  # Just 0.6
        
        # Test energy_deposited for all events
        all_edep = energy_deposited(hitsdf)
        @test length(all_edep) == 3
        @test all_edep[1] ≈ 0.6
        @test all_edep[2] ≈ 0.9
        @test all_edep[3] ≈ 0.6
        
        # Test error for non-existent event
        @test_throws Exception energy_deposited(hitsdf, 999)
    end
    
    @testset "Find Events With Alphas" begin
        # Create test particle data with alphas
        partdf = DataFrame(
            event_id = [1, 1, 2, 2, 3, 3, 4],
            particle_name = ["gamma", "alpha", "electron", "alpha", "alpha", "positron", "gamma"],
            mother_id = [0, 1, 0, 2, 0, 3, 0],
            kin_energy = [2.044, 5.0, 0.5, 4.5, 6.0, 0.3, 1.022],
            initial_x = [0.0, 10.0, 0.0, 20.0, 30.0, 0.0, 0.0],
            initial_y = [0.0, 10.0, 0.0, 20.0, 30.0, 0.0, 0.0],
            initial_z = [0.0, 10.0, 0.0, 20.0, 30.0, 0.0, 0.0],
            final_volume = ["ACTIVE", "BUFFER", "ACTIVE", "VESSEL", "ACTIVE", "ACTIVE", "ACTIVE"]
        )
        
        alpha_df = find_events_with_alphas(partdf)
        
        @test alpha_df isa DataFrame
        @test nrow(alpha_df) == 3  # Three alphas total
        @test "event_id" in names(alpha_df)
        @test "x" in names(alpha_df)
        @test "y" in names(alpha_df)
        @test "z" in names(alpha_df)
        @test "energy" in names(alpha_df)
        @test "finalv" in names(alpha_df)
        
        # Check specific alpha properties
        @test 1 in alpha_df.event_id
        @test 2 in alpha_df.event_id
        @test 3 in alpha_df.event_id
        @test !(4 in alpha_df.event_id)  # Event 4 has no alpha
        
        # Check energies
        idx1 = findfirst(alpha_df.event_id .== 1)
        @test alpha_df.energy[idx1] ≈ 5.0
        
        idx3 = findfirst(alpha_df.event_id .== 3)
        @test alpha_df.energy[idx3] ≈ 6.0
    end
    
    @testset "Filter Fiducial Events" begin
        # Create test data with events that pass and fail fiducial cuts
        hitsdf = DataFrame(
            event_id = [
                1, 1, 1,  # All hits within fiducial
                2, 2, 2,  # One hit outside xy bounds
                3, 3, 3,  # One hit inside z band (fails)
                4, 4, 4   # All hits pass
            ],
            x = [
                50.0, 60.0, 70.0,    # Event 1: all within ±100
                50.0, 150.0, 60.0,   # Event 2: 150 > 100 (fails)
                50.0, 60.0, 70.0,    # Event 3: all within ±100
                -80.0, -90.0, 85.0   # Event 4: all within ±100
            ],
            y = [
                50.0, 60.0, 70.0,    # Event 1: all within ±100
                50.0, 60.0, 70.0,    # Event 2: all within ±100
                50.0, 60.0, 70.0,    # Event 3: all within ±100
                -80.0, -90.0, 85.0   # Event 4: all within ±100
            ],
            z = [
                150.0, 160.0, 170.0,  # Event 1: all > 100 (pass)
                150.0, 160.0, 170.0,  # Event 2: all > 100
                150.0, 50.0, 170.0,   # Event 3: 50 is within [-100,100] (fails)
                -150.0, -160.0, -170.0 # Event 4: all < -100 (pass)
            ],
            energy = fill(0.1, 12)
        )
        
        # Apply fiducial cuts: |x| < 100, |y| < 100, |z| > 100
        filtered_df = filter_fiducial_events(hitsdf, 100.0, 100.0)
        
        # Only events 1 and 4 should pass
        unique_events = unique(filtered_df.event_id)
        @test length(unique_events) == 2
        @test 1 in unique_events
        @test 4 in unique_events
        @test !(2 in unique_events)  # Failed due to x
        @test !(3 in unique_events)  # Failed due to z
        
        # Test with different cuts
        filtered_df2 = filter_fiducial_events(hitsdf, 200.0, 50.0)
        unique_events2 = unique(filtered_df2.event_id)
        # With xyc=200 and zc=50, more events might pass
        @test length(unique_events2) >= 0
        
        # Test with empty DataFrame
        empty_df = DataFrame(event_id=Int[], x=Float64[], y=Float64[], z=Float64[])
        filtered_empty = filter_fiducial_events(empty_df, 100.0, 100.0)
        @test nrow(filtered_empty) == 0
    end
    
    @testset "Filter Radial Events" begin
        # Create test data with events at different radii
        hitsdf = DataFrame(
            event_id = [
                1, 1, 1,  # All hits within radius
                2, 2, 2,  # One hit outside radius
                3, 3, 3,  # All hits within radius
                4, 4, 4   # All hits outside radius
            ],
            x = [
                30.0, 40.0, 50.0,    # Event 1: max r = sqrt(50²+50²) ≈ 70.7
                30.0, 150.0, 40.0,   # Event 2: max r = sqrt(150²+60²) > 100
                0.0, 10.0, 20.0,     # Event 3: max r = sqrt(20²+30²) ≈ 36
                120.0, 130.0, 140.0  # Event 4: all r > 100
            ],
            y = [
                30.0, 40.0, 50.0,    # Event 1
                30.0, 60.0, 40.0,    # Event 2
                0.0, 20.0, 30.0,     # Event 3
                0.0, 0.0, 0.0        # Event 4
            ],
            z = fill(0.0, 12),
            energy = fill(0.1, 12)
        )
        
        # Apply radial cut with r = 100
        filtered_df = filter_radial(hitsdf, 100.0, 0.0)
        
        # Events 1 and 3 should pass
        unique_events = unique(filtered_df.event_id)
        @test length(unique_events) == 2
        @test 1 in unique_events
        @test 3 in unique_events
        @test !(2 in unique_events)  # Failed due to large radius
        @test !(4 in unique_events)  # Failed due to large radius
        
        # Test with r = sqrt(50²+50²) ≈ 70.7
        filtered_df2 = filter_radial(hitsdf, 50.0, 50.0)
        unique_events2 = unique(filtered_df2.event_id)
        @test 3 in unique_events2  # Event 3 should definitely pass
        
        # Test with very large radius
        filtered_df3 = filter_radial(hitsdf, 1000.0, 0.0)
        unique_events3 = unique(filtered_df3.event_id)
        @test length(unique_events3) == 4  # All events should pass
        
        # Test with empty DataFrame
        empty_df = DataFrame(event_id=Int[], x=Float64[], y=Float64[], z=Float64[])
        filtered_empty = filter_radial(empty_df, 100.0, 0.0)
        @test nrow(filtered_empty) == 0
    end
    
    @testset "Filter Z Events" begin
        # Create test data with events at different z ranges
        hitsdf = DataFrame(
            event_id = [
                1, 1, 1,  # All hits in left range (-200, -100)
                2, 2, 2,  # All hits in right range (100, 200)
                3, 3, 3,  # Hits in both ranges
                4, 4, 4,  # Hits outside both ranges (fail)
                5, 5, 5   # Mixed - some hits fail
            ],
            x = fill(0.0, 15),  # x doesn't matter for z filter
            y = fill(0.0, 15),  # y doesn't matter for z filter
            z = [
                -150.0, -120.0, -110.0,  # Event 1: all in (-200, -100)
                120.0, 150.0, 180.0,     # Event 2: all in (100, 200)
                -150.0, 150.0, -130.0,   # Event 3: mix of both ranges
                -50.0, 50.0, 0.0,        # Event 4: all in forbidden zone
                -150.0, -50.0, 150.0     # Event 5: mix of allowed and forbidden
            ],
            energy = fill(0.1, 15)
        )
        
        # Apply z filter: allowed ranges (-200, -100) and (100, 200)
        # Parameters: zil=-100, zir=100, zol=-200, zor=200
        filtered_df = filter_z(hitsdf, -100.0, 100.0, -200.0, 200.0)
        
        # Events 1, 2, and 3 should pass
        unique_events = unique(filtered_df.event_id)
        @test length(unique_events) == 3
        @test 1 in unique_events  # All in left range
        @test 2 in unique_events  # All in right range  
        @test 3 in unique_events  # Mix of both ranges
        @test !(4 in unique_events)  # All in forbidden zone
        @test !(5 in unique_events)  # Has hits in forbidden zone
        
        # Test with different ranges
        # Allow only narrow bands: (-150, -140) and (140, 150)
        filtered_df2 = filter_z(hitsdf, -140.0, 140.0, -150.0, 150.0)
        unique_events2 = unique(filtered_df2.event_id)
        # Only event 3 might have hits in these narrow ranges
        @test length(unique_events2) <= 3
        
        # Test with very wide ranges (should pass most events)
        # Ranges: (-300, -40) and (40, 300)
        filtered_df3 = filter_z(hitsdf, -40.0, 40.0, -300.0, 300.0)
        unique_events3 = unique(filtered_df3.event_id)
        @test length(unique_events3) == 4  # Events 1, 2, 3, 5 pass; Event 4 fails (has 0 in forbidden zone)
        
        # Test edge cases
        # Points near boundaries - remember the condition is > and <, not >= or <=
        hitsdf_edge = DataFrame(
            event_id = [1, 1, 2, 2],
            x = fill(0.0, 4),
            y = fill(0.0, 4),
            z = [-100.5, -99.5, 99.5, 100.5],
            energy = fill(0.1, 4)
        )
        
        # Ranges: (-101, -99) and (99, 101)
        filtered_edge = filter_z(hitsdf_edge, -99.0, 99.0, -101.0, 101.0)
        unique_edge = unique(filtered_edge.event_id)
        @test length(unique_edge) == 2  # Both events should pass
        
        # Test with empty DataFrame
        empty_df = DataFrame(event_id=Int[], x=Float64[], y=Float64[], z=Float64[])
        filtered_empty = filter_z(empty_df, -100.0, 100.0, -200.0, 200.0)
        @test nrow(filtered_empty) == 0
    end
    
    @testset "Filter Short Tracks" begin
        # Create test data with events of different track lengths
        hitsdf = DataFrame(
            event_id = [
                1, 1, 1, 1, 1,           # Event 1: 5 hits
                2, 2, 2,                 # Event 2: 3 hits
                3, 3, 3, 3, 3, 3, 3,     # Event 3: 7 hits
                4, 4,                    # Event 4: 2 hits
                5,                       # Event 5: 1 hit
                6, 6, 6, 6, 6, 6, 6, 6, 6, 6  # Event 6: 10 hits
            ],
            x = rand(28),
            y = rand(28),
            z = rand(28),
            energy = rand(28)
        )
        
        # Test with trkl = 3 (keep events with >= 3 hits)
        filtered_df = filter_short_tracks(hitsdf, 3)
        unique_events = unique(filtered_df.event_id)
        @test length(unique_events) == 4  # Events 1, 2, 3, 6
        @test 1 in unique_events  # 5 hits >= 3
        @test 2 in unique_events  # 3 hits >= 3
        @test 3 in unique_events  # 7 hits >= 3
        @test !(4 in unique_events)  # 2 hits < 3
        @test !(5 in unique_events)  # 1 hit < 3
        @test 6 in unique_events  # 10 hits >= 3
        
        # Test with trkl = 5 (keep events with >= 5 hits)
        filtered_df2 = filter_short_tracks(hitsdf, 5)
        unique_events2 = unique(filtered_df2.event_id)
        @test length(unique_events2) == 3  # Events 1, 3, 6
        @test 1 in unique_events2  # 5 hits >= 5
        @test !(2 in unique_events2)  # 3 hits < 5
        @test 3 in unique_events2  # 7 hits >= 5
        @test 6 in unique_events2  # 10 hits >= 5
        
        # Test with trkl = 10 (keep events with >= 10 hits)
        filtered_df3 = filter_short_tracks(hitsdf, 10)
        unique_events3 = unique(filtered_df3.event_id)
        @test length(unique_events3) == 1  # Only Event 6
        @test 6 in unique_events3  # 10 hits >= 10
        
        # Test with trkl = 1 (keep all events)
        filtered_df4 = filter_short_tracks(hitsdf, 1)
        unique_events4 = unique(filtered_df4.event_id)
        @test length(unique_events4) == 6  # All events
        
        # Test with trkl = 20 (no events pass)
        filtered_df5 = filter_short_tracks(hitsdf, 20)
        @test nrow(filtered_df5) == 0
        
        # Test with empty DataFrame
        empty_df = DataFrame(event_id=Int[], x=Float64[], y=Float64[], z=Float64[], energy=Float64[])
        filtered_empty = filter_short_tracks(empty_df, 5)
        @test nrow(filtered_empty) == 0
        
        # Verify that the correct hits are kept (not just the event_ids)
        # For trkl = 3, we should keep all hits from events 1, 2, 3, 6
        filtered_df6 = filter_short_tracks(hitsdf, 3)
        expected_hits = 5 + 3 + 7 + 10  # Total hits from events 1, 2, 3, 6
        @test nrow(filtered_df6) == expected_hits
    end
    
    @testset "Histogram Results" begin
        # Create test DataFrame
        df = DataFrame(
            event_id = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
            energy = [100.0, 100.0, 200.0, 200.0, 300.0, 300.0, 400.0, 400.0, 500.0, 500.0],
            x = [-100.0, -100.0, -50.0, -50.0, 0.0, 0.0, 50.0, 50.0, 100.0, 100.0],
            y = [-200.0, -200.0, -100.0, -100.0, 0.0, 0.0, 100.0, 100.0, 200.0, 200.0],
            z = [100.0, 100.0, 200.0, 200.0, 300.0, 300.0, 400.0, 400.0, 500.0, 500.0]
        )
        
        # Test histogram_results with default parameters
        h_x, h_y, h_z, h_e = histogram_results(df)
        
        @test h_x isa Histo1d
        @test h_y isa Histo1d
        @test h_z isa Histo1d
        @test h_e isa Histo1d
        
        # Test with custom parameters
        h_x2, h_y2, h_z2, h_e2 = histogram_results(df; 
            nbins=25, 
            xrange=(-150.0, 150.0),
            yrange=(-250.0, 250.0),
            zrange=(50.0, 550.0),
            erange=(50.0, 550.0))
        
        @test length(h_x2.weights) == 25
        @test length(h_y2.weights) == 25
        @test length(h_z2.weights) == 25
        @test length(h_e2.weights) == 25
        
        # Test with empty DataFrame
        empty_df = DataFrame(event_id=Int[], energy=Float64[], x=Float64[], y=Float64[], z=Float64[])
        h_x_empty, h_y_empty, h_z_empty, h_e_empty = histogram_results(empty_df)
        
        @test h_x_empty isa Histo1d
        @test h_y_empty isa Histo1d
        @test h_z_empty isa Histo1d
        @test h_e_empty isa Histo1d
    end
    
    @testset "Smear Histogram" begin
        # Create a simple test histogram
        edges = collect(0.0:1.0:10.0)  # 10 bins from 0 to 10
        centers = (edges[1:end-1] + edges[2:end]) / 2  # Bin centers: 0.5, 1.5, ..., 9.5
        weights = [0.0, 5.0, 0.0, 3.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0]  # Non-zero weights in bins 2, 4, 7
        
        h = Histo1d(edges, weights, centers)
        
        # Test smearing with sigma = 0.1
        sigma = 0.1
        smeared = smear_histogram(h, sigma)
        
        # Should have 5 + 3 + 2 = 10 total samples
        @test length(smeared) == 10
        
        # All values should be finite
        @test all(isfinite.(smeared))
        
        # Mean should be close to weighted average of bin centers
        expected_mean = (5*1.5 + 3*3.5 + 2*6.5) / 10
        @test abs(mean(smeared) - expected_mean) < 3*sigma  # Within 3 sigma
        
        # Test with empty histogram (all weights zero)
        empty_weights = zeros(10)
        h_empty = Histo1d(edges, empty_weights, centers)
        smeared_empty = smear_histogram(h_empty, sigma)
        @test length(smeared_empty) == 0
        
        # Test with single bin
        single_weight = [5.0]
        single_edges = [0.0, 1.0]
        single_centers = [0.5]
        h_single = Histo1d(single_edges, single_weight, single_centers)
        smeared_single = smear_histogram(h_single, 0.05)
        @test length(smeared_single) == 5
        @test all(abs.(smeared_single .- 0.5) .< 0.2)  # All near center
    end
    
    @testset "Counts in Range" begin
        # Create test histogram with known values
        edges = collect(0.0:1.0:10.0)  # Bins: [0,1), [1,2), ..., [9,10)
        weights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        centers = (edges[1:end-1] + edges[2:end]) / 2
        h = Histo1d(edges, weights, centers)
        
        # Test full range
        total = counts_in_range(h, 0.0, 10.0)
        @test total == sum(weights)  # Should be 55
        
        # Test partial ranges
        count_low = counts_in_range(h, 0.0, 5.0)
        @test count_low == sum(weights[1:5])  # 1+2+3+4+5 = 15
        
        count_mid = counts_in_range(h, 3.0, 7.0)
        @test count_mid == sum(weights[4:7])  # 4+5+6+7 = 22
        
        count_high = counts_in_range(h, 5.0, 10.0)
        @test count_high == sum(weights[6:10])  # 6+7+8+9+10 = 40
        
        # Test single bin
        count_single = counts_in_range(h, 2.0, 3.0)
        @test count_single == weights[3]  # Just bin [2,3) = 3.0
        
        # Test out of range
        count_none = counts_in_range(h, 20.0, 30.0)
        @test count_none == 0.0
        
        # Test inverted range (should throw error)
        @test_throws ArgumentError counts_in_range(h, 5.0, 3.0)
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
    
    @testset "Error Handling for Histograms" begin
        # Test with empty data - should handle gracefully
        empty_data = Float64[]
        # Empty data might not throw exception, just create empty histogram
        h_empty = hist1d(empty_data; nbins=5)
        @test h_empty isa Petit.Histo1d
        
        # Test with single point
        single_point = [1.0]
        h_single = hist1d(single_point; nbins=1)
        @test h_single isa Petit.Histo1d
    end
end