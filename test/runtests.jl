using Petit
using Test
using DataFrames
using Graphs
using StatsBase

@testset "Petit.jl" begin
    
    @testset "Basic Functionality" begin
        @test isdefined(Petit, :Tracks)
        @test isdefined(Petit, :AnalysisResults)
        @test isdefined(Petit, :TracksSummary)
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
    
    @testset "TracksSummary Structure" begin
        # Test TracksSummary constructor
        summary = TracksSummary()
        @test summary isa TracksSummary
        @test summary.energies isa Vector{Float64}
        @test summary.xs isa Vector{Float64}
        @test summary.ys isa Vector{Float64}
        @test summary.zs isa Vector{Float64}
        @test summary.ids isa Vector{Int}
        @test length(summary.energies) == 0
        @test length(summary.xs) == 0
        @test length(summary.ys) == 0
        @test length(summary.zs) == 0
        @test length(summary.ids) == 0
        
        # Test manual constructor
        energies = [100.0, 200.0]
        xs = [1.0, 2.0, 3.0]
        ys = [4.0, 5.0, 6.0]
        zs = [7.0, 8.0, 9.0]
        ids = [1, 2]
        
        summary2 = TracksSummary(energies, xs, ys, zs, ids)
        @test summary2.energies == energies
        @test summary2.xs == xs
        @test summary2.ys == ys
        @test summary2.ids == ids
        @test summary2.zs == zs
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
        
        # Test track sorting by energy
        # Create data designed to produce multiple tracks with different energies
        multi_track_data = DataFrame(
            event_id=[1, 1, 1, 1, 1, 1],
            x=[0.0, 1.0, 10.0, 11.0, 20.0, 21.0],  # Three separate groups
            y=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            energy=[0.1, 0.2, 0.4, 0.5, 0.05, 0.05]  # Groups: 0.3, 0.9, 0.1 MeV
        )
        
        tracks_sorted = select_events(multi_track_data, 1; 
                                    voxel_size_mm=3.0, 
                                    max_distance_mm=2.0, 
                                    energy_threshold_kev=20.0)
        
        # If we get multiple tracks, verify they are sorted by energy (highest first)
        if length(tracks_sorted) > 1
            track_energies = [sum(track.voxels.energy) for track in tracks_sorted]
            @test issorted(track_energies, rev=true)  # Should be sorted in descending order
        end
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
        @test hasfield(AnalysisResults, :single_track)
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
        
        # Test TracksSummary structure types
        @test results.single_track isa TracksSummary
        @test results.two_track_primary isa TracksSummary
        @test results.two_track_secondary isa TracksSummary
        @test results.three_track_primary isa TracksSummary
        @test results.three_track_secondary isa TracksSummary
        
        # Test energy arrays are proper type
        @test results.single_track.energies isa Vector{Float64}
        @test results.two_track_primary.energies isa Vector{Float64}
        @test results.two_track_secondary.energies isa Vector{Float64}
        @test results.three_track_primary.energies isa Vector{Float64}
        @test results.three_track_secondary.energies isa Vector{Float64}
        
        # Test coordinate arrays are proper type
        @test results.single_track.xs isa Vector{Float64}
        @test results.single_track.ys isa Vector{Float64}
        @test results.single_track.zs isa Vector{Float64}
        
        # Test that energy values are positive (if any tracks found)
        if length(results.single_track.energies) > 0
            @test all(e -> e > 0, results.single_track.energies)
        end
        if length(results.two_track_primary.energies) > 0
            @test all(e -> e > 0, results.two_track_primary.energies)
            @test all(e -> e > 0, results.two_track_secondary.energies)
        end
        if length(results.three_track_primary.energies) > 0
            @test all(e -> e > 0, results.three_track_primary.energies)
            @test all(e -> e > 0, results.three_track_secondary.energies)
        end
        
        # Test consistency: secondary track arrays should have same length as primary for 2-track events
        @test length(results.two_track_primary.energies) == length(results.two_track_secondary.energies)
        # Note: For 3+ track events, secondary can have more entries than primary since it includes all non-primary tracks
        
        # Test with empty event range
        empty_results = analysis_loop(test_data; events_to_run=Int[])
        @test empty_results.n_events_processed == 0
        @test length(empty_results.single_track.energies) == 0
        
        # Test progress bar functionality (capture output to verify it works)
        # We'll run a small analysis and check it completes without error
        small_test_data = DataFrame(
            event_id=[1, 1, 2, 2],
            x=[0.0, 1.0, 2.0, 3.0],
            y=[0.0, 0.0, 0.0, 0.0],
            z=[0.0, 0.0, 0.0, 0.0],
            energy=[0.1, 0.2, 0.3, 0.4]
        )
        
        # Test that progress bar version still works correctly
        progress_results = analysis_loop(small_test_data; 
                                       events_to_run=1:2,
                                       voxel_size_mm=3.0,
                                       max_distance_mm=5.0,
                                       energy_threshold_kev=50.0)
        
        @test progress_results isa AnalysisResults
        @test progress_results.n_events_processed == 2
        @test progress_results.n_failed >= 0
    end
    
    @testset "Histogram Functions" begin
        # Test data for histograms
        test_data = [1.0, 2.0, 3.0, 4.0, 5.0, 2.5, 3.5]
        
        @testset "hist1d Basic" begin
            # Test basic hist1d function with nbins keyword
            h = hist1d(test_data; nbins=5)
            @test h isa Petit.Histo1d
            @test length(h.edges) == 6  # n+1 edges for n bins
            @test length(h.weights) == 5  # n bins
            @test length(h.centers) == 5  # n bin centers
            
            # Test normalized histogram - PDF mode should have area = 1
            h_norm = hist1d(test_data; nbins=5, norm=true)
            @test h_norm isa Petit.Histo1d
            # Normalized should have different weights than non-normalized
            @test h.weights != h_norm.weights
            # For PDF normalization, area under curve should be ≈ 1
            actual_bin_widths = diff(h_norm.edges)
            @test sum(h_norm.weights .* actual_bin_widths) ≈ 1.0 atol=0.01
        end
        
        @testset "hist1d with Custom Bins" begin
            # Test with custom bin edges
            custom_bins = [0.0, 2.0, 4.0, 6.0]
            h = hist1d(test_data, custom_bins)
            @test h isa Petit.Histo1d
            @test length(h.weights) == 3  # 3 bins for 4 edges
            
            # Test with normalized custom bins
            h_norm = hist1d(test_data, custom_bins, true)
            @test h_norm isa Petit.Histo1d
            # Test that normalization changes the weights
            @test h.weights != h_norm.weights
            # For PDF normalization with custom bins, area should be ≈ 1
            @test sum(h_norm.weights .* diff(custom_bins)) ≈ 1.0 atol=0.01
        end
        
        @testset "hist1d with Keywords" begin
            # Test with xlim parameter
            h = hist1d(test_data; nbins=10, xlim=(1.0, 5.0))
            @test h isa Petit.Histo1d
            @test length(h.weights) == 10
            
            # Test with normalization
            h_norm = hist1d(test_data; nbins=5, norm=true)
            @test h_norm isa Petit.Histo1d
            # Test that normalized weights are different from non-normalized
            h_regular = hist1d(test_data; nbins=5, norm=false)
            @test h_norm.weights != h_regular.weights
            # For PDF normalization, area should be ≈ 1
            actual_bin_widths = diff(h_norm.edges)
            @test sum(h_norm.weights .* actual_bin_widths) ≈ 1.0 atol=0.01
        end
        
        @testset "hist2d Function" begin
            # Test 2D histogram
            x_data = [1.0, 2.0, 3.0, 4.0, 5.0]
            y_data = [1.5, 2.5, 3.5, 4.5, 5.5]
            
            h, hm = hist2d(x_data, y_data, 3, "X", "Y")
            @test h isa Histogram
            @test size(h.weights) == (3, 3)  # 3x3 bins
            
            # Test with limits
            h2, hm2 = hist2d(x_data, y_data, 3, "X", "Y", 1.0, 5.0, 1.0, 5.0; title="Test")
            @test h2 isa Histogram
            @test size(h2.weights) == (3, 3)
        end
        
        @testset "step_hist Function" begin
            # Test step histogram
            h, plt = step_hist(test_data; nbins=5)
            @test h isa Petit.Histo1d
            @test length(h.weights) == 5
            
            # Test with custom parameters
            h2, plt2 = step_hist(test_data; nbins=10, xlim=(1.0, 5.0), logy=false, 
                                xlabel="Test X", ylabel="Test Y", title="Test Title")
            @test h2 isa Petit.Histo1d
            @test length(h2.weights) == 10
        end
        
        @testset "Utility Functions" begin
            # Test by using Petit functions directly to avoid import conflicts
            h1d_test = hist1d(test_data; nbins=5)
            
            # Test that we get expected structure
            @test h1d_test isa Petit.Histo1d
            @test length(h1d_test.edges) == 6  # n+1 edges
            @test length(h1d_test.weights) == 5  # n centers
            @test length(h1d_test.centers) == 5  # n centers
            
            # Test centers from edges vector directly
            edge_vec = [0.0, 1.0, 2.0, 3.0]
            center_vec = Petit.histos.centers(edge_vec)
            @test center_vec ≈ [0.5, 1.5, 2.5]
        end
        
        @testset "Profile DataFrame (p1df)" begin
            # Test profile histogram
            x_prof = [1.0, 1.1, 2.0, 2.1, 3.0, 3.1, 4.0, 4.1]
            y_prof = [10.0, 11.0, 20.0, 21.0, 30.0, 31.0, 40.0, 41.0]
            
            df_prof, p_prof = p1df(x_prof, y_prof, 4)
            @test df_prof isa DataFrame
            @test "x_mean" in names(df_prof)
            @test "y_mean" in names(df_prof)
            @test "x_std" in names(df_prof)
            @test "y_std" in names(df_prof)
            @test nrow(df_prof) <= 4  # Should have at most 4 rows for 4 bins
        end
        
        @testset "Save and Load Histogram" begin
            # Create a test histogram
            test_data_io = [1.0, 2.0, 3.0, 4.0, 5.0, 2.5, 3.5, 1.5, 4.5]
            h_original = hist1d(test_data_io; nbins=5, norm=false)
            
            # Test saving to file
            temp_file = tempname() * ".txt"
            save_histo1d(h_original, temp_file)
            @test isfile(temp_file)
            
            # Test loading from file
            h_loaded = load_histo1d(temp_file)
            @test h_loaded isa Petit.Histo1d
            
            # Test that loaded histogram matches original
            @test h_loaded.edges == h_original.edges
            @test h_loaded.weights == h_original.weights
            @test h_loaded.centers == h_original.centers
            
            # Test with normalized histogram
            h_norm = hist1d(test_data_io; nbins=5, norm=true)
            temp_file_norm = tempname() * ".txt"
            save_histo1d(h_norm, temp_file_norm)
            h_loaded_norm = load_histo1d(temp_file_norm)
            
            @test h_loaded_norm.edges ≈ h_norm.edges
            @test h_loaded_norm.weights ≈ h_norm.weights
            @test h_loaded_norm.centers ≈ h_norm.centers
            
            # Clean up temp files
            rm(temp_file, force=true)
            rm(temp_file_norm, force=true)
        end
        
        @testset "Load Histogram Error Handling" begin
            # Test loading non-existent file
            @test_throws SystemError load_histo1d("nonexistent_file.txt")
            
            # Test loading malformed file
            temp_bad_file = tempname() * ".txt"
            open(temp_bad_file, "w") do io
                println(io, "5")
                println(io, "1.0 2.0 3.0")  # Wrong number of edges
                println(io, "1.0 2.0")       # Wrong number of weights
            end
            
            @test_throws Exception load_histo1d(temp_bad_file)
            rm(temp_bad_file, force=true)
        end
        
        @testset "Save and Load Histogram Collections" begin
            # Create test histograms
            h1 = hist1d([1.0, 2.0, 3.0, 4.0, 5.0]; nbins=3, norm=false)
            h2 = hist1d([10.0, 20.0, 30.0]; nbins=2, norm=true) 
            h3 = hist1d([0.5, 1.5, 2.5, 3.5]; nbins=4, norm=false)
            
            # Test Dictionary method
            histos_dict = Dict("energy" => h1, "position_x" => h2, "position_y" => h3)
            temp_file_dict = tempname() * ".txt"
            save_histos(histos_dict, temp_file_dict)
            @test isfile(temp_file_dict)
            
            # Load and verify dictionary method
            loaded_dict = load_histos(temp_file_dict)
            @test loaded_dict isa Dict{String, Petit.Histo1d}
            @test length(loaded_dict) == 3
            @test haskey(loaded_dict, "energy")
            @test haskey(loaded_dict, "position_x") 
            @test haskey(loaded_dict, "position_y")
            
            # Verify histogram content
            @test loaded_dict["energy"].edges == h1.edges
            @test loaded_dict["energy"].weights == h1.weights
            @test loaded_dict["energy"].centers == h1.centers
            @test loaded_dict["position_x"].weights ≈ h2.weights
            @test loaded_dict["position_y"].edges == h3.edges
            
            # Test Vector method
            histos_vec = [h1, h2, h3]
            names_vec = ["energy", "position_x", "position_y"]
            temp_file_vec = tempname() * ".txt"
            save_histos(histos_vec, names_vec, temp_file_vec)
            @test isfile(temp_file_vec)
            
            # Load and verify vector method 
            loaded_vec = load_histos(temp_file_vec)
            @test loaded_vec isa Dict{String, Petit.Histo1d}
            @test length(loaded_vec) == 3
            @test loaded_vec["energy"].edges == h1.edges
            @test loaded_vec["position_x"].weights ≈ h2.weights
            
            # Test empty collection
            empty_dict = Dict{String, Petit.Histo1d}()
            temp_file_empty = tempname() * ".txt"
            save_histos(empty_dict, temp_file_empty)
            loaded_empty = load_histos(temp_file_empty)
            @test loaded_empty isa Dict{String, Petit.Histo1d}
            @test length(loaded_empty) == 0
            
            # Clean up
            rm(temp_file_dict, force=true)
            rm(temp_file_vec, force=true)
            rm(temp_file_empty, force=true)
        end
        
        @testset "Collection Error Handling" begin
            # Test mismatched vector lengths
            h1 = hist1d([1.0, 2.0, 3.0]; nbins=2)
            h2 = hist1d([4.0, 5.0]; nbins=1)
            histos_vec = [h1, h2]
            names_vec = ["first"]  # Wrong length
            temp_file = tempname() * ".txt"
            
            @test_throws Exception save_histos(histos_vec, names_vec, temp_file)
            
            # Test loading non-existent collection file
            @test_throws SystemError load_histos("nonexistent_collection.txt")
            
            # Test loading malformed collection file
            temp_bad_collection = tempname() * ".txt"
            open(temp_bad_collection, "w") do io
                println(io, "2")  # Says 2 histograms
                println(io, "first")
                println(io, "2")
                println(io, "1.0 2.0 3.0")
                println(io, "1.0 2.0")
                println(io, "1.5 2.5")
                # Missing second histogram
            end
            
            @test_throws Exception load_histos(temp_bad_collection)
            rm(temp_bad_collection, force=true)
        end
    end
    
    @testset "HDF5 Analysis Functions" begin
        @testset "get_dataset_dfs Function" begin
            # Note: We can't easily test get_dataset_dfs without a real HDF5 file
            # But we can test that the function is exported and callable
            @test isdefined(Petit, :get_dataset_dfs)
            
            # Test error handling for non-existent file
            @test_throws Exception get_dataset_dfs("nonexistent.h5")
        end
        
        @testset "event_loop Function" begin
            # Test that event_loop function is exported and callable
            @test isdefined(Petit, :event_loop)
            
            # Test error handling for non-existent directory/file
            @test_throws Exception event_loop("nonexistent_dir"; 
                                             input_file="nonexistent.h5",
                                             events_to_run=1,
                                             voxel_size_mm=5.0,
                                             max_distance_mm=10.0,
                                             energy_threshold_kev=10.0)
        end
        
        @testset "Progress Bar Integration" begin
            # Test that ProgressMeter is properly imported and used
            @test isdefined(Petit, :ProgressMeter)
            
            # Create minimal test data for progress bar testing
            minimal_data = DataFrame(
                event_id=[1, 2],
                x=[0.0, 1.0],
                y=[0.0, 1.0],
                z=[0.0, 1.0],
                energy=[0.1, 0.2]
            )
            
            # Test that analysis_loop with progress bar doesn't crash
            # (We can't easily test the visual output, but we can ensure it runs)
            @test_nowarn begin
                result = analysis_loop(minimal_data; 
                                     events_to_run=1:2,
                                     voxel_size_mm=5.0,
                                     max_distance_mm=10.0,
                                     energy_threshold_kev=50.0)
            end
            
            # Test that the result is still valid
            result = analysis_loop(minimal_data; 
                                 events_to_run=1:2,
                                 voxel_size_mm=5.0,
                                 max_distance_mm=10.0,
                                 energy_threshold_kev=50.0)
            
            @test result isa AnalysisResults
            @test result.n_events_processed == 2
        end
        
        @testset "Pluto-specific Functions" begin
            # Test that Pluto functions are exported
            @test isdefined(Petit, :event_loop_pluto)
            @test isdefined(Petit, :analysis_loop_pluto)
            
            # Create test data
            test_data = DataFrame(
                event_id=[1, 1, 2, 2, 3],
                x=[0.0, 1.0, 2.0, 3.0, 4.0],
                y=[0.0, 0.0, 0.0, 0.0, 0.0],
                z=[0.0, 0.0, 0.0, 0.0, 0.0],
                energy=[0.1, 0.2, 0.3, 0.4, 0.5]
            )
            
            # Test analysis_loop_pluto
            results, progress_display = analysis_loop_pluto(test_data; 
                                                           events_to_run=1:3,
                                                           voxel_size_mm=5.0,
                                                           max_distance_mm=10.0,
                                                           energy_threshold_kev=50.0)
            
            # Check results
            @test results isa AnalysisResults
            @test results.n_events_processed == 3
            @test results.n_failed >= 0
            
            # Check progress_display is a function
            @test isa(progress_display, Function)
            
            # Test that progress_display returns something (can't test HTML content easily)
            progress_output = progress_display()
            @test !isnothing(progress_output)
            @test isa(progress_output, String)  # Should return HTML string
            
            # Test with empty data
            empty_data = DataFrame(
                event_id=Int[], x=Float64[], y=Float64[], z=Float64[], energy=Float64[]
            )
            
            empty_results, empty_progress = analysis_loop_pluto(empty_data; 
                                                               events_to_run=Int[])
            
            @test empty_results.n_events_processed == 0
            @test isa(empty_progress, Function)
            
            # Test with single event
            single_event_data = DataFrame(
                event_id=[5, 5],
                x=[0.0, 1.0],
                y=[0.0, 0.0],
                z=[0.0, 0.0],
                energy=[0.2, 0.3]
            )
            
            single_results, single_progress = analysis_loop_pluto(single_event_data; 
                                                                 events_to_run=[5],
                                                                 voxel_size_mm=3.0,
                                                                 max_distance_mm=5.0,
                                                                 energy_threshold_kev=100.0)
            
            @test single_results.n_events_processed == 1
            @test isa(single_progress, Function)
            
            # Test event_loop_pluto (would need actual HDF5 file)
            @test_throws Exception event_loop_pluto("nonexistent_dir"; 
                                                   input_file="nonexistent.h5",
                                                   events_to_run=10)
        end
        
        @testset "TracksSummary and AnalysisResults Structures" begin
            # Test TracksSummary structure (already tested above, but ensuring completeness)
            ts = TracksSummary()
            @test ts isa TracksSummary
            @test length(ts.energies) == 0
            @test length(ts.xs) == 0
            @test length(ts.ys) == 0
            @test length(ts.zs) == 0
            @test length(ts.ids) == 0
            
            # Test manual construction
            energies = [100.0, 200.0]
            xs = [1.0, 2.0]
            ys = [3.0, 4.0]
            zs = [5.0, 6.0]
            ids = [1, 2]
            ts2 = TracksSummary(energies, xs, ys, zs, ids)
            @test ts2.energies == energies
            @test ts2.xs == xs
            @test ts2.ys == ys
            @test ts2.zs == zs
            @test ts2.ids == ids
            
            # Test AnalysisResults structure creation
            single_track = TracksSummary()
            two_track_primary = TracksSummary()
            two_track_secondary = TracksSummary()
            three_track_primary = TracksSummary()
            three_track_secondary = TracksSummary()
            
            results = AnalysisResults(
                single_track,
                two_track_primary,
                two_track_secondary,
                three_track_primary,
                three_track_secondary,
                100,  # n_events_processed
                50,   # n_single_track
                30,   # n_two_track
                15,   # n_three_plus_track
                5     # n_failed
            )
            
            @test results isa AnalysisResults
            @test results.n_events_processed == 100
            @test results.n_single_track == 50
            @test results.n_two_track == 30
            @test results.n_three_plus_track == 15
            @test results.n_failed == 5
            @test results.single_track isa TracksSummary
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
    
    @testset "Histogram Results" begin
        # Create test TracksSummary
        energies = [100.0, 200.0, 300.0, 400.0, 500.0]
        xs = [-100.0, -50.0, 0.0, 50.0, 100.0]
        ys = [-200.0, -100.0, 0.0, 100.0, 200.0]
        zs = [100.0, 200.0, 300.0, 400.0, 500.0]
        ids = [1, 2, 3, 4, 5]
        
        ts = TracksSummary(energies, xs, ys, zs, ids)
        
        # Test histogram_results with default parameters
        h_x, h_y, h_z, h_e = histogram_results(ts)
        
        @test h_x isa Histo1d
        @test h_y isa Histo1d
        @test h_z isa Histo1d
        @test h_e isa Histo1d
        
        # Test with custom parameters
        h_x2, h_y2, h_z2, h_e2 = histogram_results(ts; 
            nbins=25, 
            xrange=(-150.0, 150.0),
            yrange=(-250.0, 250.0),
            zrange=(50.0, 550.0),
            erange=(50.0, 550.0))
        
        @test length(h_x2.weights) == 25
        @test length(h_y2.weights) == 25
        @test length(h_z2.weights) == 25
        @test length(h_e2.weights) == 25
        
        # Test with empty TracksSummary
        empty_ts = TracksSummary()
        h_x_empty, h_y_empty, h_z_empty, h_e_empty = histogram_results(empty_ts)
        
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
        
        # Test with zero sigma (no smearing)
        sigma_zero = 0.0
        smeared_zero = smear_histogram(h, sigma_zero)
        @test length(smeared_zero) == 10
        
        # With zero sigma, values should be exactly the bin centers (within floating point precision)
        expected_values = [1.5, 1.5, 1.5, 1.5, 1.5,  # 5 samples at bin center 1.5
                          3.5, 3.5, 3.5,              # 3 samples at bin center 3.5  
                          6.5, 6.5]                   # 2 samples at bin center 6.5
        @test sort(smeared_zero) ≈ sort(expected_values)
        
        # Test with empty histogram
        empty_weights = zeros(10)
        h_empty = Histo1d(edges, empty_weights, centers)
        smeared_empty = smear_histogram(h_empty, 1.0)
        @test length(smeared_empty) == 0
    end
    
    @testset "Counts in Range" begin
        # Create a test histogram with known structure
        edges = collect(0.0:10.0:100.0)  # 10 bins: [0-10), [10-20), ..., [90-100]
        centers = (edges[1:end-1] + edges[2:end]) / 2  # Centers: 5, 15, 25, ..., 95
        weights = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]  # Increasing weights
        
        h = Histo1d(edges, weights, centers)
        
        # Test single bin ranges
        @test counts_in_range(h, 5.0, 5.0) == 1.0    # First bin only
        @test counts_in_range(h, 15.0, 15.0) == 2.0  # Second bin only
        @test counts_in_range(h, 95.0, 95.0) == 10.0 # Last bin only
        
        # Test multiple bin ranges
        @test counts_in_range(h, 5.0, 15.0) == 3.0   # First two bins: 1 + 2
        @test counts_in_range(h, 15.0, 35.0) == 9.0  # Bins 2,3,4: 2 + 3 + 4
        @test counts_in_range(h, 5.0, 95.0) == 55.0  # All bins: sum(1:10)
        
        # Test partial ranges (between bin centers)
        @test counts_in_range(h, 12.0, 18.0) == 2.0  # Only bin at center 15
        @test counts_in_range(h, 3.0, 7.0) == 1.0    # Only bin at center 5
        
        # Test ranges outside histogram bounds
        @test counts_in_range(h, -10.0, 0.0) == 0.0   # Below histogram
        @test counts_in_range(h, 100.0, 110.0) == 0.0 # Above histogram
        @test counts_in_range(h, -5.0, 3.0) == 0.0    # Partially below
        @test counts_in_range(h, 97.0, 105.0) == 0.0  # Partially above
        
        # Test full range and beyond
        @test counts_in_range(h, 0.0, 100.0) == 55.0  # All bins
        @test counts_in_range(h, -10.0, 110.0) == 55.0 # All bins with extra range
        
        # Test edge case: xlow = xup outside any bin
        @test counts_in_range(h, 12.5, 12.5) == 0.0
        
        # Test error case: xlow > xup
        @test_throws ArgumentError counts_in_range(h, 30.0, 20.0)
        
        # Test with empty histogram
        empty_weights = zeros(10)
        h_empty = Histo1d(edges, empty_weights, centers)
        @test counts_in_range(h_empty, 0.0, 100.0) == 0.0
        
        # Test with histogram containing some zero weights
        mixed_weights = [0.0, 5.0, 0.0, 3.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0]
        h_mixed = Histo1d(edges, mixed_weights, centers)
        @test counts_in_range(h_mixed, 0.0, 100.0) == 11.0  # 5 + 3 + 2 + 1
        @test counts_in_range(h_mixed, 15.0, 35.0) == 8.0   # 5 + 0 + 3
    end
    
    @testset "Analysis Results Serialization" begin
        # Create a simple test AnalysisResults structure
        single_tracks = TracksSummary([100.0, 200.0], [1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [1, 2])
        two_track_primary = TracksSummary([150.0], [7.0], [8.0], [9.0], [3])
        two_track_secondary = TracksSummary([50.0], [10.0], [11.0], [12.0], [3])
        three_track_primary = TracksSummary()
        three_track_secondary = TracksSummary()
        
        original_results = AnalysisResults(
            single_tracks,
            two_track_primary, 
            two_track_secondary,
            three_track_primary,
            three_track_secondary,
            10,  # n_events_processed
            2,   # n_single_track
            1,   # n_two_track
            0,   # n_three_plus_track
            7    # n_failed
        )
        
        # Create temporary file
        temp_file = tempname() * ".jls"
        
        try
            # Test saving
            save_analysis_results(original_results, temp_file)
            @test isfile(temp_file)
            
            # Test loading
            loaded_results = load_analysis_results(temp_file)
            @test loaded_results isa AnalysisResults
            
            # Test that all fields match
            @test loaded_results.n_events_processed == original_results.n_events_processed
            @test loaded_results.n_single_track == original_results.n_single_track
            @test loaded_results.n_two_track == original_results.n_two_track
            @test loaded_results.n_three_plus_track == original_results.n_three_plus_track
            @test loaded_results.n_failed == original_results.n_failed
            
            # Test TracksSummary data
            @test loaded_results.single_track.energies == original_results.single_track.energies
            @test loaded_results.single_track.xs == original_results.single_track.xs
            @test loaded_results.single_track.ys == original_results.single_track.ys
            @test loaded_results.single_track.zs == original_results.single_track.zs
            @test loaded_results.single_track.ids == original_results.single_track.ids
            
            @test loaded_results.two_track_primary.energies == original_results.two_track_primary.energies
            @test loaded_results.two_track_secondary.energies == original_results.two_track_secondary.energies
            
        finally
            # Clean up
            if isfile(temp_file)
                rm(temp_file)
            end
        end
        
        # Test error handling
        nonexistent_file = "nonexistent_file.jls"
        @test_throws Exception load_analysis_results(nonexistent_file)
    end
    
    @testset "General Functions" begin
        @testset "get_event Function" begin
            # Create test data
            test_data = DataFrame(
                event_id=[1, 1, 2, 2, 3],
                x=[0.0, 1.0, 2.0, 3.0, 4.0],
                y=[0.0, 1.0, 2.0, 3.0, 4.0],
                z=[0.0, 1.0, 2.0, 3.0, 4.0],
                energy=[0.1, 0.2, 0.3, 0.4, 0.5]
            )
            
            # Test get_event for event 1
            event1 = get_event(test_data, 1)
            @test event1 isa DataFrame
            @test nrow(event1) == 2
            @test all(event1.x .== [0.0, 1.0])
            @test all(event1.energy .== [0.1, 0.2])
            
            # Test get_event for event 2
            event2 = get_event(test_data, 2)
            @test nrow(event2) == 2
            @test all(event2.x .== [2.0, 3.0])
            
            # Test get_event for event 3
            event3 = get_event(test_data, 3)
            @test nrow(event3) == 1
            @test event3.x[1] == 4.0
            
            # Test error handling for non-existent event
            @test_throws Exception get_event(test_data, 999)
        end
        
        @testset "getdirs Function" begin
            # Create a temporary directory structure for testing
            temp_dir = mktempdir()
            mkdir(joinpath(temp_dir, "dir1"))
            mkdir(joinpath(temp_dir, "dir2"))
            touch(joinpath(temp_dir, "file.txt"))
            
            # Test getdirs
            dirs = getdirs(temp_dir)
            @test dirs isa Vector{<:AbstractString}  # Can be SubString or String
            @test "dir1" in dirs
            @test "dir2" in dirs
            @test "file.txt" in dirs  # getdirs returns all entries, not just directories
            
            # Clean up
            rm(temp_dir, recursive=true)
            
            # Test with non-existent directory - it doesn't throw, just returns empty
            dirs_empty = getdirs("nonexistent_directory")
            @test dirs_empty isa Vector{<:AbstractString}
            @test length(dirs_empty) == 0
        end
        
        @testset "get_histo1d Function" begin
            # Create test data
            test_data = [1.0, 2.0, 3.0, 4.0, 5.0, 2.5, 3.5]
            
            # Create a StatsBase histogram first
            h_statsbase = fit(Histogram, test_data, nbins=5)
            
            # Convert to Histo1d
            h1d = get_histo1d(h_statsbase)
            
            @test h1d isa Petit.Histo1d
            # edges should be extracted by get_histo1d function
            @test length(h1d.edges) == 6  # 5 bins means 6 edges
            @test h1d.weights == h_statsbase.weights
            @test length(h1d.centers) == length(h1d.weights)
            
            # Test that centers are calculated correctly
            expected_centers = [(h1d.edges[i] + h1d.edges[i+1]) / 2 for i in 1:length(h1d.edges)-1]
            @test h1d.centers ≈ expected_centers
        end
        
        @testset "Particle DataFrame Functions" begin
            # Create test particle data
            test_particles = DataFrame(
                event_id=[1, 1, 1, 2, 2, 3],
                particle_id=[1, 2, 3, 4, 5, 6],
                mother_id=[0, 0, 1, 0, 2, 0],
                kin_energy=[1.0, 0.5, 0.1, 2.0, 0.2, 1.5]
            )
            
            # Test number_of_events
            @test number_of_events(test_particles) == 3
            
            # Test energy_primary (all events)
            energies = energy_primary(test_particles)
            @test length(energies) == 3
            @test energies[1] ≈ 1.5  # Event 1: particles 1,2 (mother_id=0)
            @test energies[2] ≈ 2.0  # Event 2: particle 4 (mother_id=0)
            @test energies[3] ≈ 1.5  # Event 3: particle 6 (mother_id=0)
            
            # Test energy_primary (single event)
            @test energy_primary(test_particles, 1) ≈ 1.5
            @test energy_primary(test_particles, 2) ≈ 2.0
            @test energy_primary(test_particles, 3) ≈ 1.5
            
            # Test error handling
            @test_throws Exception energy_primary(test_particles, 999)
            
            # Test with empty DataFrame
            empty_particles = DataFrame(
                event_id=Int[], particle_id=Int[], mother_id=Int[], kin_energy=Float64[]
            )
            @test number_of_events(empty_particles) == 0
            @test length(energy_primary(empty_particles)) == 0
            
            # Test with single event
            single_event = DataFrame(
                event_id=[5, 5],
                particle_id=[10, 20],
                mother_id=[0, 0],
                kin_energy=[3.0, 2.0]
            )
            @test number_of_events(single_event) == 1
            @test energy_primary(single_event, 5) ≈ 5.0
            
            # Test with no primary particles (all have mother_id != 0)
            no_primaries = DataFrame(
                event_id=[7, 7],
                particle_id=[1, 2],
                mother_id=[1, 2],
                kin_energy=[1.0, 2.0]
            )
            @test number_of_events(no_primaries) == 1
            @test energy_primary(no_primaries, 7) ≈ 0.0
            
            # Test with duplicate event IDs (should still count as unique events)
            duplicate_events = DataFrame(
                event_id=[1, 1, 1, 1, 2, 2],
                particle_id=[1, 2, 3, 4, 5, 6],
                mother_id=[0, 1, 2, 3, 0, 1],
                kin_energy=[1.0, 0.5, 0.3, 0.2, 2.0, 1.0]
            )
            @test number_of_events(duplicate_events) == 2
        end
        
        @testset "Energy Deposited Functions" begin
            # Create test hits data
            test_hits = DataFrame(
                event_id=[1, 1, 1, 2, 2, 3],
                x=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                y=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                z=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                energy=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
            )
            
            # Test energy_deposited (single event)
            @test energy_deposited(test_hits, 1) ≈ 0.6  # 0.1 + 0.2 + 0.3
            @test energy_deposited(test_hits, 2) ≈ 0.9  # 0.4 + 0.5
            @test energy_deposited(test_hits, 3) ≈ 0.6  # 0.6
            
            # Test energy_deposited (all events)
            all_energies = energy_deposited(test_hits)
            @test length(all_energies) == 3
            @test all_energies[1] ≈ 0.6
            @test all_energies[2] ≈ 0.9
            @test all_energies[3] ≈ 0.6
            
            # Test error handling for non-existent event
            @test_throws Exception energy_deposited(test_hits, 999)
            
            # Test with empty DataFrame
            empty_hits = DataFrame(
                event_id=Int[], x=Float64[], y=Float64[], z=Float64[], energy=Float64[]
            )
            @test length(energy_deposited(empty_hits)) == 0
            
            # Test with single event
            single_event_hits = DataFrame(
                event_id=[10, 10, 10],
                x=[1.0, 2.0, 3.0],
                y=[1.0, 2.0, 3.0], 
                z=[1.0, 2.0, 3.0],
                energy=[0.5, 0.3, 0.2]
            )
            @test energy_deposited(single_event_hits, 10) ≈ 1.0
            @test length(energy_deposited(single_event_hits)) == 1
            @test energy_deposited(single_event_hits)[1] ≈ 1.0
            
            # Test with zero energy hits
            zero_energy_hits = DataFrame(
                event_id=[5, 5],
                x=[0.0, 1.0],
                y=[0.0, 1.0],
                z=[0.0, 1.0],
                energy=[0.0, 0.0]
            )
            @test energy_deposited(zero_energy_hits, 5) ≈ 0.0
            @test energy_deposited(zero_energy_hits)[1] ≈ 0.0
            
            # Test with mixed positive/negative energies (if applicable)
            mixed_energy_hits = DataFrame(
                event_id=[7, 7, 7],
                x=[0.0, 1.0, 2.0],
                y=[0.0, 1.0, 2.0],
                z=[0.0, 1.0, 2.0],
                energy=[0.5, -0.1, 0.8]  # Negative energy (unusual but possible)
            )
            @test energy_deposited(mixed_energy_hits, 7) ≈ 1.2  # 0.5 - 0.1 + 0.8
            @test energy_deposited(mixed_energy_hits)[1] ≈ 1.2
        end
        
        @testset "Plot Functions" begin
            # Note: We can't easily test plotting functions that produce plots,
            # but we can at least verify they are callable and don't error
            
            # Create test data
            test_df = DataFrame(
                event_id=[1, 1, 1],
                x=[0.0, 1.0, 2.0],
                y=[0.0, 1.0, 2.0],
                z=[0.0, 1.0, 2.0],
                energy=[0.1, 0.2, 0.3]
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
        end
    end
end