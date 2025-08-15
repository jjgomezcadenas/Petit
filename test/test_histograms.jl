# Test histogram functionality
using Petit
using Test
using DataFrames
using StatsBase

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