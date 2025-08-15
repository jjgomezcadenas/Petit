# Test I/O functionality including JSON and serialization
using Petit
using Test
using DataFrames
using Serialization
using OrderedCollections

@testset "I/O Functions" begin
    
    @testset "JSON Analysis Summary I/O" begin
        # Test data
        data_type = "znubb"
        trklm = 50
        xr = 1300.0
        zil = -100.0
        zir = 100.0
        zol = -2000.0
        zor = 2000.0
        erex = 12.5
        roi_xi = 2440.0
        roi_xu = 2500.0
        eff_contained = 0.95
        eff_1trk = 0.8
        eff_radial = 0.9
        eff_trkl = 0.85
        eff_zfid = 0.92
        eff_roi = 0.88
        eff_2e = 0.8
        eff_total = 0.45
        
        # Create temporary file for testing
        temp_file = tempname() * ".json"
        
        try
            # Test writing
            json_analysis_summary(temp_file, data_type, trklm, xr,
                                zil, zir, zol, zor,
                                erex, roi_xi, roi_xu,
                                eff_contained, eff_1trk, eff_radial,
                                eff_trkl, eff_zfid, eff_roi,
                                eff_2e, eff_total)
            
            # Check file exists
            @test isfile(temp_file)
            
            # Test reading
            summary = read_json_analysis_summary(temp_file)
            
            # Check that it's an OrderedDict
            @test summary isa OrderedDict
            
            # Check all values are correctly read
            @test summary["data_type"] == data_type
            @test summary["cut_track_length_nhits"] == trklm
            @test summary["cut_radial_mm"] ≈ xr
            @test summary["cut_z_inner_left_mm"] ≈ zil
            @test summary["cut_z_inner_right_mm"] ≈ zir
            @test summary["cut_z_outer_left_mm"] ≈ zol
            @test summary["cut_z_outer_right_mm"] ≈ zor
            @test summary["energy_resolution_keV"] ≈ erex
            @test summary["cut_ROI_left_keV"] ≈ roi_xi
            @test summary["cut_ROI_right_keV"] ≈ roi_xu
            @test summary["eff_contained"] ≈ eff_contained
            @test summary["eff_1trk"] ≈ eff_1trk
            @test summary["eff_radial"] ≈ eff_radial
            @test summary["eff_trkl"] ≈ eff_trkl
            @test summary["eff_zfid"] ≈ eff_zfid
            @test summary["eff_roi"] ≈ eff_roi
            @test summary["eff_2e"] ≈ eff_2e
            @test summary["eff_total"] ≈ eff_total
            
            # Test that field order is preserved
            keys_list = collect(keys(summary))
            expected_order = ["data_type", "energy_resolution_keV", "cut_track_length_nhits",
                            "cut_radial_mm", "cut_z_inner_left_mm", "cut_z_inner_right_mm",
                            "cut_z_outer_left_mm", "cut_z_outer_right_mm",
                            "cut_ROI_left_keV", "cut_ROI_right_keV",
                            "eff_contained", "eff_1trk", "eff_radial", "eff_trkl",
                            "eff_zfid", "eff_roi", "eff_2e", "eff_total"]
            @test keys_list == expected_order
            
        finally
            # Clean up
            rm(temp_file, force=true)
        end
    end
    
    
end