### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 349825ff-7ffe-4fa1-ba26-a772041f0323
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Plots
	using Printf
	using HDF5
	using Markdown
	using InteractiveUtils
	using Statistics
	using StatsBase
	using Distributions
	using StatsPlots
	
end

# ╔═╡ c9fc0547-0e73-4629-9909-e59c3d75169d
begin
	using Petit
end

# ╔═╡ 04b446d6-f34f-11ed-2565-0b15d65b6781
PlutoUI.TableOfContents(title="HD5t Post-Analysis", indent=true)


# ╔═╡ 871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
begin
	cmdir=joinpath(ENV["DATA"], "HD5t")
	pdir =joinpath(ENV["PROJECTS"], "Petit")
end

# ╔═╡ 947c237c-9852-40e9-a83f-c23666db90aa
begin
      using Pkg
      Pkg.activate(pdir)
      Pkg.instantiate()
  end

# ╔═╡ 3ae5d298-b9df-4f4e-85d9-d75e4aee1654
md"""
  | File | Particles | Hits |
  |:-----|----------:|-----:|
  | 0nubb.h5 | 8,995 | 8,995 |
  | bi214_copper.h5 | 2,977 | 2,977 |
  | bi214_surface.h5 | 11,150 | 11,150 |
  | electron_0nubb_energy.h5 | 8,483 | 8,483 |
  | tl208_copper.h5 | 3,514 | 3,514 |
  | tl208_surface.h5 | 5,121 | 5,121 |
  """

# ╔═╡ 6c59aeae-7990-4b43-8378-0de210a3291a
md"""
# Load Analysis Results
This notebook loads and analyzes the saved results from a previous analysis run.
"""

# ╔═╡ ae89f5dc-e958-496a-91ac-0bd977355563
md"""
## Select Results File
"""

# ╔═╡ 4af9a3ef-e883-4bc3-a2f1-212102e4951b
begin
	# Default results file location
	results_dir = joinpath(pdir, "znubb")
	results_znubb = joinpath(results_dir, "0nubb.jls")
	
	md"""
	Loading results from: $(results_znubb)
	"""
end

# ╔═╡ 5caa7c21-82e5-420b-b4e7-a0e33543b74b
md"Load znubb results? $(@bind load_znubb CheckBox(default=true))"

# ╔═╡ 8a55b4a3-5cbf-48c3-b150-2bd4ad73f440
if load_znubb
	rznubb = Petit.load_analysis_results(results_znubb)
	md"""
	✅ **znubb Results!**
	
	- **Events processed:** $(rznubb.n_events_processed)
	- **Single track events:** $(rznubb.n_single_track)
	- **Two track events:** $(rznubb.n_two_track)
	- **Three+ track events:** $(rznubb.n_three_plus_track)
	- **Failed events:** $(rznubb.n_failed)
	"""
else
	md"""
	⚠️ Check the box above to load the results.
	"""
end

# ╔═╡ 2b3c4d5e-f6a7-4b8c-9d1e-2f3a4b5c6d7e
md"""
## Select Background File
Choose which background analysis results to load:
"""

# ╔═╡ 9e2f3d4e-a5b6-4c7d-8e9f-1a2b3c4d5f6a
@bind background_file Select([
	"bi214_copper" => "Bi214 Copper",
	"bi214_surface" => "Bi214 Surface",
	"tl208_copper" => "Tl208 Copper",
	"tl208_surface" => "Tl208 Surface",
	"electron_0nubb_energy" => "Electron 0νββ Energy"
], default="bi214_copper")

# ╔═╡ 9f3e1cdb-758c-4223-ac86-1d63ac8f03da
begin
	# Construct the results file path based on selection
	results_bkg = joinpath(results_dir, "$(background_file).jls")
	# Extract display name from the file name
	bkg_display_name = replace(background_file, "_" => " ") |> titlecase
	md"""
	Load $(bkg_display_name) results? $(@bind load_bkg CheckBox(default=true))
	
	Results file: $(results_bkg)
	"""
end

# ╔═╡ 4d5e6f7a-b8c9-4d1e-2f3a-4b5c6d7e8f9a
if load_bkg
	if isfile(results_bkg)
		rbi214 = Petit.load_analysis_results(results_bkg)
		md"""
		✅ **$(bkg_display_name) Results Loaded!**
		
		- **Events processed:** $(rbi214.n_events_processed)
		- **Single track events:** $(rbi214.n_single_track)
		- **Two track events:** $(rbi214.n_two_track)
		- **Three+ track events:** $(rbi214.n_three_plus_track)
		- **Failed events:** $(rbi214.n_failed)
		"""
	else
		md"""
		⚠️ Results file not found: $(results_bkg)
		
		Please ensure the analysis has been run for $(bkg_display_name).
		"""
	end
else
	md"""
	⚠️ Check the box above to load the $(bkg_display_name) results.
	"""
end

# ╔═╡ 07b4e4c1-0469-41ca-9890-bb4f990b4645
md"""
## Analysis Summary
"""

# ╔═╡ f60557f4-113e-44ab-ab53-56e967f8fda8
if load_znubb
	md"""
	### Event Classification Summary (0νββ)
	
	| Event Type | Count | Percentage |
	|:-----------|------:|----------:|
	| Single Track | $(rznubb.n_single_track) | $(round(100*rznubb.n_single_track/rznubb.n_events_processed, digits=1))% |
	| Two Tracks | $(rznubb.n_two_track) | $(round(100*rznubb.n_two_track/rznubb.n_events_processed, digits=1))% |
	| Three+ Tracks | $(rznubb.n_three_plus_track) | $(round(100*rznubb.n_three_plus_track/rznubb.n_events_processed, digits=1))% |
	| Failed | $(rznubb.n_failed) | $(round(100*rznubb.n_failed/rznubb.n_events_processed, digits=1))% |
	| **Total** | **$(rznubb.n_events_processed)** | **100.0%** |
	"""
end

# ╔═╡ 5e6f7a8b-c9d1-4e2f-3a4b-5c6d7e8f9a1b
if load_bkg
	md"""
	### Event Classification Summary (Bi214)
	
	| Event Type | Count | Percentage |
	|:-----------|------:|----------:|
	| Single Track | $(rbi214.n_single_track) | $(round(100*rbi214.n_single_track/rbi214.n_events_processed, digits=1))% |
	| Two Tracks | $(rbi214.n_two_track) | $(round(100*rbi214.n_two_track/rbi214.n_events_processed, digits=1))% |
	| Three+ Tracks | $(rbi214.n_three_plus_track) | $(round(100*rbi214.n_three_plus_track/rbi214.n_events_processed, digits=1))% |
	| Failed | $(rbi214.n_failed) | $(round(100*rbi214.n_failed/rbi214.n_events_processed, digits=1))% |
	| **Total** | **$(rbi214.n_events_processed)** | **100.0%** |
	"""
end

# ╔═╡ 144d60a3-d70f-442b-a252-76178fecdbf7
md"""
## Track Energy Distributions
"""

# ╔═╡ 0fa67f2a-50a5-4c82-b42b-07f45f14e914
md"""
### Single Track Events (0νββ): e = $(round(mean(rznubb.single_track.energies), digits=2))
"""

# ╔═╡ eba21011-2c2a-4540-b096-d39141d9abb2
0.5 * 2457.5/100

# ╔═╡ 6624fb61-e8ac-41e5-a602-c8765e21cede
if load_znubb && length(rznubb.single_track.energies) > 0
	eht1, pht1 = Petit.step_hist(rznubb.single_track.energies;
         nbins = 40,
         xlim   = (2300.0, 2700.0),
         xlabel = "E (keV)",
         ylabel = "Frequency",
         title="Single Track Energy Distribution")
	plot(pht1)
else
	md"No single track events found."
end

# ╔═╡ 6f7a8b9c-d1e2-4f3a-4b5c-6d7e8f9a1b2c
md"""
### Single Track Events (Bi214) : e = $(round(mean(rbi214.single_track.energies), digits=2))
"""

# ╔═╡ 7a8b9c1d-e2f3-4a5b-6c7d-8e9f1a2b3c4d
if load_bkg && length(rbi214.single_track.energies) > 0
	eht1_bi214, pht1_bi214 = Petit.step_hist(rbi214.single_track.energies;
         nbins = 40,
         xlim   = (2300.0, 2700.0),
         xlabel = "E (keV)",
         ylabel = "Frequency",
         title="Single Track Energy Distribution (Bi214)")
	plot(pht1_bi214)
else
	md"No single track events found for Bi214."
end

# ╔═╡ 474b25f8-bd95-4389-867a-bb753dc77d45
md"""
### Energy Resolution Analysis
"""

# ╔═╡ 6c58a746-da45-4e39-a178-91e060b2f34b
md"Energy resolution (keV): $(@bind erex NumberField(0.0:30.0, default=12.5))"

# ╔═╡ aa71c155-9bed-4ad6-963e-575100094a9c
md"ROI low (keV): $(@bind roi_low NumberField(2400.0:2500.0, default=2400.0))"

# ╔═╡ 6c148761-f9e2-4c7d-8d50-ee78bc0a8baf
md"ROI up (keV): $(@bind roi_up NumberField(2475.0:2550.0, default=2500.0))"

# ╔═╡ d1c1a8cb-a70e-480b-8f65-2cfa1a7022c3
if load_znubb && length(rznubb.single_track.energies) > 0
	eres = Petit.smear_histogram(eht1, erex)
	ehrx, phrx = Petit.step_hist(eres;
         nbins = 40,
         xlim   = (2300.0, 2700.0),
         xlabel = "E (keV)",
         ylabel = "Frequency",
         title="bb0nu Energy (with $(erex) keV resolution)")
	plot(phrx)
else
	md"No data to apply resolution."
end

# ╔═╡ 272c88a5-c6c0-41e1-b782-f8e23007eefc
if load_znubb && length(rznubb.single_track.energies) > 0
	counts_roi = Petit.counts_in_range(ehrx, roi_low, roi_up)
	md"""
	**Counts in ROI [$(roi_low), $(roi_up)] keV:** $(counts_roi)
	
	**Efficiency:** $(round(100*counts_roi/rznubb.n_single_track, digits=2))%
	"""
else
	md"No data for ROI analysis."
end

# ╔═╡ 7fde840d-0c00-4593-b0ae-52d1990e8c4e
if load_bkg && length(rbi214.single_track.energies) > 0
	eres_bi214 = Petit.smear_histogram(eht1_bi214, erex)
	ehrx_bi214, phrx_bi214 = Petit.step_hist(eres_bi214;
         nbins = 40,
         xlim   = (2300.0, 2700.0),
         xlabel = "E (keV)",
         ylabel = "Frequency",
         title="Bi-214 Energy (with $(erex) keV resolution)")
	plot(phrx_bi214)
else
	md"No data to apply resolution."
end

# ╔═╡ 7c20e495-0c39-44a5-9c40-8cd9c5d8e3de
md"""
### Signal Efficiency vs ROI
"""

# ╔═╡ cd3ef36e-7150-47a5-9cec-9fa4d0e5779e
 

# ╔═╡ 37d7c197-0a10-4145-ab85-b5a22eae273e
md"""
## Multi-Track Events Analysis
"""

# ╔═╡ 8c929794-64a4-48eb-bd85-ccfad6133b1e
md"""
## Track Position Distributions
"""

# ╔═╡ c6bdf2b5-0b45-4e6a-af2a-b8f1c1ed4b78
md"""
## Statistical Summary
"""

# ╔═╡ ee47c998-1649-4df2-bdc2-33bf53818e62
if load_znubb
	single_energies = rznubb.single_track.energies
	two_primary_energies = rznubb.two_track_primary.energies
	
	stats_df = DataFrame(
		Statistic = ["Mean", "Median", "Std Dev", "Min", "Max", "Count"],
		Single_Track = if length(single_energies) > 0
			[
				round(mean(single_energies), digits=1),
				round(median(single_energies), digits=1),
				round(std(single_energies), digits=1),
				round(minimum(single_energies), digits=1),
				round(maximum(single_energies), digits=1),
				length(single_energies)
			]
		else
			fill("N/A", 6)
		end,
		Two_Track_Primary = if length(two_primary_energies) > 0
			[
				round(mean(two_primary_energies), digits=1),
				round(median(two_primary_energies), digits=1),
				round(std(two_primary_energies), digits=1),
				round(minimum(two_primary_energies), digits=1),
				round(maximum(two_primary_energies), digits=1),
				length(two_primary_energies)
			]
		else
			fill("N/A", 6)
		end
	)
	
	md"""
	### Energy Statistics (keV)
	
	$(stats_df)
	"""
else
	md"Load results to see statistics."
end

# ╔═╡ 1d2e3f4a-b5c6-7d8e-9f1a-2b3c4d5e6f7a
if load_bkg
	single_energies_bi214 = rbi214.single_track.energies
	two_primary_energies_bi214 = rbi214.two_track_primary.energies
	
	stats_df_bi214 = DataFrame(
		Statistic = ["Mean", "Median", "Std Dev", "Min", "Max", "Count"],
		Single_Track = if length(single_energies_bi214) > 0
			[
				round(mean(single_energies_bi214), digits=1),
				round(median(single_energies_bi214), digits=1),
				round(std(single_energies_bi214), digits=1),
				round(minimum(single_energies_bi214), digits=1),
				round(maximum(single_energies_bi214), digits=1),
				length(single_energies_bi214)
			]
		else
			fill("N/A", 6)
		end,
		Two_Track_Primary = if length(two_primary_energies_bi214) > 0
			[
				round(mean(two_primary_energies_bi214), digits=1),
				round(median(two_primary_energies_bi214), digits=1),
				round(std(two_primary_energies_bi214), digits=1),
				round(minimum(two_primary_energies_bi214), digits=1),
				round(maximum(two_primary_energies_bi214), digits=1),
				length(two_primary_energies_bi214)
			]
		else
			fill("N/A", 6)
		end
	)
	
	md"""
	### Energy Statistics (keV) - $(bkg_display_name)
	
	$(stats_df_bi214)
	"""
else
	md"Load $(bkg_display_name) results to see statistics."
end

# ╔═╡ 21f70e07-ef56-4b84-87df-dd4e0d468ed1
md"""
## Export Options
"""

# ╔═╡ 9a89ffdc-3889-48c0-a952-ad2c36094b4c
md"Export histograms to CSV? $(@bind export_csv CheckBox(default=false))"

# ╔═╡ a60ea1b2-1c4c-4565-ac34-d2580dd016e3
if export_csv && load_znubb
	# Export single track histogram
	if length(rznubb.single_track.energies) > 0
		hist_df = DataFrame(
			energy_keV = ehrx.edges[1:end-1],
			counts = ehrx.weights
		)
		csv_file = joinpath(results_dir, "single_track_histogram.csv")
		CSV.write(csv_file, hist_df)
		md"""
		✅ Exported histogram to: `$(csv_file)`
		"""
	else
		md"No data to export."
	end
else
	md"Check the box above to export histogram data."
end

# ╔═╡ 995555f1-dd99-4903-a3b7-9e48b62d675d
md"""
## Debug Information
"""

# ╔═╡ 1362e7bf-0555-4551-9fcd-6047b7d9b555
md"Show debug info? $(@bind show_debug CheckBox(default=false))"

# ╔═╡ 37110c55-a96e-4c4e-9ca1-9fb464865af0
if show_debug && load_znubb
	md"""
	### Results Structure Fields
	
	- `single_track`: $(fieldnames(typeof(rznubb.single_track)))
	- `two_track_primary`: $(fieldnames(typeof(rznubb.two_track_primary)))
	- Number of single track entries: $(length(rznubb.single_track.energies))
	- Number of two track primary entries: $(length(rznubb.two_track_primary.energies))
	- Number of three track secondary entries: $(length(rznubb.three_track_secondary.energies))
	"""
else
	md"Check the box above to see debug information."
end

# ╔═╡ d115f04b-1cf2-47b4-81fd-d15aa8d6cd97
md"""
---
## Summary

This notebook successfully loads and analyzes the pre-computed results from the batch analysis. The `results` object contains all the track information and statistics from the analyzed events.

Key features:
- ✅ Loads saved analysis results from `.jls` files
- ✅ Displays event classification statistics
- ✅ Shows energy distributions for different track types
- ✅ Applies energy resolution smearing
- ✅ Calculates ROI efficiencies
- ✅ Visualizes position distributions
- ✅ Provides statistical summaries
- ✅ Optional CSV export functionality
"""

# ╔═╡ d0f9e26b-2a3c-42a4-a826-7c8694f5d470
md"""
## Functions
"""

# ╔═╡ 1e6adccf-b640-4fd7-92cc-02cf61f27a24
function signal_eff(ehrx, rlow, rup; step=10.0) 
	countx = []
	norm = sum(ehrx.weights)
	for rx in rlow:step:rup
		push!(countx, Petit.counts_in_range(ehrx, rx, rup))
	end
	countx = countx /norm
end

# ╔═╡ 9c7ca8ed-5d3f-436f-bd82-675a03e6a5e9
if load_znubb && length(rznubb.single_track.energies) > 0
	step = 5.0
	xx = roi_low:step:roi_up
	countx = signal_eff(ehrx, roi_low, roi_up; step=step)
	plot(xx, countx,
		xlabel="ROI Lower Bound (keV)",
		ylabel="Signal Efficiency",
		title="Signal Efficiency vs ROI Lower Bound",
		marker=:circle,
		legend=false)
else
	md"No data for efficiency plot."
end

# ╔═╡ 1963d9dc-0cb9-4b8e-98ee-41491c6c784b
if load_bkg && length(rbi214.single_track.energies) > 0
	#step = 5.0
	#xx = roi_low:step:roi_up
	countx_bi214 = signal_eff(ehrx_bi214, roi_low, roi_up; step=step)
	plot(xx, countx_bi214,
		xlabel="ROI Lower Bound (keV)",
		ylabel="Signal Efficiency",
		title="Signal Efficiency vs ROI Lower Bound",
		marker=:circle,
		legend=false)
else
	md"No data for efficiency plot."
end

# ╔═╡ 705979ca-69e9-48e1-9440-ef0ab8c66154
function integrated_signal_eff(ehrx, rlow, rup, rmin, rmax) 
	n1 = counts_in_range(ehrx, rlow, rup)
	n2 = counts_in_range(ehrx, rmin, rmax)
	return n1, n2, n2/n1
end

# ╔═╡ 742b17ab-27fc-418c-9265-522fde5235bb
integrated_signal_eff(ehrx, 2400.0, 2500.0, 2430.0, 2500.0)

# ╔═╡ b3634a2f-6600-4480-b316-2a79af18653e
integrated_signal_eff(ehrx_bi214, 2400.0, 2500.0, 2430.0, 2500.0)

# ╔═╡ f301f5cf-0dc6-49d9-bd4d-ff117ea56a2e
begin
	signal_2e = 0.8
	signal_contained = 0.9
	signal_1trk = rznubb.n_single_track/rznubb.n_events_processed
	signal_roi = integrated_signal_eff(ehrx, 2400.0, 2500.0, 2430.0, 2500.0)[3]
	signal_teff = signal_2e * signal_contained * signal_1trk * signal_roi
	md"""
	### Signal efficiency ($(erex) keV FWHM resolution)
	- Fraction of events contained $(signal_contained)
	- Fraction of events 1 track $(round(signal_1trk, digits=2))
	- Fraction of events ROIS $(round(signal_1trk, digits=2))
	- Fraction of events 1/2 electron $(signal_2e)
	- Total signal efficiency: $(round(signal_teff, digits=2))
	"""
end

# ╔═╡ 06bac1d3-c340-49ce-bd7c-92164e412590
begin
	bi214_cu_2e = 0.04
	bi214_cu_contained = 3e-6
	bi214_cu_1trk = rbi214.n_single_track/rbi214.n_events_processed
	bi214_cu_roi = integrated_signal_eff(ehrx_bi214, 2400.0, 2500.0, 2430.0, 2500.0)[3]
	bi214_cu_teff = bi214_cu_2e * bi214_cu_contained * bi214_cu_1trk * bi214_cu_roi
	md"""
	### Signal efficiency ($(erex) keV FWHM resolution)
	- Fraction of events contained $(bi214_cu_contained)
	- Fraction of events 1 track $(round(bi214_cu_1trk, digits=2))
	- Fraction of events ROIS $(round(bi214_cu_1trk, digits=2))
	- Fraction of events 1/2 electron $(bi214_cu_2e)
	- Total signal efficiency: $(bi214_cu_teff)
	"""
end

# ╔═╡ ea1afc2f-9467-4b8b-a3e3-237925ac9fa1
function histogram_energies_trks(rznubb)
	
	p1 = if length(rznubb.single_track.energies) > 0
		_, p = Petit.step_hist(rznubb.single_track.energies;
	         nbins = 50,
	         xlim   = (0.0, 2700.0),
	         xlabel = "E (keV)",
	         ylabel = "Frequency",
	         title="Single track")
		p
	else
		plot(title="No single track data")
	end
	
	p2 = if length(rznubb.two_track_primary.energies) > 0
		_, p = Petit.step_hist(rznubb.two_track_primary.energies;
	         nbins = 50,
	         xlim   = (0.0, 2700.0),
	         xlabel = "E (keV)",
	         ylabel = "Frequency",
	         title="Two tracks (primary)")
		p
	else
		plot(title="No two track data")
	end
	
	p3 = if length(rznubb.two_track_secondary.energies) > 0
		_, p = Petit.step_hist(rznubb.two_track_secondary.energies;
	         nbins = 50,
	         xlim   = (0.0, 250.0),
	         xlabel = "E (keV)",
	         ylabel = "Frequency",
	         title="Two tracks (secondary)")
		p
	else
		plot(title="No secondary data")
	end
	
	p4 = if length(rznubb.three_track_secondary.energies) > 0
		_, p = Petit.step_hist(rznubb.three_track_secondary.energies;
	         nbins = 50,
	         xlim   = (0.0, 250.0),
	         xlabel = "E (keV)",
	         ylabel = "Frequency",
	         title="Three+ tracks (secondaries)")
		p
	else
		plot(title="No 3+ track data")
	end
	
	plot(p1, p2, p3, p4, layout=(2,2), size=(800,600))
end

# ╔═╡ 5248add6-0d0c-446d-8490-0ccc009fd6e9
if load_znubb
	histogram_energies_trks(rznubb)
else
	md"Load results to see multi-track energy distributions."
end

# ╔═╡ 8b9c1d2e-f3a4-5b6c-7d8e-9f1a2b3c4d5e
if load_bkg
	histogram_energies_trks(rbi214)
else
	md"Load $(bkg_display_name) results to see multi-track energy distributions."
end

# ╔═╡ de61fd85-48cf-4c49-8abc-f929bdc4c3fc
function plot_position_distributions(rznubb)
	# Single track positions
	p1 = if length(rznubb.single_track.xs) > 0
		scatter(rznubb.single_track.xs, rznubb.single_track.ys,
			xlabel="X (mm)", ylabel="Y (mm)",
			title="Single Track XY",
			markersize=2, alpha=0.5, legend=false)
	else
		plot(title="No single track data")
	end
	
	p2 = if length(rznubb.single_track.zs) > 0
		histogram(rznubb.single_track.zs,
			xlabel="Z (mm)", ylabel="Frequency",
			title="Single Track Z Distribution",
			bins=50, legend=false)
	else
		plot(title="No Z data")
	end
	
	# Two track positions
	p3 = if length(rznubb.two_track_primary.xs) > 0
		scatter(rznubb.two_track_primary.xs, rznubb.two_track_primary.ys,
			xlabel="X (mm)", ylabel="Y (mm)",
			title="Two Track Primary XY",
			markersize=2, alpha=0.5, legend=false)
	else
		plot(title="No two track data")
	end
	
	p4 = if length(rznubb.two_track_primary.zs) > 0
		histogram(rznubb.two_track_primary.zs,
			xlabel="Z (mm)", ylabel="Frequency",
			title="Two Track Primary Z",
			bins=50, legend=false)
	else
		plot(title="No Z data")
	end
	
	plot(p1, p2, p3, p4, layout=(2,2), size=(800,600))
end

# ╔═╡ d4892a5a-4488-4244-9513-8ece86d59661
if load_znubb
	plot_position_distributions(rznubb)
else
	md"Load results to see position distributions."
end

# ╔═╡ 9c1d2e3f-a4b5-6c7d-8e9f-1a2b3c4d5e6f
if load_bkg
	plot_position_distributions(rbi214)
else
	md"Load $(bkg_display_name) results to see position distributions."
end

# ╔═╡ Cell order:
# ╠═04b446d6-f34f-11ed-2565-0b15d65b6781
# ╠═871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
# ╠═947c237c-9852-40e9-a83f-c23666db90aa
# ╠═349825ff-7ffe-4fa1-ba26-a772041f0323
# ╠═c9fc0547-0e73-4629-9909-e59c3d75169d
# ╠═3ae5d298-b9df-4f4e-85d9-d75e4aee1654
# ╟─6c59aeae-7990-4b43-8378-0de210a3291a
# ╠═ae89f5dc-e958-496a-91ac-0bd977355563
# ╠═4af9a3ef-e883-4bc3-a2f1-212102e4951b
# ╠═5caa7c21-82e5-420b-b4e7-a0e33543b74b
# ╠═8a55b4a3-5cbf-48c3-b150-2bd4ad73f440
# ╠═2b3c4d5e-f6a7-4b8c-9d1e-2f3a4b5c6d7e
# ╠═9e2f3d4e-a5b6-4c7d-8e9f-1a2b3c4d5f6a
# ╠═9f3e1cdb-758c-4223-ac86-1d63ac8f03da
# ╠═4d5e6f7a-b8c9-4d1e-2f3a-4b5c6d7e8f9a
# ╠═07b4e4c1-0469-41ca-9890-bb4f990b4645
# ╠═f60557f4-113e-44ab-ab53-56e967f8fda8
# ╠═5e6f7a8b-c9d1-4e2f-3a4b-5c6d7e8f9a1b
# ╟─144d60a3-d70f-442b-a252-76178fecdbf7
# ╠═0fa67f2a-50a5-4c82-b42b-07f45f14e914
# ╠═eba21011-2c2a-4540-b096-d39141d9abb2
# ╠═6624fb61-e8ac-41e5-a602-c8765e21cede
# ╠═6f7a8b9c-d1e2-4f3a-4b5c-6d7e8f9a1b2c
# ╠═7a8b9c1d-e2f3-4a5b-6c7d-8e9f1a2b3c4d
# ╟─474b25f8-bd95-4389-867a-bb753dc77d45
# ╠═6c58a746-da45-4e39-a178-91e060b2f34b
# ╠═aa71c155-9bed-4ad6-963e-575100094a9c
# ╠═6c148761-f9e2-4c7d-8d50-ee78bc0a8baf
# ╠═d1c1a8cb-a70e-480b-8f65-2cfa1a7022c3
# ╟─272c88a5-c6c0-41e1-b782-f8e23007eefc
# ╠═7fde840d-0c00-4593-b0ae-52d1990e8c4e
# ╟─7c20e495-0c39-44a5-9c40-8cd9c5d8e3de
# ╠═9c7ca8ed-5d3f-436f-bd82-675a03e6a5e9
# ╠═742b17ab-27fc-418c-9265-522fde5235bb
# ╠═1963d9dc-0cb9-4b8e-98ee-41491c6c784b
# ╠═b3634a2f-6600-4480-b316-2a79af18653e
# ╠═f301f5cf-0dc6-49d9-bd4d-ff117ea56a2e
# ╠═06bac1d3-c340-49ce-bd7c-92164e412590
# ╠═cd3ef36e-7150-47a5-9cec-9fa4d0e5779e
# ╟─37d7c197-0a10-4145-ab85-b5a22eae273e
# ╟─5248add6-0d0c-446d-8490-0ccc009fd6e9
# ╟─8b9c1d2e-f3a4-5b6c-7d8e-9f1a2b3c4d5e
# ╟─8c929794-64a4-48eb-bd85-ccfad6133b1e
# ╟─d4892a5a-4488-4244-9513-8ece86d59661
# ╟─9c1d2e3f-a4b5-6c7d-8e9f-1a2b3c4d5e6f
# ╟─c6bdf2b5-0b45-4e6a-af2a-b8f1c1ed4b78
# ╟─ee47c998-1649-4df2-bdc2-33bf53818e62
# ╟─1d2e3f4a-b5c6-7d8e-9f1a-2b3c4d5e6f7a
# ╟─21f70e07-ef56-4b84-87df-dd4e0d468ed1
# ╟─9a89ffdc-3889-48c0-a952-ad2c36094b4c
# ╟─a60ea1b2-1c4c-4565-ac34-d2580dd016e3
# ╟─995555f1-dd99-4903-a3b7-9e48b62d675d
# ╟─1362e7bf-0555-4551-9fcd-6047b7d9b555
# ╟─37110c55-a96e-4c4e-9ca1-9fb464865af0
# ╟─d115f04b-1cf2-47b4-81fd-d15aa8d6cd97
# ╟─d0f9e26b-2a3c-42a4-a826-7c8694f5d470
# ╠═1e6adccf-b640-4fd7-92cc-02cf61f27a24
# ╠═705979ca-69e9-48e1-9440-ef0ab8c66154
# ╠═ea1afc2f-9467-4b8b-a3e3-237925ac9fa1
# ╠═de61fd85-48cf-4c49-8abc-f929bdc4c3fc
