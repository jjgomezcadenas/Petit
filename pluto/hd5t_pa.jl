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
	#using JSON
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

# ╔═╡ 137d5f5b-e20b-407c-8137-48ba08c12bd2
function ingredients(path::String)
    # this is from the Julia source code (evalfile in base/loading.jl)
    # but with the modification that it returns the module instead of the last object
    name = Symbol(basename(path))
    m = Module(name)
    Core.eval(m,
        Expr(:toplevel,
                :(eval(x) = $(Expr(:core, :eval))($name, x)),
                :(include(x) = $(Expr(:top, :include))($name, x)),
                :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
                :(include($path))))
    m
end

# ╔═╡ 007b3b41-1290-45e3-8f3c-ea27069ee1db
begin
	jn = ingredients(string(pdir,"/src/Petit.jl"))
end

# ╔═╡ eb8abbaf-6d7a-45e6-85ce-5ac5bc6dcfe7
md"""
# Pre-selection
"""

# ╔═╡ b97eb6e0-40d9-4de8-bc0b-c2b1469cabab
dfstats =Dict("0nubb"=>1e+4,
	"bi214_copper_endcaps" => 5e+6,
	"bi214_copper_shell" => 5e+6,
	"bi214_cathode_volume" => 1e+6,
	"bi214_cathode_surface" => 5e+5,
	"bi214_ptfe_barrel" => 5e+6,
	"bi214_ptfe_endcap" => 5e+6,
	"tl208_copper_shell" => 5e+5,
	"tl208_copper_endcaps" => 5e+5,
	"tl208_cathode_volume" => 1e+5,
	"tl208_ptfe_barrel"=> 1e+5,
	"tl208_ptfe_endcap"=> 1e+5,
	"electron_0nubb_energy" => 1e+4)

# ╔═╡ e093fcd9-5c60-4308-b2df-04151f40b2c3
md"""
| File | Size data set |
|:-----|----------:|
| 0nubb.h5 | $(dfstats["0nubb"]) | 
| bi214\_copper\_endcaps.h5 | $(dfstats["bi214_copper_endcaps"])|
| bi214\_ptfe\_endcap.h5 | $(dfstats["bi214_ptfe_endcap"]) |
| bi214\_copper\_shell.h5 | $(dfstats["bi214_copper_shell"]) |
| bi214\_ptfe\_barrel.h5 | $(dfstats["bi214_ptfe_barrel"]) |
| bi214\_cathode\_volume.h5: | $(dfstats["bi214_cathode_volume"]) |
| bi214\_cathode\_surface.h5 | $(dfstats["bi214_cathode_surface"]) |
| tl208\_copper\_shell.h5    | $(dfstats["tl208_copper_shell"]) |
| tl208\_copper\_endcaps.h5  | $(dfstats["tl208_copper_endcaps"]) |
| tl208\_cathode\_volume.h5  | $(dfstats["tl208_cathode_volume"]) |
| tl208\_ptfe\_barrel.h5  | $(dfstats["tl208_ptfe_barrel"]) |
| tl208\_ptfe\_endcap.h5  | $(dfstats["tl208_ptfe_endcap"]) |
| electron\_0nubb\_energy.h5 | $(dfstats["electron_0nubb_energy"]) |
    """

# ╔═╡ e71f7afb-ec8a-4c69-ab46-74813b9aaf7c
#| bi214_copper.h5 | $(dfstats["bi214_copper"]) |
#| bi214_surface.h5 | $(dfstats["bi214_surface"])| 
#| tl208_copper.h5 | $(dfstats["tl208_copper"])|
#| tl208_surface.h5 | $(dfstats["tl208_surface"]) |

# ╔═╡ ce643771-65d5-45bf-9bd7-a8509ffa508f
begin
	energy_window_low = 2300.0 # keV 
	energy_window_up =2700.0
	md"""
	- Analysis performed in energy window: =($(energy_window_low), $(energy_window_up)) keV
	"""
end

# ╔═╡ 6c59aeae-7990-4b43-8378-0de210a3291a
md"""
## Load Pre-selection data

"""

# ╔═╡ 4af9a3ef-e883-4bc3-a2f1-212102e4951b
begin
	# Default results base directory (matches new script structure)
	results_base_dir = joinpath(pdir, "AnalysisResults")
	results_znubb = joinpath(results_base_dir, "0nubb")
	
	md"""
	Loading results from: $(results_znubb)
	"""
end

# ╔═╡ 5caa7c21-82e5-420b-b4e7-a0e33543b74b
md"Load znubb results? $(@bind load_znubb CheckBox(default=true))"

# ╔═╡ 8a55b4a3-5cbf-48c3-b150-2bd4ad73f440
if load_znubb
	rznubb = jn.Petit.load_analysis_results(results_znubb)
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

# ╔═╡ 0104c472-a737-422c-b345-a4b11ef1fc76
#"bi214_copper" => "Bi214 Copper",
#"bi214_surface" => "Bi214 Surface",
#"tl208_copper" => "Tl208 Copper",
#"tl208_surface" => "Tl208 Surface",

# ╔═╡ 9e2f3d4e-a5b6-4c7d-8e9f-1a2b3c4d5f6a
@bind background_file Select([
	"bi214_copper_shell" => "Bi214 Copper Shell",
	"bi214_copper_endcaps" => "Bi214 Copper End Cap",
	"bi214_cathode_volume" => "Bi214 Cathode Volume",
	"bi214_cathode_surface" => "Bi214 Cathode Surface",
	"bi214_ptfe_barrel" => "Bi214 PTFE Barrel",
	"bi214_ptfe_endcap" => "Bi214 PTFE End Cap",
	"tl208_copper_shell" => "Tl208 Copper Shell",
	"tl208_copper_endcaps" => "Tl208 Copper End Cap",
	"tl208_cathode_volume" => "Tl208 Cathode Volume",
	"tl208_ptfe_barrel"=> "Tl208 PTFE Barrel",
	"tl208_ptfe_endcap"=> "Tl208 PTFE End Cap",
	"electron_0nubb_energy" => "Electron 0νββ Energy"
	
], default="tl208_copper_endcaps")

# ╔═╡ 9f3e1cdb-758c-4223-ac86-1d63ac8f03da
begin
	# Construct the results directory path based on selection
	results_bkg = joinpath(results_base_dir, background_file)
	# Extract display name from the file name
	bkg_display_name = replace(background_file, "_" => " ") |> titlecase
	md"""
	Load $(bkg_display_name) results? $(@bind load_bkg CheckBox(default=true))
	
	Results directory: $(results_bkg)
	"""
end

# ╔═╡ 4d5e6f7a-b8c9-4d1e-2f3a-4b5c6d7e8f9a
if load_bkg
	if isdir(results_bkg)
		rbi214 = jn.Petit.load_analysis_results(results_bkg)
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
		⚠️ Results directory not found: $(results_bkg)
		
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
## Pre-selection Summary
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
	### Event Classification Summary ($(bkg_display_name))
	
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
# Analysis
"""

# ╔═╡ 9a53fd38-82b4-486c-96e7-d756a9eb3ffa
begin
	evt1t = rznubb.single_track[1,2]
	plot_hits_evt(rznubb.single_track, evt1t; nbins=100)
end

# ╔═╡ 1eeaa3a9-2892-4fed-8ad6-70d98d6f9fe5
md"Radial cut (mm): $(@bind xr NumberField(1000.0:50.0:2000.0, default=1300.0))"

# ╔═╡ 9c8bf44c-a996-4d6b-997c-f680e933621f
md"""
### Fiducial Radial cut
- Cut: All hits radial position $\sqrt{(x^2+y^2)}$ below $(xr) mm
"""

# ╔═╡ 7cdc518c-cd96-4315-97ff-ab4ced327469
begin
	znust = jn.Petit.filter_radial(rznubb.single_track, xr, xr)
	nznust = number_of_events(znust)
	md"""
	- Radial cut > $(xr)
	- Cut efficency =$(nznust/rznubb.n_single_track)
	"""
end

# ╔═╡ 0fa67f2a-50a5-4c82-b42b-07f45f14e914
md"""
## Single Track Events (0νββ): e = $(round(mean(znust.energy), digits=2))
- Cut: A single track reconstructed in the event, no floating energy
"""

# ╔═╡ 3fbe8c31-9638-451a-83e3-de66af0af6c3
md"""
### Cut in track length and fiducial cut in Z
"""

# ╔═╡ 3ac9c660-51cc-452f-bacd-58cdc3915c23
md"Track length cut (mm): $(@bind trklm NumberField(25:200, default=50))"

# ╔═╡ 60d6bd0e-fa00-4312-9c7c-a859e2232e2b
begin
	znuxs = jn.Petit.filter_short_tracks(znust, trklm)
	nznuxs = number_of_events(znuxs)
	md"""
	- Cut track length > $(trklm)
	- Cut efficency =$(nznuxs/nznust)
	"""
end

# ╔═╡ bf6e00e2-2f81-46ab-bd6b-d4c8c2afbf60
md"""
- Z cut
"""

# ╔═╡ 1bf41951-4db8-4c28-ae18-1946ebf8489b
md"Z inner left: $(@bind zil NumberField(-500.0:0.0, default=0.0))"

# ╔═╡ 68c093a5-c555-4dda-bea4-e8ff79b11e3e
md"Z inner right: $(@bind zir NumberField(0.0:500.0, default=0.0))"

# ╔═╡ 01c544b2-a655-4d68-880f-7bdb2944ca0d
md"Z outer left: $(@bind zol NumberField(-2000.0:-1500.0, default=-1950.0))"

# ╔═╡ 42f671cf-7d2c-4efb-9c2b-6a90226a8323
md"Z outer right: $(@bind zor NumberField(1500.0:2000.0, default=1950.0))"

# ╔═╡ 90bc6dd6-2904-4b16-9419-f9fc8147772c
begin
	znuxz = jn.Petit.filter_z(znuxs, zil, zir, zol, zor)
	nznuxz = number_of_events(znuxz)
	md"""
	- Z cut:
	- Cut efficency =$(nznuxz/nznuxs)
	"""
end

# ╔═╡ 474b25f8-bd95-4389-867a-bb753dc77d45
md"""
### Energy Resolution Analysis
"""

# ╔═╡ 6624fb61-e8ac-41e5-a602-c8765e21cede
if load_znubb && length(znuxz.energy) > 0
	znue = combine(groupby(znuxz, :event_id), :energy => first => :energy)
	eht1, pht1 = jn.Petit.step_hist(znue.energy;
         nbins = 40,
         xlim   = (energy_window_low, energy_window_up),
         xlabel = "E (keV)",
         ylabel = "Frequency",
         title="Single Track Energy Distribution")
	plot(pht1)
else
	md"No single track events found."
end

# ╔═╡ 6c58a746-da45-4e39-a178-91e060b2f34b
md"Energy resolution (keV): $(@bind erex NumberField(0.0:0.5:30.0, default=12.5))"

# ╔═╡ aa71c155-9bed-4ad6-963e-575100094a9c
md"ROI low (keV): $(@bind roi_low NumberField(2400.0:2500.0, default=2400.0))"

# ╔═╡ 6c148761-f9e2-4c7d-8d50-ee78bc0a8baf
md"ROI up (keV): $(@bind roi_up NumberField(2475.0:2550.0, default=2500.0))"

# ╔═╡ d542d7d6-01ad-4b57-a515-da6569609a8f


# ╔═╡ 7c20e495-0c39-44a5-9c40-8cd9c5d8e3de
md"""
### Signal Efficiency vs ROI
"""

# ╔═╡ 87705166-150b-4dea-903f-d2467dfe8d3b
md" Energy ROI (low) -keV: $(@bind roi_xi NumberField(2400.0:2500.0, default=2440.0))"

# ╔═╡ b426a406-6c4a-4965-8f79-6557ed624e4b
md" Energy ROI (up) -keV: $(@bind roi_xu NumberField(2400.0:2500.0, default=2480.0))"

# ╔═╡ d1c1a8cb-a70e-480b-8f65-2cfa1a7022c3
if load_znubb && length(znuxz.energy) > 0
	eres = jn.Petit.smear_histogram(eht1, erex)
	ehrx, phrx = jn.Petit.step_hist(eres;
         nbins = 80,
         xlim   = (energy_window_low, energy_window_up),
         xlabel = "E (keV)",
         ylabel = "Frequency",
         title="bb0nu Energy (with $(erex) keV resolution)")
	plot(phrx)
	vline!([roi_xi], color=:red, linewidth=1, linestyle=:dash, label="")
	vline!([roi_xu], color=:red, linewidth=1, linestyle=:dash, label="")
else
	md"No data to apply resolution."
end

# ╔═╡ 0253f041-218e-48be-acab-edeff1a53510
begin
	step = 5.0
	xx = roi_low:step:roi_up
end

# ╔═╡ 272c88a5-c6c0-41e1-b782-f8e23007eefc
if load_znubb && length(znuxz.energy) > 0
	counts_roi = jn.Petit.counts_in_range(ehrx, roi_xi, roi_xu)
	eff_roi = counts_roi/nznuxz
	md"""
	**Counts in ROI [$(roi_xi), $(roi_xu)] keV:** $(counts_roi)
	
	**Efficiency:** $(round(eff_roi, digits=3))%
	"""
else
	md"No data for ROI analysis."
end

# ╔═╡ f301f5cf-0dc6-49d9-bd4d-ff117ea56a2e
begin
	signal_contained = rznubb.n_events_processed/dfstats["0nubb"]
	signal_1trk = rznubb.n_single_track/rznubb.n_events_processed
	signal_radial = nznust/rznubb.n_single_track
	signal_trkl = nznuxs/nznust
	signal_zfid = nznuxz/nznuxs
	signal_roi = eff_roi
	signal_2e = 0.8 # From Josh paper
	#signal_roi = integrated_signal_eff(ehrx, energy_window_low, energy_window_up, roi_xi, roi_xu)[3]
	signal_teff = signal_contained * signal_1trk * signal_radial * signal_trkl * signal_zfid * signal_roi * signal_2e 
	md"""
	### Signal efficiency ($(erex) keV FWHM resolution)
	- Fraction of events contained $(signal_contained)
	- Fraction of events 1 track $(round(signal_1trk, digits=2))
	- Fraction of events Radial cut $(round(signal_radial, digits=2))
	- Fraction of events Track Length cut $(round(signal_trkl, digits=2))
	- Fraction of events Z Length cut $(round(signal_zfid, digits=2))
	- Fraction of events ROIS $(round(signal_roi, digits=2))
	- Fraction of events 2 electron ID $(signal_2e)
	- **Total signal efficiency: $(round(signal_teff, digits=2))**
	"""
end

# ╔═╡ 9df845b9-5608-4624-94bc-a8e5dab589c9
begin
	path_signal = "AnalysisSummary/znubb.js"
	jn.Petit.json_analysis_summary(path_signal, "znubb", trklm, xr,
						  zil, zir, zol, zor,
						  erex, roi_xi, roi_xu,
						  signal_contained,
						  signal_1trk,
						  signal_radial,
						  signal_trkl,
						  signal_zfid,
						  signal_roi,
						  signal_2e,  
						  signal_teff)
	md"""
	#### Save summary in $(path_signal)
	"""
end

# ╔═╡ 6f7a8b9c-d1e2-4f3a-4b5c-6d7e8f9a1b2c
md"""
## Single Track Events ($(bkg_display_name)) : e = $(round(mean(rbi214.single_track.energy), digits=2))
"""

# ╔═╡ 7697451d-851d-4805-b7e2-ef6e47ca4cb0
begin
	evtbkg = rbi214.single_track[1,2]
	plot_hits_evt(rbi214.single_track, evtbkg; nbins=100)
end

# ╔═╡ cd8e462f-307e-4f86-94e9-1d01ebd4c642
md"""
### Fiducial Radial cut
- Cut: All hits radial position $\sqrt{(x^2+y^2)}$ below $(xr) mm
"""

# ╔═╡ e9ec8fad-d7e2-47b2-9e95-bb432510b8fe
begin
	bkgnd = jn.Petit.filter_radial(rbi214.single_track, xr, xr)
	nbkgnd = number_of_events(bkgnd)
	md"""
	- Radial cut > $(xr)
	- Cut efficency =$(nbkgnd/rbi214.n_single_track)
	"""
end

# ╔═╡ dbd08aa8-3b1a-4620-aea5-38391fd86632
md"""
### Cut in track length and fiducial cut in Z
"""

# ╔═╡ 5401d3be-650f-464a-bd53-47d5c1e03d35
begin
	bkxs = jn.Petit.filter_short_tracks(bkgnd, trklm)
	nbkxs = number_of_events(bkxs)
	md"""
	- Cut track length > $(trklm)
	- Cut efficency =$(nbkxs/nbkgnd)
	"""
end

# ╔═╡ e317c79e-d1ad-4146-8614-6223598fdd7f
begin
	bkxz = jn.Petit.filter_z(bkxs, zil, zir, zol, zor)
	nbkxz = number_of_events(bkxz)
	md"""
	- Z cut:
	- Cut efficency =$(nbkxz/nbkxs)
	"""
end

# ╔═╡ 7a8b9c1d-e2f3-4a5b-6c7d-8e9f1a2b3c4d
if load_bkg && length(Set(bkxz.energy)) > 0
	bkge = combine(groupby(bkxz, :event_id), :energy => first => :energy)
	eht1_bi214, pht1_bi214 = jn.Petit.step_hist(bkge.energy;
         nbins = 40,
         xlim   = (energy_window_low, energy_window_up),
         xlabel = "E (keV)",
         ylabel = "Frequency",
         title="Single Track Energy Distribution (Bi214)")
	plot(pht1_bi214)
else
	md"No single track events found for Bi214."
end

# ╔═╡ 93b19995-01d4-4178-803d-f380c0812099
sum(eht1_bi214.weights)

# ╔═╡ 7fde840d-0c00-4593-b0ae-52d1990e8c4e
if load_bkg && length(bkxz.energy) > 0
	eres_bi214 = jn.Petit.smear_histogram(eht1_bi214, erex)
	ehrx_bi214, phrx_bi214 = jn.Petit.step_hist(eres_bi214;
         nbins = 60,
         xlim   = (energy_window_low, energy_window_up),
         xlabel = "E (keV)",
         ylabel = "Frequency",
         title="Bi-214 Energy (with $(erex) keV resolution)")
	plot(phrx_bi214)
	vline!([roi_xi], color=:red, linewidth=1, linestyle=:dash, label="")
	vline!([roi_xu], color=:red, linewidth=1, linestyle=:dash, label="")
else
	md"No data to apply resolution."
end

# ╔═╡ b3634a2f-6600-4480-b316-2a79af18653e
#if load_bkg && length(bkgnd.energy) > 0
#	eroi_bi214 = integrated_signal_eff(ehrx_bi214, energy_window_low, energy_window_up, roi_xi, roi_xu)[3]
#	md"""
#	
#	- Selection efficiency ROI ($(roi_xi), $(roi_xu)) = $(round(eroi_bi214, digits=2))
	
#	"""
#end

# ╔═╡ d1325000-6afd-4eb3-a2ab-6edae9a21c57
if load_bkg && length(bkxz.energy) > 0
	counts_roi_bk = jn.Petit.counts_in_range(ehrx_bi214, roi_xi, roi_xu)
	
	if counts_roi_bk == 0.0  # Take a Poisson approach for an upper limit!!! 
		counts_roi_bk = 1.0
	end
	eff_roi_bk = counts_roi_bk/nbkxz
	md"""
	**Counts in ROI [$(roi_xi), $(roi_xu)] keV:** $(counts_roi_bk)
	
	**Efficiency:** $(round(eff_roi_bk, digits=3))
	"""
else
	md"No data for ROI analysis."
end

# ╔═╡ 06bac1d3-c340-49ce-bd7c-92164e412590
begin
	bkgnd_contained = rbi214.n_events_processed/dfstats[background_file]
	bkgnd_1trk = rbi214.n_single_track/rbi214.n_events_processed
	bkgnd_radial = nbkgnd/rbi214.n_single_track
	bkgnd_trkl = nbkxs/nbkgnd
	bkgnd_zfid = nbkxz/nbkxs
	bkgnd_roi = eff_roi_bk
	bkgnd_2e    = 0.04 # from Josh
	bkgnd_teff = bkgnd_contained * bkgnd_1trk * bkgnd_radial * bkgnd_trkl *     bkgnd_zfid * bkgnd_roi * bkgnd_2e 
	
	md"""
	### Background efficiency ($(erex) keV FWHM resolution)
	- Fraction of events contained $(bkgnd_contained)
	- Fraction of events 1 track $(round(bkgnd_1trk, digits=3))
	- Fraction of events Radial cut $(round(bkgnd_radial, digits=2))
	- Fraction of events Track Length cut $(round(bkgnd_trkl, digits=2))
	- Fraction of events Z Length cut $(round(bkgnd_zfid, digits=2))
	- Fraction of events ROIS $(round(bkgnd_roi, digits=3))
	- Fraction of events 2 electron ID $(bkgnd_2e)
	- **Total signal efficiency: $(round(bkgnd_teff, digits=8))**
	"""
end

# ╔═╡ 438ce4fe-a1ca-4b1a-b120-b169f508728d
begin
	path_bkg = string("AnalysisSummary/", background_file, ".js")
	json_analysis_summary(path_bkg, background_file, trklm, xr,
						  zil, zir, zol, zor,
						  erex, roi_xi, roi_xu,
						  bkgnd_contained,
						  bkgnd_1trk,
						  bkgnd_radial,
						  bkgnd_trkl,
						  bkgnd_zfid,
						  bkgnd_roi,
						  bkgnd_2e,  
						  bkgnd_teff)
	md"""
	#### Save summary in $(path_bkg)
	"""
end

# ╔═╡ 37d7c197-0a10-4145-ab85-b5a22eae273e
md"""
## Track and Multi track Plots
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
		push!(countx, jn.Petit.counts_in_range(ehrx, rx, rup))
	end
	countx = countx /norm
end

# ╔═╡ 9c7ca8ed-5d3f-436f-bd82-675a03e6a5e9
if load_znubb && length(znuxz.energy) > 0
	countx = signal_eff(ehrx, roi_low, roi_up; step=step)
	plot(xx, countx,
		xlabel="ROI Lower Bound (keV)",
		ylabel="Efficiency",
		title="Signal Efficiency vs ROI Lower Bound",
		marker=:circle,
		legend=false)
else
	md"No data for efficiency plot."
end

# ╔═╡ 1963d9dc-0cb9-4b8e-98ee-41491c6c784b
if load_bkg && length(bkxz.energy) > 0
	countx_bi214 = signal_eff(ehrx_bi214, roi_low, roi_up; step=step)
	plot(xx, countx_bi214,
		xlabel="ROI Lower Bound (keV)",
		ylabel="Efficiency",
		title="Background Efficiency vs ROI Lower Bound",
		marker=:circle,
		legend=false)
else
	md"No data for efficiency plot."
end

# ╔═╡ 705979ca-69e9-48e1-9440-ef0ab8c66154
function integrated_signal_eff(ehrx, rlow, rup, rmin, rmax) 
	n1 = jn.Petit.counts_in_range(ehrx, rlow, rup)
	n2 = jn.Petit.counts_in_range(ehrx, rmin, rmax)
	return n1, n2, n2/n1
end

# ╔═╡ 742b17ab-27fc-418c-9265-522fde5235bb
if load_znubb && length(bkxz.energy) > 0
	eroi = integrated_signal_eff(ehrx, roi_low, roi_up, roi_xi, roi_xu)[3]
	
	md"""
	- Selection efficiency ROI ($(roi_xi), $(roi_xu)) = $(round(eroi, digits=3))
	
	"""
end

# ╔═╡ 696f320d-0454-401d-ac65-12150467080b
function histogram_track_length(hitsdf::DataFrame; trkl=200.0)
    # Group by event_id and count the number of rows (hits) in each group
    grouped_df = groupby(hitsdf, :event_id)
    
    # Calculate track length for each event
    result = combine(grouped_df, nrow => :track_length)
	htl, ptl = jn.Petit.step_hist(result.track_length;
	         nbins = 50,
	         xlim   = (0.0, trkl),
	         xlabel = "Track length (in mm)",
	         ylabel = "Frequency",
	         title="track length")
end

# ╔═╡ 2f05b56f-9f10-429a-93d9-1c0dd3eb26f0
begin
	htl, ptl = histogram_track_length(znust; trkl=200.0)
	plot(ptl)
end

# ╔═╡ b012a0c2-5de4-441e-b638-da85546a0c72
begin
	htlb, ptlb = histogram_track_length(bkgnd; trkl=200.0)
	plot(ptlb)
end

# ╔═╡ ea1afc2f-9467-4b8b-a3e3-237925ac9fa1
function histogram_energies_trks(rznubb; esec=500.0)
	energies = collect(Set(rznubb.single_track.energy))
	p1 = if length(energies) > 0
		_, p = jn.Petit.step_hist(energies;
	         nbins = 50,
	         xlim   = (0.0, 2700.0),
	         xlabel = "E (keV)",
	         ylabel = "Frequency",
	         title="Single track")
		p
	else
		plot(title="No single track data")
	end
	energies = collect(Set(rznubb.two_track_primary.energy))
	p2 = if length(energies) > 0
		_, p = jn.Petit.step_hist(energies;
	         nbins = 50,
	         xlim   = (0.0, 2700.0),
	         xlabel = "E (keV)",
	         ylabel = "Frequency",
	         title="Two tracks (primary)")
		p
	else
		plot(title="No two track data")
	end
	energies = collect(Set(rznubb.two_track_secondary.energy))
	p3 = if length(energies) > 0
		_, p = jn.Petit.step_hist(energies;
	         nbins = 50,
	         xlim   = (0.0, esec),
	         xlabel = "E (keV)",
	         ylabel = "Frequency",
	         title="Two tracks (secondary)")
		p
	else
		plot(title="No secondary data")
	end
	energies = collect(Set(rznubb.three_track_secondary.energy))
	p4 = if length(energies) > 0		
		_, p = jn.Petit.step_hist(energies;
	         nbins = 50,
	         xlim   = (0.0, esec),
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
function plot_position_distributions(rznubb; xc=0.0, yc=0.0, zil=0.0, zir=0.0, 
                                     zol=-2000.0, zor=2000.0)
      # Calculate radius for the circle
      r = sqrt(xc^2 + yc^2)

      # Generate circle points
      θ = range(0, 2π, length=100)
      circle_x = r .* cos.(θ)
      circle_y = r .* sin.(θ)

      # Single track positions
      p1 = if length(rznubb.single_track.x) > 0
          scatter(rznubb.single_track.x, rznubb.single_track.y,
              xlabel="X (mm)", ylabel="Y (mm)",
              title="Single Track XY",
              markersize=2, alpha=0.5, legend=false)
          plot!(circle_x, circle_y, color=:red, linewidth=2, label="")
      else
          plot(title="No single track data")
      end

      p2 = if length(rznubb.single_track.z) > 0
          histogram(rznubb.single_track.z,
              xlabel="Z (mm)", ylabel="Frequency",
              title="Single Track Z Distribution",
              bins=50, legend=false)
          # Add vertical lines for z boundaries
          vline!([zol], color=:red, linewidth=2, linestyle=:solid, label="")
          vline!([zil], color=:red, linewidth=2, linestyle=:solid, label="")
          vline!([zir], color=:red, linewidth=2, linestyle=:solid, label="")
          vline!([zor], color=:red, linewidth=2, linestyle=:solid, label="")
      else
          plot(title="No Z data")
      end

      # Two track positions
      p3 = if length(rznubb.two_track_primary.x) > 0
          scatter(rznubb.two_track_primary.x, rznubb.two_track_primary.y,
              xlabel="X (mm)", ylabel="Y (mm)",
              title="Two Track Primary XY",
              markersize=2, alpha=0.5, legend=false)
          plot!(circle_x, circle_y, color=:red, linewidth=2, label="")
      else
          plot(title="No two track data")
      end

      p4 = if length(rznubb.two_track_primary.z) > 0
          histogram(rznubb.two_track_primary.z,
              xlabel="Z (mm)", ylabel="Frequency",
              title="Two Track Primary Z",
              bins=50, legend=false)
          # Add vertical lines for z boundaries
          vline!([zol], color=:red, linewidth=2, linestyle=:solid, label="")
          vline!([zil], color=:red, linewidth=2, linestyle=:solid, label="")
          vline!([zir], color=:red, linewidth=2, linestyle=:solid, label="")
          vline!([zor], color=:red, linewidth=2, linestyle=:solid, label="")
      else
          plot(title="No Z data")
      end

      plot(p1, p2, p3, p4, layout=(2,2), size=(800,600))
  end

# ╔═╡ d4892a5a-4488-4244-9513-8ece86d59661
if load_znubb
	plot_position_distributions(rznubb, xc=1300.0, yc=1300.0, zil=zil, zir=zir, 
                                     zol=zol, zor=zor)
else
	md"Load results to see position distributions."
end

# ╔═╡ 9c1d2e3f-a4b5-6c7d-8e9f-1a2b3c4d5e6f
if load_bkg
	plot_position_distributions(rbi214, xc=1300.0, yc=1300.0, zil=zil, zir=zir, 
                                     zol=zol, zor=zor)
else
	md"Load $(bkg_display_name) results to see position distributions."
end

# ╔═╡ Cell order:
# ╠═04b446d6-f34f-11ed-2565-0b15d65b6781
# ╠═871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
# ╠═947c237c-9852-40e9-a83f-c23666db90aa
# ╠═349825ff-7ffe-4fa1-ba26-a772041f0323
# ╠═c9fc0547-0e73-4629-9909-e59c3d75169d
# ╠═137d5f5b-e20b-407c-8137-48ba08c12bd2
# ╠═007b3b41-1290-45e3-8f3c-ea27069ee1db
# ╟─eb8abbaf-6d7a-45e6-85ce-5ac5bc6dcfe7
# ╠═b97eb6e0-40d9-4de8-bc0b-c2b1469cabab
# ╠═e093fcd9-5c60-4308-b2df-04151f40b2c3
# ╟─e71f7afb-ec8a-4c69-ab46-74813b9aaf7c
# ╟─ce643771-65d5-45bf-9bd7-a8509ffa508f
# ╟─6c59aeae-7990-4b43-8378-0de210a3291a
# ╠═4af9a3ef-e883-4bc3-a2f1-212102e4951b
# ╠═5caa7c21-82e5-420b-b4e7-a0e33543b74b
# ╟─8a55b4a3-5cbf-48c3-b150-2bd4ad73f440
# ╟─2b3c4d5e-f6a7-4b8c-9d1e-2f3a4b5c6d7e
# ╟─0104c472-a737-422c-b345-a4b11ef1fc76
# ╠═9e2f3d4e-a5b6-4c7d-8e9f-1a2b3c4d5f6a
# ╟─9f3e1cdb-758c-4223-ac86-1d63ac8f03da
# ╟─4d5e6f7a-b8c9-4d1e-2f3a-4b5c6d7e8f9a
# ╠═07b4e4c1-0469-41ca-9890-bb4f990b4645
# ╟─f60557f4-113e-44ab-ab53-56e967f8fda8
# ╟─5e6f7a8b-c9d1-4e2f-3a4b-5c6d7e8f9a1b
# ╠═144d60a3-d70f-442b-a252-76178fecdbf7
# ╠═0fa67f2a-50a5-4c82-b42b-07f45f14e914
# ╠═9a53fd38-82b4-486c-96e7-d756a9eb3ffa
# ╠═9c8bf44c-a996-4d6b-997c-f680e933621f
# ╠═1eeaa3a9-2892-4fed-8ad6-70d98d6f9fe5
# ╠═7cdc518c-cd96-4315-97ff-ab4ced327469
# ╠═d4892a5a-4488-4244-9513-8ece86d59661
# ╠═3fbe8c31-9638-451a-83e3-de66af0af6c3
# ╠═2f05b56f-9f10-429a-93d9-1c0dd3eb26f0
# ╠═3ac9c660-51cc-452f-bacd-58cdc3915c23
# ╠═60d6bd0e-fa00-4312-9c7c-a859e2232e2b
# ╠═bf6e00e2-2f81-46ab-bd6b-d4c8c2afbf60
# ╠═1bf41951-4db8-4c28-ae18-1946ebf8489b
# ╠═68c093a5-c555-4dda-bea4-e8ff79b11e3e
# ╠═01c544b2-a655-4d68-880f-7bdb2944ca0d
# ╠═42f671cf-7d2c-4efb-9c2b-6a90226a8323
# ╠═90bc6dd6-2904-4b16-9419-f9fc8147772c
# ╟─474b25f8-bd95-4389-867a-bb753dc77d45
# ╠═6624fb61-e8ac-41e5-a602-c8765e21cede
# ╠═6c58a746-da45-4e39-a178-91e060b2f34b
# ╠═aa71c155-9bed-4ad6-963e-575100094a9c
# ╠═6c148761-f9e2-4c7d-8d50-ee78bc0a8baf
# ╠═d542d7d6-01ad-4b57-a515-da6569609a8f
# ╠═d1c1a8cb-a70e-480b-8f65-2cfa1a7022c3
# ╠═7c20e495-0c39-44a5-9c40-8cd9c5d8e3de
# ╠═87705166-150b-4dea-903f-d2467dfe8d3b
# ╠═b426a406-6c4a-4965-8f79-6557ed624e4b
# ╠═0253f041-218e-48be-acab-edeff1a53510
# ╠═9c7ca8ed-5d3f-436f-bd82-675a03e6a5e9
# ╠═272c88a5-c6c0-41e1-b782-f8e23007eefc
# ╠═f301f5cf-0dc6-49d9-bd4d-ff117ea56a2e
# ╠═9df845b9-5608-4624-94bc-a8e5dab589c9
# ╠═6f7a8b9c-d1e2-4f3a-4b5c-6d7e8f9a1b2c
# ╠═7697451d-851d-4805-b7e2-ef6e47ca4cb0
# ╠═cd8e462f-307e-4f86-94e9-1d01ebd4c642
# ╠═e9ec8fad-d7e2-47b2-9e95-bb432510b8fe
# ╠═9c1d2e3f-a4b5-6c7d-8e9f-1a2b3c4d5e6f
# ╠═dbd08aa8-3b1a-4620-aea5-38391fd86632
# ╠═b012a0c2-5de4-441e-b638-da85546a0c72
# ╠═5401d3be-650f-464a-bd53-47d5c1e03d35
# ╠═e317c79e-d1ad-4146-8614-6223598fdd7f
# ╠═7a8b9c1d-e2f3-4a5b-6c7d-8e9f1a2b3c4d
# ╠═93b19995-01d4-4178-803d-f380c0812099
# ╠═7fde840d-0c00-4593-b0ae-52d1990e8c4e
# ╠═742b17ab-27fc-418c-9265-522fde5235bb
# ╠═1963d9dc-0cb9-4b8e-98ee-41491c6c784b
# ╟─b3634a2f-6600-4480-b316-2a79af18653e
# ╠═d1325000-6afd-4eb3-a2ab-6edae9a21c57
# ╠═06bac1d3-c340-49ce-bd7c-92164e412590
# ╠═438ce4fe-a1ca-4b1a-b120-b169f508728d
# ╠═37d7c197-0a10-4145-ab85-b5a22eae273e
# ╠═5248add6-0d0c-446d-8490-0ccc009fd6e9
# ╠═8b9c1d2e-f3a4-5b6c-7d8e-9f1a2b3c4d5e
# ╟─d0f9e26b-2a3c-42a4-a826-7c8694f5d470
# ╠═1e6adccf-b640-4fd7-92cc-02cf61f27a24
# ╠═705979ca-69e9-48e1-9440-ef0ab8c66154
# ╠═696f320d-0454-401d-ac65-12150467080b
# ╠═ea1afc2f-9467-4b8b-a3e3-237925ac9fa1
# ╠═de61fd85-48cf-4c49-8abc-f929bdc4c3fc
