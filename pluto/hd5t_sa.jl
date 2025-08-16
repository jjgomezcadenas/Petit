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
	using DataFrames
	using CSV
	using Printf
	using Markdown
	using InteractiveUtils
	using OrderedCollections
end

# ╔═╡ c9fc0547-0e73-4629-9909-e59c3d75169d
begin
	using Petit
end

# ╔═╡ 16a23b13-537b-45dd-b90f-7220e3d969d9
 include("hd5t_sa_functions.jl")

# ╔═╡ 04b446d6-f34f-11ed-2565-0b15d65b6781
PlutoUI.TableOfContents(title="HD5t Summary Analysis", indent=true)

# ╔═╡ 871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
begin
	pdir = joinpath(ENV["PROJECTS"], "Petit")
	summary_dir = joinpath(pdir, "AnalysisSummary")
end

# ╔═╡ c17d78d1-39a1-40da-a3bc-d0a8717a4e10
begin
      using Pkg
      Pkg.activate(pdir)
      Pkg.instantiate()
  end

# ╔═╡ eb8abbaf-6d7a-45e6-85ce-5ac5bc6dcfe7
md"""
# Analysis Summary

This notebook reads and displays the analysis summaries from JSON files.
"""

# ╔═╡ 1629b0c1-64a3-4b6e-b24b-45f423f8a66b
summary_dir

# ╔═╡ b97eb6e0-40d9-4de8-bc0b-c2b1469cabab
md"""
## Available Summary Files
"""

# ╔═╡ 3f2262a1-33e5-4fe2-a1ef-c03df3d1aa4b
summary_files = find_summary_files(summary_dir)

# ╔═╡ 8a55b4a3-5cbf-48c3-b150-2bd4ad73f440
md"""
## Select Summary File

Choose which analysis summary to display:
"""

# ╔═╡ 9e2f3d4e-a5b6-4c7d-8e9f-1a2b3c4d5f6a
if length(summary_files) > 0
	@bind selected_file Select(summary_files, default=summary_files[1])
else
	@bind selected_file Select(["No files available"])
end

# ╔═╡ f60557f4-113e-44ab-ab53-56e967f8fda8
if length(summary_files) > 0 && selected_file != "No files available"
	begin
		# Read the selected JSON file
		file_path = joinpath(summary_dir, selected_file)
		summary_data = read_json_analysis_summary(file_path)
		
		md"""
		## File: $(selected_file)
		"""
	end
else
	md"""
	No file selected or available.
	"""
end

# ╔═╡ 0fa67f2a-50a5-4c82-b42b-07f45f14e914
if length(summary_files) > 0 && selected_file != "No files available"
	begin
		# Create cuts parameters table
		cuts_df = DataFrame(
			Parameter = [
				"Data Type",
				"Energy Resolution (keV)",
				"Track Length Cut (nhits)",
				"Radial Cut (mm)",
				"Z Inner Left (mm)",
				"Z Inner Right (mm)",
				"Z Outer Left (mm)",
				"Z Outer Right (mm)",
				"ROI Lower Bound (keV)",
				"ROI Upper Bound (keV)"
			],
			Value = [
				summary_data["data_type"],
				summary_data["energy_resolution_keV"],
				summary_data["cut_track_length_nhits"],
				summary_data["cut_radial_mm"],
				summary_data["cut_z_inner_left_mm"],
				summary_data["cut_z_inner_right_mm"],
				summary_data["cut_z_outer_left_mm"],
				summary_data["cut_z_outer_right_mm"],
				summary_data["cut_ROI_left_keV"],
				summary_data["cut_ROI_right_keV"]
			]
		)
		
		md"""
		### Selection Cuts
		
		$(cuts_df)
		"""
	end
end

# ╔═╡ 7cdc518c-cd96-4315-97ff-ab4ced327469
if length(summary_files) > 0 && selected_file != "No files available"
	begin
		# Create efficiency table
		eff_df = DataFrame(
			"Efficiency Component" => [
				"Contained Events",
				"Single Track",
				"Radial Cut",
				"Track Length Cut",
				"Z Fiducial Cut",
				"ROI Selection",
				"Two Electron ID",
				"**Total**"
			],
			"Value" => [
				round(summary_data["eff_contained"], digits=4),
				round(summary_data["eff_1trk"], digits=4),
				round(summary_data["eff_radial"], digits=4),
				round(summary_data["eff_trkl"], digits=4),
				round(summary_data["eff_zfid"], digits=4),
				round(summary_data["eff_roi"], digits=4),
				round(summary_data["eff_2e"], digits=4),
				round(summary_data["eff_total"], digits=8)
			],
			"Percentage" => [
				"$(round(100*summary_data["eff_contained"], digits=2))%",
				"$(round(100*summary_data["eff_1trk"], digits=2))%",
				"$(round(100*summary_data["eff_radial"], digits=2))%",
				"$(round(100*summary_data["eff_trkl"], digits=2))%",
				"$(round(100*summary_data["eff_zfid"], digits=2))%",
				"$(round(100*summary_data["eff_roi"], digits=2))%",
				"$(round(100*summary_data["eff_2e"], digits=2))%",
				"**$(round(100*summary_data["eff_total"], digits=6))%**"
			]
		)
		
		md"""
		### Efficiency Breakdown
		
		$(eff_df)
		"""
	end
end

# ╔═╡ 60d6bd0e-fa00-4312-9c7c-a859e2232e2b
md"""
## Compare All Summaries
"""

# ╔═╡ 0e156010-ce22-40b0-a6fd-e57f8ca503a3
if length(summary_files) > 1
      all_summaries = load_all_summaries(summary_dir, summary_files)
      comparison_df = create_comparison_table(all_summaries)

      md"""
      ### Efficiency Comparison Across All Files
      
      $(comparison_df)
      """
  else
      md"""
      Only one summary file available - no comparison possible.
      """
  end

# ╔═╡ 90bc6dd6-2904-4b16-9419-f9fc8147772c
comparison_df

# ╔═╡ f301f5cf-0dc6-49d9-bd4d-ff117ea56a2e
md"""
## Compute Backgrounds
"""

# ╔═╡ 8d10058a-ad7e-411a-898a-5910450ccd78
md"""
### Obtain selection efficiencies

"""

# ╔═╡ a1b2c3d4-5678-9abc-def0-123456789abc
if @isdefined(comparison_df)
	selection_eff = DataFrame(
		Component = comparison_df.DataType,
		sel_eff = comparison_df.Total
	)
	
	md"""
	### Selection Efficiencies
	
	$(selection_eff)
	"""
else
	md"No comparison_df available yet"
end

# ╔═╡ 7e81a3ed-4888-45f7-9d61-551aa11d64a4
md"""
### Read activity file

"""

# ╔═╡ 74c204fa-eb02-47b5-95fe-39da7e406bdb
hd5t_activities = read_activity_file()

# ╔═╡ 81d08904-154e-4234-a81b-eb7118838133
md"""
### Compute events per year: produced and selected. 

"""

# ╔═╡ a6cc62d9-2dff-452e-bac2-df5711097bc0
evt_year = nof_year(selection_eff, hd5t_activities)

# ╔═╡ b2c3d4e5-6789-abcd-ef01-23456789abcd
n_bi214_y_ton_roi = evt_year[13,"events_year_sel_bi214"]/4.5

# ╔═╡ c3d4e5f6-789a-bcde-f012-3456789abcde
n_tl208_y = evt_year[13,"events_year_sel_tl208"]/4.5

# ╔═╡ Cell order:
# ╠═04b446d6-f34f-11ed-2565-0b15d65b6781
# ╠═871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
# ╠═c17d78d1-39a1-40da-a3bc-d0a8717a4e10
# ╠═349825ff-7ffe-4fa1-ba26-a772041f0323
# ╠═c9fc0547-0e73-4629-9909-e59c3d75169d
# ╠═16a23b13-537b-45dd-b90f-7220e3d969d9
# ╠═eb8abbaf-6d7a-45e6-85ce-5ac5bc6dcfe7
# ╠═1629b0c1-64a3-4b6e-b24b-45f423f8a66b
# ╟─b97eb6e0-40d9-4de8-bc0b-c2b1469cabab
# ╠═3f2262a1-33e5-4fe2-a1ef-c03df3d1aa4b
# ╟─8a55b4a3-5cbf-48c3-b150-2bd4ad73f440
# ╠═9e2f3d4e-a5b6-4c7d-8e9f-1a2b3c4d5f6a
# ╠═f60557f4-113e-44ab-ab53-56e967f8fda8
# ╟─0fa67f2a-50a5-4c82-b42b-07f45f14e914
# ╟─7cdc518c-cd96-4315-97ff-ab4ced327469
# ╠═60d6bd0e-fa00-4312-9c7c-a859e2232e2b
# ╠═0e156010-ce22-40b0-a6fd-e57f8ca503a3
# ╠═90bc6dd6-2904-4b16-9419-f9fc8147772c
# ╠═f301f5cf-0dc6-49d9-bd4d-ff117ea56a2e
# ╠═8d10058a-ad7e-411a-898a-5910450ccd78
# ╠═a1b2c3d4-5678-9abc-def0-123456789abc
# ╠═7e81a3ed-4888-45f7-9d61-551aa11d64a4
# ╠═74c204fa-eb02-47b5-95fe-39da7e406bdb
# ╠═81d08904-154e-4234-a81b-eb7118838133
# ╠═a6cc62d9-2dff-452e-bac2-df5711097bc0
# ╠═b2c3d4e5-6789-abcd-ef01-23456789abcd
# ╠═c3d4e5f6-789a-bcde-f012-3456789abcde
