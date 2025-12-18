### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ d530f8a8-944e-4630-87d5-8d1c0cc587e7
begin
	using Markdown
	using CSV
	using DataFrames
	using Plots
	using Printf
	using HDF5
	using Statistics
	using StatsBase
end

# ╔═╡ edfdcde8-d9f5-11f0-a56f-97e2dd08f857
begin
	cmdir=joinpath(ENV["DATA"], "HD5t/itaca")
	pdir =joinpath(ENV["PROJECTS"], "Petit")
end

# ╔═╡ 9fb56c53-b3a0-4e74-8d90-e465b18583c2
begin
      using Pkg
      Pkg.activate(pdir)
      Pkg.instantiate()
  end

# ╔═╡ 6d9c26d2-3136-44fb-b3a6-1eceaac63ad4
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

# ╔═╡ 6d7a198c-e66b-4e2d-a8eb-85c917b04765
begin
	jn = ingredients(string(pdir,"/src/Petit.jl"))
	
end

# ╔═╡ e0e333f4-e857-41a0-85e5-3444fedd3f97
#it = ingredients(string(pdir,"/pluto/blob_analysis_functions.jl"))


# ╔═╡ 077cd0a7-f1f8-4a1c-a4f3-78b525c02c56
"""
    load_analysis_results(filepath::String) -> DataFrame

Load analysis_results.csv into a DataFrame.
"""
function load_analysis_results(filepath::String)
    CSV.read(filepath, DataFrame)
end

# ╔═╡ 423a098c-57e2-4f9c-99d1-8247669c1563
"""
    get_efficiencies_at_cut(eff_data, cut_value) -> (signal_eff, background_eff)

Get signal and background efficiencies at a given cut value by linear interpolation.
"""
function get_efficiencies_at_cut(eff_data, cut_value)
    cuts = eff_data.cuts
    
    # Find bracketing indices
    idx = searchsortedlast(cuts, cut_value)
    
    if idx == 0
        return (eff_data.signal_eff[1], eff_data.background_eff[1])
    elseif idx >= length(cuts)
        return (eff_data.signal_eff[end], eff_data.background_eff[end])
    end
    
    # Linear interpolation
    t = (cut_value - cuts[idx]) / (cuts[idx+1] - cuts[idx])
    sig_eff = eff_data.signal_eff[idx] + t * (eff_data.signal_eff[idx+1] - eff_data.signal_eff[idx])
    bkg_eff = eff_data.background_eff[idx] + t * (eff_data.background_eff[idx+1] - eff_data.background_eff[idx])
    
    return (signal_eff=sig_eff, background_eff=bkg_eff)
end

# ╔═╡ 6d359001-1f6a-43dc-a748-9d3753cbc7e3
function plot_vary(df::DataFrame, var;
                        title="Blob Energy Asymmetry",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50,
                        kwargs...)
    histogram(df[:, var];
              xlabel=xlabel,
              ylabel=ylabel,
              title=title,
              legend=false,
              bins=bins,
              kwargs...)
end

# ╔═╡ 2c2fd14f-39e5-4cff-ba66-1ab8549175a3
md"""
# Analysis
"""

# ╔═╡ f4a2a61e-c2c9-4592-a2d0-bd927a2afbdb
begin
	path_xe = "/Users/jjgomezcadenas/Data/HD5t/itaca/xe137/MST2_L100/itaca_analysis.csv"
	path_bb = "/Users/jjgomezcadenas/Data/HD5t/itaca/bb0nu/MST2_L100/itaca_analysis.csv"
	bbdf = load_analysis_results(path_bb)
	xedf = load_analysis_results(path_xe)
	path_xe_el = "/Users/jjgomezcadenas/Data/HD5t/itaca/xe137/MST_ELE_L100/itaca_analysis.csv"
	path_bb_el = "/Users/jjgomezcadenas/Data/HD5t/itaca/bb0nu/MST_ELE_L100/itaca_analysis.csv"
	bbeldf = load_analysis_results(path_bb_el)
	xeeldf = load_analysis_results(path_xe_el)
	path_xe_rb5 = "/Users/jjgomezcadenas/Data/HD5t/itaca/xe137/MST_ION_RB5_L100/itaca_analysis.csv"
	path_bb_rb5 = "/Users/jjgomezcadenas/Data/HD5t/itaca/bb0nu/MST_ION_RB5_L100/itaca_analysis.csv"
	bbrb5df = load_analysis_results(path_bb_rb5)
	xerb5df = load_analysis_results(path_xe_rb5)
end

# ╔═╡ effa65dd-b661-4757-b9a3-1081184c0c1d
md"""
## Pure Monte Carlo
"""

# ╔═╡ 32176bab-5df5-431e-9421-ca4102ff8261
md"""
### KDT method
"""

# ╔═╡ bcb8091a-6c3e-40e1-b71a-42aa73cb8b5c
let
	sbb =jn.Petit.plot_eb1_vs_eb2(bbdf;
                         Eb1="Eb1_mc_keV",
                         Eb2="Eb2_mc_keV",
                         title="Eb1 vs Eb2 bb0nu (MC/COMB)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	sxe =jn.Petit.plot_eb1_vs_eb2(xedf;
                         Eb1="Eb1_mc_keV",
                         Eb2="Eb2_mc_keV",
                         title="Eb1 vs Eb2 xe137 (MC/COMB)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	plot(sbb, sxe; 
         layout=(1,2), 
         size=(900, 400),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ fc04a87e-1f0b-4d34-bfbc-972dba3c1e04
jn.Petit.plot_eb2_cut_eff_and_fom(bbdf, xedf;
                                  Eb2="Eb2_mc_keV",
                                  length=100,
                                  rinf = 100,
                                  rsup = 800,
								  fom_bias=0.001,
                                  suptitle="MC/COMB")

# ╔═╡ ae765de2-1523-412a-a5d8-c375b954b6bf
begin
	eb2cut_mc_kdt = 350.0
	xeff_mc_kdt = jn.Petit.compute_efficiencies(bbdf[:, "Eb2_mc_keV"],
                              xedf[:, "Eb2_mc_keV"];
                              cuts = range(0, 801, length=81))
	eff_mc_kdt = get_efficiencies_at_cut(xeff_mc_kdt, eb2cut_mc_kdt)
	println("At cut=$(eb2cut_mc_kdt) keV: signal_eff=$(eff_mc_kdt.signal_eff), bkg_eff=$(eff_mc_kdt.background_eff)")
end

# ╔═╡ 84f13196-6015-49aa-8780-42c6e4fe8e34
let
	sbb =jn.Petit.plot_eb1_vs_eb2(bbrb5df;
                         Eb1="Eb1_mc_keV",
                         Eb2="Eb2_mc_keV",
                         title="Eb1 vs Eb2 bb0nu (MC/RB5)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	sxe =jn.Petit.plot_eb1_vs_eb2(xerb5df;
                         Eb1="Eb1_mc_keV",
                         Eb2="Eb2_mc_keV",
                         title="Eb1 vs Eb2 xe137 (MC/RB5)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	plot(sbb, sxe; 
         layout=(1,2), 
         size=(900, 400),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ 44f2c55a-4ca1-4126-b990-d067b55919ff
jn.Petit.plot_eb2_cut_eff_and_fom(bbrb5df, xerb5df;
                                  Eb2="Eb2_mc_keV",
                                  length=100,
                                  rinf = 100,
                                  rsup = 800,
								  fom_bias=0.001,
                                  suptitle="MC/RB5")

# ╔═╡ eea5716b-4491-4db6-ba6c-2f6c899b08fd
begin
	eb2cut_mc_rb5 = 220.0
	xeff_mc_rb5 = jn.Petit.compute_efficiencies(bbrb5df[:, "Eb2_mc_keV"],
                              xerb5df[:, "Eb2_mc_keV"];
                              cuts = range(0, 801, length=81))
	eff_mc_rb5 = get_efficiencies_at_cut(xeff_mc_rb5, eb2cut_mc_rb5)
	println("At cut=$(eb2cut_mc_rb5) keV: signal_eff=$(eff_mc_rb5.signal_eff), bkg_eff=$(eff_mc_rb5.background_eff)")
end

# ╔═╡ abc3cffa-2fb6-44de-aec4-648a44546b31
md"""
### MST method
"""

# ╔═╡ 188803cc-01b9-4155-84bd-a31519237c49
let
	sbb =jn.Petit.plot_eb1_vs_eb2(bbdf;
                         Eb1="Eb1_mc_mst_keV",
                         Eb2="Eb2_mc_mst_keV",
                         title="Eb1 vs Eb2 bb0nu (MC/MST)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	sxe =jn.Petit.plot_eb1_vs_eb2(xedf;
                         Eb1="Eb1_mc_mst_keV",
                         Eb2="Eb2_mc_mst_keV",
                         title="Eb1 vs Eb2 xe137 (MC/MST)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	plot(sbb, sxe; 
         layout=(1,2), 
         size=(900, 400),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ 5714aa54-9e5d-42a6-b08e-b1a4efefb88c
jn.Petit.plot_eb2_cut_eff_and_fom(bbdf, xedf;
                                  Eb2="Eb2_mc_mst_keV",
                                  length=100,
                                  rinf = 100,
                                  rsup = 800,
								  fom_bias=0.001,
                                  suptitle="MC/MST")

# ╔═╡ c1550a3b-3c0f-4394-9868-0797888dd8b6
begin
	eb2cut_mc_mst = 320.0
	xeff_mc_mst = jn.Petit.compute_efficiencies(bbdf[:, "Eb2_mc_mst_keV"],
                              xedf[:, "Eb2_mc_mst_keV"];
                              cuts = range(0, 801, length=81))
	eff_mc_mst = get_efficiencies_at_cut(xeff_mc_mst, eb2cut_mc_mst)
	println("At cut=$(eb2cut_mc_mst) keV: signal_eff=$(eff_mc_mst.signal_eff), bkg_eff=$(eff_mc_mst.background_eff)")
end

# ╔═╡ 5b0c30af-123d-402b-9aba-ac6953a42ca2
md"""
## Ions
"""

# ╔═╡ d478f35c-147a-4f5e-b273-495db6d600fa
md"""
### KDT method
"""

# ╔═╡ 426a0755-ef2e-4254-a653-4266f80f2160
let
	sbb =jn.Petit.plot_eb1_vs_eb2(bbdf;
                         Eb1="Eb1_keV",
                         Eb2="Eb2_keV",
                         title="Eb1 vs Eb2 bb0nu (RECO/COMB)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	sxe =jn.Petit.plot_eb1_vs_eb2(xedf;
                         Eb1="Eb1_keV",
                         Eb2="Eb2_keV",
                         title="Eb1 vs Eb2 xe137 (RECO/COMB)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	plot(sbb, sxe; 
         layout=(1,2), 
         size=(900, 400),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ 743c14ad-a4c8-49b2-ab12-50a73f0e636e
fom_kdt_ions = jn.Petit.plot_eb2_cut_eff_and_fom(bbdf, xedf;
                                  Eb2="Eb2_keV",
                                  length=100,
                                  rinf = 100,
                                  rsup = 800,
								  fom_bias=0.001,
                                  suptitle="RECO/IONS")

# ╔═╡ 56363146-177a-439a-ad6a-76432d007d2d
begin
	eb2cut_kdt = 305.0
	xeff_kdt = jn.Petit.compute_efficiencies(bbdf[:, "Eb2_keV"],
                              xedf[:, "Eb2_keV"];
                              cuts = range(0, 801, length=81))
	eff_kdt = get_efficiencies_at_cut(xeff_kdt, eb2cut_kdt)
	println("At cut=$(eb2cut_kdt) keV: signal_eff=$(eff_kdt.signal_eff), bkg_eff=$(eff_kdt.background_eff)")
end

# ╔═╡ 418d522d-9aa7-4d11-9771-7b7f69e6d1f2
md"""
### MST method
"""

# ╔═╡ 390e5934-c125-44fb-bcba-048ba8941939
let
	sbb =jn.Petit.plot_eb1_vs_eb2(bbdf;
                         Eb1="Eb1_mst_keV",
                         Eb2="Eb2_mst_keV",
                         title="Eb1 vs Eb2 bb0nu (RECO/MST)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	sxe =jn.Petit.plot_eb1_vs_eb2(xedf;
                         Eb1="Eb1_keV",
                         Eb2="Eb2_keV",
                         title="Eb1 vs Eb2 xe137 (RECO/MST)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	plot(sbb, sxe; 
         layout=(1,2), 
         size=(900, 400),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ ef3510ad-fb85-4177-92cf-820261279c26
jn.Petit.plot_eb2_cut_eff_and_fom(bbdf, xedf;
                                  Eb2="Eb2_mst_keV",
                                  length=100,
                                  rinf = 100,
                                  rsup = 800,
								  fom_bias=0.001,
                                  suptitle="RECO/COMB")

# ╔═╡ 23682cfc-b2c7-4609-b452-ab71d124c55e
begin
	eb2cut_mst = 332.0
	xeff_mst = jn.Petit.compute_efficiencies(bbdf[:, "Eb2_mst_keV"],
                              xedf[:, "Eb2_mst_keV"];
                              cuts = range(0, 801, length=81))
	eff_mst = get_efficiencies_at_cut(xeff_mst, eb2cut_mst)
	println("At cut=$(eb2cut_mst) keV: signal_eff=$(eff_mst.signal_eff), bkg_eff=$(eff_mst.background_eff)")
end

# ╔═╡ 4f38a6f1-3d17-4026-bacc-127934a8e495
let
	pabb = plot_vary(bbdf, "asymmetry_mc";
                        title="Asymmetry bb MC/COMB",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	paxe = plot_vary(xedf, "asymmetry_mc";
                        title="Asymmetry xe137 MC/COMB",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	pabbm = plot_vary(bbdf, "asymmetry_mc_mst";
                        title="Asymmetry bb MC/MST",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	paxem = plot_vary(xedf, "asymmetry_mc_mst";
                        title="Asymmetry xe137 MC/MST",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	plot(pabb, paxe, pabbm, paxem; 
         layout=(2,2), 
         size=(1000, 800),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ a6392fdf-42e7-4ddf-bafd-ea9280fb3b80
let
	pabb = plot_vary(bbdf, "asymmetry";
                        title="Asymmetry bb RECO/COMB",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	paxe = plot_vary(xedf, "asymmetry";
                        title="Asymmetry xe137 RECON/COMB",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	pabbm = plot_vary(bbdf, "asymmetry_mst";
                        title="Asymmetry bb RECO/MST",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	paxem = plot_vary(xedf, "asymmetry_mst";
                        title="Asymmetry xe137 RECO/MST",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	plot(pabb, paxe, pabbm, paxem; 
         layout=(2,2), 
         size=(1000, 800),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ 3aea06e0-b477-4ab7-8c98-27843a42bb13
md"""
## Electrons
"""

# ╔═╡ 269fc891-3684-45b4-b15f-009a5a6c57d9
md"""
### KDT method
"""

# ╔═╡ a051bfbb-463b-48ad-8606-bdc9cfb6e54f
let
	sbb =jn.Petit.plot_eb1_vs_eb2(bbeldf;
                         Eb1="Eb1_keV",
                         Eb2="Eb2_keV",
                         title="Eb1 vs Eb2 bb0nu (RECO/COMB)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	sxe =jn.Petit.plot_eb1_vs_eb2(xeeldf;
                         Eb1="Eb1_keV",
                         Eb2="Eb2_keV",
                         title="Eb1 vs Eb2 xe137 (RECO/COMB)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	plot(sbb, sxe; 
         layout=(1,2), 
         size=(900, 400),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ 6348f10b-c3dd-4850-8a71-b87a099ea187
fom_kdt_elec =jn.Petit.plot_eb2_cut_eff_and_fom(bbeldf, xeeldf;
                                  Eb2="Eb2_keV",
                                  length=100,
                                  rinf = 100,
                                  rsup = 800,
								  fom_bias=0.001,
                                  suptitle="RECO/ELEC")

# ╔═╡ ee87c4df-21cf-4e3c-acdc-459ac60c500f
begin
	eb2cut_el = 312.0
	xeff_el = jn.Petit.compute_efficiencies(bbeldf[:, "Eb2_keV"],
                              xeeldf[:, "Eb2_keV"];
                              cuts = range(0, 801, length=81))
	eff_el = get_efficiencies_at_cut(xeff_el, eb2cut_el)
	println("At cut=$(eb2cut_el) keV: signal_eff=$(eff_el.signal_eff), bkg_eff=$(eff_el.background_eff)")
end

# ╔═╡ 0438068c-7a04-43ae-845e-215e28003b76
md"""
### MST method
"""

# ╔═╡ 35b0cc37-addc-4cfa-af6b-d8cc7d7f1e1f
let
	sbb =jn.Petit.plot_eb1_vs_eb2(bbeldf;
                         Eb1="Eb1_mst_keV",
                         Eb2="Eb2_mst_keV",
                         title="Eb1 vs Eb2 bb0nu (RECO/MST)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	sxe =jn.Petit.plot_eb1_vs_eb2(xeeldf;
                         Eb1="Eb1_keV",
                         Eb2="Eb2_keV",
                         title="Eb1 vs Eb2 xe137 (RECO/MST)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	plot(sbb, sxe; 
         layout=(1,2), 
         size=(900, 400),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ e77a9cad-aea6-4103-a594-661f16c437e4
jn.Petit.plot_eb2_cut_eff_and_fom(bbeldf, xeeldf;
                                  Eb2="Eb2_mst_keV",
                                  length=100,
                                  rinf = 100,
                                  rsup = 800,
								  fom_bias=0.001,
                                  suptitle="RECO/COMB")

# ╔═╡ 8fc67da6-a287-4c02-bc7b-3aefd8676234
begin
	eb2cut_el_mst = 376.0
	xeff_el_mst = jn.Petit.compute_efficiencies(bbeldf[:, "Eb2_mst_keV"],
                              xeeldf[:, "Eb2_mst_keV"];
                              cuts = range(0, 801, length=81))
	eff_el_mst = get_efficiencies_at_cut(xeff_el_mst, eb2cut_el_mst)
	println("At cut=$(eb2cut_el_mst) keV: signal_eff=$(eff_el_mst.signal_eff), bkg_eff=$(eff_el_mst.background_eff)")
end

# ╔═╡ 65506270-0714-403e-963f-adf95c498f0e
let
	pabb = plot_vary(bbeldf, "asymmetry";
                        title="Asymmetry bb RECO/COMB",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	paxe = plot_vary(xeeldf, "asymmetry";
                        title="Asymmetry xe137 RECO/COMB",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	pabbm = plot_vary(bbeldf, "asymmetry_mst";
                        title="Asymmetry bb RECO/MST",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	paxem = plot_vary(xeeldf, "asymmetry_mst";
                        title="Asymmetry xe137 RECO/MST",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	plot(pabb, paxe, pabbm, paxem; 
         layout=(2,2), 
         size=(1000, 800),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ 323ab641-5108-4c8c-922e-73d07d3df524
md"""
## Comparisons: KDT method
"""

# ╔═╡ 1ff09b06-10dc-499d-8b6f-baf353bb9132
plot(fom_kdt_ions, fom_kdt_elec; 
         layout=(2,1), 
         size=(900, 800),
         margin=5Plots.mm,
         dpi=150)

# ╔═╡ 4dd6e818-0df1-4e60-8453-926ed420aebd
let
	pabb = plot_vary(bbdf, "asymmetry";
                        title="Asymmetry bb RECO/IONS",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	paxe = plot_vary(xedf, "asymmetry";
                        title="Asymmetry xe137 RECON/IONS",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	pabbm = plot_vary(bbeldf, "asymmetry";
                        title="Asymmetry bb RECO/ELEC",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	paxem = plot_vary(xeeldf, "asymmetry";
                        title="Asymmetry xe137 ELEC",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	plot(pabb, paxe, pabbm, paxem; 
         layout=(2,2), 
         size=(1000, 800),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ eedb27d1-6557-4604-a5d8-ee6375416344
begin
	bbdf2, xedf2 =jn.Petit.select_blob(bbdf, xedf; eb2cut=eb2cut_kdt)
	eff_bb = size(bbdf2)[1] / size(bbdf)[1]
	eff_xe = size(xedf2)[1] / size(xedf)[1]
	md"""
	- EB2 cut = $(eb2cut_kdt)
	- Selection efficiency for bb: $(eff_bb)
	- Selection efficiency for Xe137: $(eff_xe)
	"""
end

# ╔═╡ 0c66443b-cb05-4627-8312-ef6cb27818ec
begin
	bbeldf2, xeeldf2 =jn.Petit.select_blob(bbeldf, xeeldf; eb2cut=eb2cut_el)
	eff_bb_el = size(bbeldf2)[1] / size(bbeldf)[1]
	eff_xe_el = size(xeeldf2)[1] / size(xeeldf)[1]
	md"""
	- EB2 cut = $(eb2cut_el)
	- Selection efficiency for bb: $(eff_bb_el)
	- Selection efficiency for Xe137: $(eff_xe_el)
	"""
end

# ╔═╡ 27542af0-9acb-4e9f-b266-c42934376486
let
	pabb = plot_vary(bbdf2, "asymmetry";
                        title="Asymmetry bb RECO/IONS",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	paxe = plot_vary(xedf2, "asymmetry";
                        title="Asymmetry xe137 RECON/IONS",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	pabbm = plot_vary(bbeldf2, "asymmetry";
                        title="Asymmetry bb RECO/ELEC",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	paxem = plot_vary(xeeldf2, "asymmetry";
                        title="Asymmetry xe137 ELEC",
                        xlabel="Asymmetry (Eb1-Eb2)/(Eb1+Eb2)",
                        ylabel="Counts",
                        bins=50)
	plot(pabb, paxe, pabbm, paxem; 
         layout=(2,2), 
         size=(1000, 800),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ 5d7d922a-127d-456d-b98e-565178234d83
let
	sbb =jn.Petit.plot_eb1_vs_eb2(bbdf2;
                         Eb1="Eb1_keV",
                         Eb2="Eb2_keV",
                         title="Eb1 vs Eb2 bb0nu (RECO)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	sxe =jn.Petit.plot_eb1_vs_eb2(xedf2;
                         Eb1="Eb1_keV",
                         Eb2="Eb2_keV",
                         title="Eb1 vs Eb2 xe137 (RECO)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	sbb_mst =jn.Petit.plot_eb1_vs_eb2(bbdf2;
                         Eb1="Eb1_mst_keV",
                         Eb2="Eb2_mst_keV",
                         title="Eb1 vs Eb2 bb0nu (RECO/MST)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	sxe_mst =jn.Petit.plot_eb1_vs_eb2(xedf2;
                         Eb1="Eb1_keV",
                         Eb2="Eb2_keV",
                         title="Eb1 vs Eb2 xe137 (RECO/MST)",
                         xlabel="Eb1 (keV)",
                         ylabel="Eb2 (keV)")
	plot(sbb, sxe, sbb_mst, sxe_mst; 
         layout=(2,2), 
         size=(900, 800),
         margin=5Plots.mm,
         dpi=150)
end

# ╔═╡ Cell order:
# ╠═edfdcde8-d9f5-11f0-a56f-97e2dd08f857
# ╠═9fb56c53-b3a0-4e74-8d90-e465b18583c2
# ╠═6d9c26d2-3136-44fb-b3a6-1eceaac63ad4
# ╠═6d7a198c-e66b-4e2d-a8eb-85c917b04765
# ╠═e0e333f4-e857-41a0-85e5-3444fedd3f97
# ╠═d530f8a8-944e-4630-87d5-8d1c0cc587e7
# ╠═077cd0a7-f1f8-4a1c-a4f3-78b525c02c56
# ╠═423a098c-57e2-4f9c-99d1-8247669c1563
# ╠═6d359001-1f6a-43dc-a748-9d3753cbc7e3
# ╠═2c2fd14f-39e5-4cff-ba66-1ab8549175a3
# ╠═f4a2a61e-c2c9-4592-a2d0-bd927a2afbdb
# ╠═effa65dd-b661-4757-b9a3-1081184c0c1d
# ╠═32176bab-5df5-431e-9421-ca4102ff8261
# ╠═bcb8091a-6c3e-40e1-b71a-42aa73cb8b5c
# ╠═fc04a87e-1f0b-4d34-bfbc-972dba3c1e04
# ╠═ae765de2-1523-412a-a5d8-c375b954b6bf
# ╠═84f13196-6015-49aa-8780-42c6e4fe8e34
# ╠═44f2c55a-4ca1-4126-b990-d067b55919ff
# ╠═eea5716b-4491-4db6-ba6c-2f6c899b08fd
# ╠═abc3cffa-2fb6-44de-aec4-648a44546b31
# ╠═188803cc-01b9-4155-84bd-a31519237c49
# ╠═5714aa54-9e5d-42a6-b08e-b1a4efefb88c
# ╠═c1550a3b-3c0f-4394-9868-0797888dd8b6
# ╠═5b0c30af-123d-402b-9aba-ac6953a42ca2
# ╠═d478f35c-147a-4f5e-b273-495db6d600fa
# ╠═426a0755-ef2e-4254-a653-4266f80f2160
# ╠═743c14ad-a4c8-49b2-ab12-50a73f0e636e
# ╠═56363146-177a-439a-ad6a-76432d007d2d
# ╠═418d522d-9aa7-4d11-9771-7b7f69e6d1f2
# ╠═390e5934-c125-44fb-bcba-048ba8941939
# ╠═ef3510ad-fb85-4177-92cf-820261279c26
# ╠═23682cfc-b2c7-4609-b452-ab71d124c55e
# ╠═4f38a6f1-3d17-4026-bacc-127934a8e495
# ╠═a6392fdf-42e7-4ddf-bafd-ea9280fb3b80
# ╠═3aea06e0-b477-4ab7-8c98-27843a42bb13
# ╠═269fc891-3684-45b4-b15f-009a5a6c57d9
# ╠═a051bfbb-463b-48ad-8606-bdc9cfb6e54f
# ╠═6348f10b-c3dd-4850-8a71-b87a099ea187
# ╠═ee87c4df-21cf-4e3c-acdc-459ac60c500f
# ╠═0438068c-7a04-43ae-845e-215e28003b76
# ╠═35b0cc37-addc-4cfa-af6b-d8cc7d7f1e1f
# ╠═e77a9cad-aea6-4103-a594-661f16c437e4
# ╠═8fc67da6-a287-4c02-bc7b-3aefd8676234
# ╠═65506270-0714-403e-963f-adf95c498f0e
# ╠═323ab641-5108-4c8c-922e-73d07d3df524
# ╠═1ff09b06-10dc-499d-8b6f-baf353bb9132
# ╠═4dd6e818-0df1-4e60-8453-926ed420aebd
# ╠═eedb27d1-6557-4604-a5d8-ee6375416344
# ╠═0c66443b-cb05-4627-8312-ef6cb27818ec
# ╠═27542af0-9acb-4e9f-b266-c42934376486
# ╠═5d7d922a-127d-456d-b98e-565178234d83
