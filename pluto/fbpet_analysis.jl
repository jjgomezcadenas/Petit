### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 12057718-f493-11ed-2917-5d302ee966d9
using Pkg; Pkg.activate(ENV["JPetit"])

# ╔═╡ bd5bb73b-e9cc-4b90-bc42-46a21479583e
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Plots
	using Printf
	using Markdown
	using InteractiveUtils
	using Statistics
	using Chain
	using StatsBase
	using Distributions
	using StatsPlots
	using DataFramesMeta
	using HDF5
	import Glob
end

# ╔═╡ 362369f1-000a-47cd-822d-5ae5354aae40
PlutoUI.TableOfContents(title="FBPET analysis", indent=true)


# ╔═╡ 66742f20-d564-42aa-a674-aeef595b085a
md"""
# Analysis
"""

# ╔═╡ 60a77652-ad1c-4395-b448-0769a8dfb830
md"""
## Read the Data Frame and fix the data 
"""

# ╔═╡ b1c200eb-35cf-4a70-a89b-1bff382fcef2
cmdir="/Users/jjgomezcadenas/Data"

# ╔═╡ 2edfc61c-f3e5-4308-a438-eb60aa6b8f17
md"Read the DF? $(@bind readdf CheckBox(default=true))"

# ╔═╡ d0175ffb-0a4c-42fb-86b8-16dafde94f45
md""" 
select number of events for analysis: $(@bind ndf NumberField(10000:1000000, default=10000))"""

# ╔═╡ 9933fccf-ca15-4321-a6a9-b91514b08254
md"""
## Detector
"""

# ╔═╡ 95e939b4-af2d-4470-8e86-36cf210ccc68
md"""
### Data features
"""

# ╔═╡ 2b0e3527-9685-4d28-b213-18c7fe55f108
md"""
### Voxelize the detector
- We now divide the detector in boxes around phi and along z
"""

# ╔═╡ b85ab963-0abf-4317-a7d8-a12ed8d0ac07
md""" 
select transverse side of the boxes (in mm): $(@bind xyb NumberField(3.0:15.0, default=3.0))"""

# ╔═╡ 204620ea-98ef-4e5f-a7cb-61ea98697d5e
md"""
- Define voxelization functions in z and phi
- Define inverse voxelization function (give a position in z and phi from index)
"""

# ╔═╡ 1bf702e3-ff19-4df5-b656-13dac11a94b1
md""" 
Select energy cut: $(@bind ecut NumberField(300.0:450.0, default=350.0))"""

# ╔═╡ d6b46558-969b-43fd-9683-222a8385fd8c

md"Run the analysis? $(@bind arun CheckBox(default=false))"


# ╔═╡ dba902c6-de2c-46e2-8a90-fd895d39e9dc
if arun
md"""
- Photoelectric (0) vs Compton (1)
"""
end

# ╔═╡ b305aca5-7daa-4fa1-873a-6e32d2fdba9a
md"""
# Functions
"""

# ╔═╡ b4d6c180-05db-48b9-88a8-93271dbcfcdb
function getdirs(bdir::AbstractString)
	fdrs = Glob.glob("*", bdir)
	[split(f,"/")[end] for f in fdrs]
end

# ╔═╡ 36f22181-7cde-4556-95f1-9c7643931e2f
let
	readdir(cmdir)
	dirs = getdirs(cmdir)
	md""" Select set  : $(@bind sdata Select(dirs))"""
end

# ╔═╡ e00078f2-a4d8-437a-a832-2b48e77ec46a
let
	ddir    = joinpath(cmdir, sdata)
	dfiles = getdirs(ddir)
	md""" Select data  : $(@bind xdata Select(dfiles))"""
end

# ╔═╡ 5f2df8c0-5ef1-4791-ab91-cd5c5e28c670
begin
	xfile    = joinpath(cmdir, sdata, xdata)
end

# ╔═╡ e7ec0d87-6f27-4939-a3ae-6d67f9223efa
if occursin("LYSO", xfile)
	fwhm = 0.15
else
	fwhm=0.05
end

# ╔═╡ dbbc9c52-7ae3-41de-ad5f-97818965689c
function get_dfs(filename::String)

	fid     = h5open(filename, "r")
	dset =fid["MC"]["vertices"]
	dda = read(dset)
	ddic = Dict()
    for (i, key) in enumerate(keys(dda[1]))
        ddic[key] = [ddi[i] for ddi in dda]
    end
	close(fid)
    DataFrame(ddic)
end

# ╔═╡ b69d177a-e687-42d1-a200-80c2c3d0b047
if readdf
	df = get_dfs(xfile)
	df.event_id = Int.(df.event_id)
	df.parent_id = Int.(df.parent_id)
	df.process_id = Int.(df.process_id)
	df.track_id = Int.(df.track_id)
	df.volume_id = Int.(df.volume_id)
	df.post_KE   = Float64.(df.post_KE)
	df.pre_KE   = Float64.(df.pre_KE)
	df.t   = Float64.(df.t)
	df.x   = Float64.(df.x)
	df.y   = Float64.(df.y)
	df.z   = Float64.(df.z)
	@chain df begin
		@select!  $(Not(:deposited)) 
		@select!  $(Not(:parent_id)) 
		@select!  $(Not(:moved))
		@rsubset! :volume_id == 0 
		@rsubset! :process_id !=2 
		@transform! :edep = :pre_KE - :post_KE
		@rtransform! :phi = atan(:y,:x) 
		@rtransform! :rho = sqrt(:y^2 + :x^2) 
		@select!  $(Not(:pre_KE))
		@select!  $(Not(:post_KE))
		@select!  $(Not(:x))
		@select!  $(Not(:y))
		@select!  $(Not(:volume_id))
	end
	select!(df, :event_id, :process_id, :track_id, :z, :phi, :rho, :edep, :t)
		
end

# ╔═╡ 20076f64-6513-41f1-aacc-90d5a0a26f07
if readdf
	rhomin = round(minimum(df.rho), sigdigits=2)
	rhomax = round(maximum(df.rho), sigdigits=2)
	phimin = round(minimum(df.phi), sigdigits=4)
	phimax = round(maximum(df.phi), sigdigits=4)
	zmin = round(minimum(df.z), sigdigits=2)
	zmax = round(maximum(df.z), sigdigits=2)
	dr = rhomax - rhomin
	dz = zmax - zmin
	sigma = fwhm/2.3
	md"""
	- Detector (all dimensiones in mm )
	
	- (rmin, rmx) = ($rhomin, $rhomax)
	
	- (zmin, zmax) = ($zmin, $zmax)
	
	- (phimn, phimax) = ($phimin, $phimax)

	- Scintillator thickness = $dr

	- Scanner length = $dz

	- Energy resolution (FWHM) = $(100 * fwhm) %
	
	"""

end

# ╔═╡ 2b079a15-3270-42fa-8680-aecdad7700ac
if readdf
md"""
	- Box size (mm) = ($xyb, $xyb, $dr)
	"""
end

# ╔═╡ a137b9f5-806b-4468-9ecf-541e69113827
if readdf
	rl = 2 * π * rhomin
	nbr = Int(floor(rl / xyb))
	nbz = Int(floor(dz / xyb))
	nbox = nbz * nbr
	dphi = 2π / nbr
	dzb = dz / nbz
	md"""
	- number of boxes along radius = $nbr
	- number of boxes along z = $nbz
	- total number of boxes = $nbox
	- phi segment = $dphi
	- z segment = $dzb
	"""
end

# ╔═╡ 6df76b2d-df44-4972-bcad-26c4be83d5f4
function g_xxz_xxphi(xyb::Float64, dphi::Float64, zmax::Float64)
	function xxz(x::Vector{Float64}) 
		Int.(floor.((x .+ zmax) ./xyb) .+1)
	end
	function xxphi(x::Vector{Float64}) 
		Int.(floor.((x .+ π)./dphi) .+1)
	end
	xxz, xxphi
end

# ╔═╡ 435af28d-0917-4af4-b580-8e104d470914
function g_idxz_idxph(xyb::Float64, dphi::Float64, zmax::Float64)
	function zfri(iz::Vector{Int64})  
		iz * xyb  .- zmax .- 0.5 * xyb
	end
	function phifri(ip::Vector{Int64}) 
		ip * dphi  .- π .- 0.5 * dphi
	end
	zfri, phifri
end

# ╔═╡ d7626314-2e53-478e-b594-c6a40082f12c
if readdf
	zindx, phidx = g_xxz_xxphi(xyb, dphi, zmax)
	idxz, idxph = g_idxz_idxph(xyb, dphi, zmax)
end

# ╔═╡ 67679808-c562-4884-b985-da002cf2dcfe
function select_gammas(df, ndf::Int64, ecut::Float64, zindx, phidx)

	# 1. Select a subset of dataframe
	dfs = df[1:ndf, :]

	# 2. Compute box indexes for z and phi and add to DF
	dfs2 = @chain dfs begin 
    transform(_, [:z]  => zindx => :iz)
	transform(_, [:phi]  => phidx => :iphi)
	end

	3. # There must be one gamma per hemisphere
	dfc2 = combine(groupby(dfs2, :event_id)) do sdf
    n2 = length(unique(sdf.track_id)) 
  	n2 != 2 ? DataFrame() : DataFrame(event_id = sdf.event_id,
		                              process_id =sdf.process_id,
		                              track_id =sdf.track_id,
		                              z =sdf.z,
		                              iz = sdf.iz,
		                              phi= sdf.phi, 
		                              iphi = sdf.iphi,
		                              rho = sdf.rho,
		                              edep=sdf.edep,
									  t =sdf.t) 
	end

	# 4. Prepare the calculation of the R barycenter in the case that there are more than one gamma in a given box. 
	dfc = transform(dfc2, [:rho, :edep] => (.*) => :rhodep)

	# 5. Group gammas in boxes and add the energies of the gammas in the same box.
	# Also compute the barycenter
	dfbx = combine(groupby(dfc, :event_id)) do sdf
			combine(groupby(sdf, [:iphi, :iz])) do zpdf
				nb = nrow(zpdf) 
				ebox = sum(zpdf.edep)
				rbox = sum(zpdf.rhodep) / ebox
				DataFrame(event_id = zpdf.event_id, 
		                              track_id =zpdf.track_id,
		                              process_id =zpdf.process_id,
									  t =zpdf.t,
		                              z =zpdf.z, iz = zpdf.iz,
		                              phi= zpdf.phi, iphi = zpdf.iphi,
		    						  rho = zpdf.rho,
	                                  edep=zpdf.edep,
				ebox = ebox, rbox=rbox, nb=nb)
           end
		 
		end
	
	#6.  obtain z and phi (zb, phib) from the indexes (iz, iphi) 
	dfb = @chain dfbx begin 
	    transform(_, [:iz]  => idxz => :zb)
		transform(_, [:iphi]  => idxph => :phib)
	end

	#7. drop z, phi and rho (use the discretized values)
	dfx =@chain dfb begin
		@select  $(Not(:z)) 
		@select  $(Not(:phi)) 
		@select  $(Not(:rho))
	end

	#8. The next step is to treat the events in which there are more than one gamma in the box (they are duplicated and need to be cleaned up). It is only woth to consider nb=1 to 4

	dfy = @subset dfx :nb .< 5

	#8.1 In the case of nb =1 there is nothing to do.
	dfb1 = @subset dfy :nb .== 1

	#8.2 We consider also the cases with 2, 3 and 4 photons
	dfb2 = @subset dfy :nb .== 2
	dfb3 = @subset dfy :nb .== 3
	dfb4 = @subset dfy :nb .== 4

	#8.2.1 add an index to the dataframes
	dfb2[!,"Index"] = 1:size(dfb2)[1]
	dfb3[!,"Index"] = 1:size(dfb3)[1]
	dfb4[!,"Index"] = 1:size(dfb4)[1]

	#8.2.1 filter odd events for 2 photons (remove one of each two lines)
	dfb2f = filter(r -> isodd(r.Index), dfb2)

	#8.2.2 filter odd events, then odd events again for 3 photons
	# (remove two of each three lines)
	dfb3x = filter(r -> isodd(r.Index), dfb3)
	dfb3x[!,"Index"] = 1:size(dfb3x)[1]
	dfb3f = filter(r -> isodd(r.Index), dfb3x)

	#8.2.3 filter odd events, then odd (or even) events for 4 photons
	# (remove three of each four lines)
	dfb4x = filter(r -> isodd(r.Index), dfb4)
	dfb4x[!,"Index"] = 1:size(dfb4x)[1]
	dfb4f = filter(r -> isodd(r.Index), dfb4x)
	
	#8.2.4 drop the Index
	dfb2fx = @select dfb2f  $(Not(:Index))
	dfb3fx = @select dfb3f  $(Not(:Index))
	dfb4fx = @select dfb4f  $(Not(:Index))

	#8.3 vcat all the DF
	dff = reduce(vcat, [dfb1,dfb2fx,dfb3fx,dfb4fx]) 

	#9. drop nb and edep 
	dffx =@chain dff begin
		@select  $(Not(:nb))
		@select!  $(Not(:edep))
	end

	#10.0 Energy cut
	dfec = @subset dffx :ebox .> ecut

	#11.0 Impose that there is at least one gamma per hemisphere after energy cut
	dfec2 = combine(groupby(dfec, :event_id)) do sdf
	n2 = length(unique(sdf.track_id)) 
	n2 != 2 ? DataFrame() : DataFrame(event_id = sdf.event_id,
									  process_id =sdf.process_id,
									  track_id =sdf.track_id,
									  zb =sdf.zb,
									  iz = sdf.iz,
									  phib= sdf.phib, 
									  iphi = sdf.iphi,
									  rbox = sdf.rbox,
									  ebox=sdf.ebox,
									  t =sdf.t) 
	end
	dfec2
end

# ╔═╡ cc27d2c0-cca7-4513-a107-b6bba4f974f8
if arun
	dfs = select_gammas(df, ndf, ecut, zindx, phidx)
end

# ╔═╡ a089cd35-7e69-4e9c-b2e0-820f73fb5a2d
if arun
	xf = size(dfs)[1] / ndf
	md"""
	- Selected events fraction -> $(round(xf, sigdigits=4))
	"""
end

# ╔═╡ 1fb54cdb-12d9-4699-be3a-54ac388dede6
if arun
	dfp = @subset dfs :process_id .== 0
	xp = size(dfp)[1] / size(dfs)[1]
	md"""
	- Photoelectric fraction -> $(round(xp, sigdigits=4))
	"""
end
	

# ╔═╡ bf18c958-e124-4721-851e-ded2ec15e1c4
if arun
	dfep = @subset dfs :ebox .> 510.5 
	xe = size(dfep)[1] / size(dfs)[1]
	dfeg =@rsubset(dfs, :ebox > 511.0 - 2.0* sigma * 511.0, :ebox < 510.5)
	xe2 = size(dfeg)[1] / size(dfs)[1]
	stn = xe/xe2
	md"""
	- Events in peak -> $(round(xe, sigdigits=4))
	- Events in gaussian -> $(round(xe2, sigdigits=4))
	- S/N -> $(round(stn, sigdigits=4))
	"""
end

# ╔═╡ 0e24ca7c-1cd2-4b92-9663-c810456c413d
if arun
	ng_range = range(0,2, length=8)
	histogram(dfs.process_id, label="Process id ", bins=ng_range, color=:gray)
end

# ╔═╡ 6691e391-0372-48ee-a707-93cb30a233ac
if arun
	z_range = range(-1000.0 ,1000.0, length=100)
	histogram(dfs.zb, label=" z (mm)", bins=z_range, color=:gray)
end

# ╔═╡ 7a7bc24e-2364-4328-9ef2-dd820c04e145
if arun
	phi_range = range(-π ,π, length=100)
	histogram(dfs.phib, label=" phi (rad)", bins=phi_range, color=:gray)
end

# ╔═╡ acffb10a-ac1c-46ff-b2b4-da4134ced9d5
if arun
	r_range = range(350.0 ,400.0, length=100)
	histogram(dfs.rbox, label=" rho (mm) ", bins=r_range, color=:gray)
end

# ╔═╡ 69751ad1-56a3-4d13-a2e7-cf2b2db0fe21
if arun
	e_range = range(ecut ,550.0, length=100)
	histogram(dfs.ebox, label=" E (MeV) ", bins=e_range, legend=:topleft, color=:gray)
end

# ╔═╡ Cell order:
# ╠═12057718-f493-11ed-2917-5d302ee966d9
# ╠═bd5bb73b-e9cc-4b90-bc42-46a21479583e
# ╠═362369f1-000a-47cd-822d-5ae5354aae40
# ╠═66742f20-d564-42aa-a674-aeef595b085a
# ╠═60a77652-ad1c-4395-b448-0769a8dfb830
# ╠═b1c200eb-35cf-4a70-a89b-1bff382fcef2
# ╠═36f22181-7cde-4556-95f1-9c7643931e2f
# ╠═e00078f2-a4d8-437a-a832-2b48e77ec46a
# ╠═5f2df8c0-5ef1-4791-ab91-cd5c5e28c670
# ╠═2edfc61c-f3e5-4308-a438-eb60aa6b8f17
# ╠═b69d177a-e687-42d1-a200-80c2c3d0b047
# ╠═d0175ffb-0a4c-42fb-86b8-16dafde94f45
# ╠═9933fccf-ca15-4321-a6a9-b91514b08254
# ╠═95e939b4-af2d-4470-8e86-36cf210ccc68
# ╠═e7ec0d87-6f27-4939-a3ae-6d67f9223efa
# ╠═20076f64-6513-41f1-aacc-90d5a0a26f07
# ╠═2b0e3527-9685-4d28-b213-18c7fe55f108
# ╠═b85ab963-0abf-4317-a7d8-a12ed8d0ac07
# ╠═2b079a15-3270-42fa-8680-aecdad7700ac
# ╠═a137b9f5-806b-4468-9ecf-541e69113827
# ╠═204620ea-98ef-4e5f-a7cb-61ea98697d5e
# ╠═d7626314-2e53-478e-b594-c6a40082f12c
# ╠═1bf702e3-ff19-4df5-b656-13dac11a94b1
# ╠═d6b46558-969b-43fd-9683-222a8385fd8c
# ╠═cc27d2c0-cca7-4513-a107-b6bba4f974f8
# ╠═a089cd35-7e69-4e9c-b2e0-820f73fb5a2d
# ╠═dba902c6-de2c-46e2-8a90-fd895d39e9dc
# ╠═1fb54cdb-12d9-4699-be3a-54ac388dede6
# ╠═bf18c958-e124-4721-851e-ded2ec15e1c4
# ╠═0e24ca7c-1cd2-4b92-9663-c810456c413d
# ╠═6691e391-0372-48ee-a707-93cb30a233ac
# ╠═7a7bc24e-2364-4328-9ef2-dd820c04e145
# ╠═acffb10a-ac1c-46ff-b2b4-da4134ced9d5
# ╠═69751ad1-56a3-4d13-a2e7-cf2b2db0fe21
# ╠═b305aca5-7daa-4fa1-873a-6e32d2fdba9a
# ╠═67679808-c562-4884-b985-da002cf2dcfe
# ╠═b4d6c180-05db-48b9-88a8-93271dbcfcdb
# ╠═dbbc9c52-7ae3-41de-ad5f-97818965689c
# ╠═6df76b2d-df44-4972-bcad-26c4be83d5f4
# ╠═435af28d-0917-4af4-b580-8e104d470914
