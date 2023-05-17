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

# ╔═╡ d1d74a85-8c66-41c2-9bf0-1d2dc18945de
using Pkg; Pkg.activate(ENV["JPetit"])

# ╔═╡ 349825ff-7ffe-4fa1-ba26-a772041f0323
begin
	using PlutoUI
	using CSV
	using DataFrames
	#using Images
	#using ImageBinarization
	#using Colors
	using Plots
	using Printf
	#using Interpolations
	#using QuadGK
	using Markdown
	using InteractiveUtils
	#using LsqFit
	using Statistics
	using Chain
	using StatsBase
	using Distributions
	#using Unitful 
	using StatsPlots
	using DataFramesMeta
	#using UnitfulEquivalences 
	#using PhysicalConstants
	#using Peaks
	#using FFTW
	#using DSP
	#using Clustering
	using HDF5
	#using Roots
	#using BayesHistogram
	import Glob
end

# ╔═╡ 04b446d6-f34f-11ed-2565-0b15d65b6781
PlutoUI.TableOfContents(title="FBPET analysis", indent=true)


# ╔═╡ 6c59aeae-7990-4b43-8378-0de210a3291a
md"""
# Read the Data Frame and fix the data 
"""

# ╔═╡ 871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
cmdir="/Users/jjgomezcadenas/Data"

# ╔═╡ 5caa7c21-82e5-420b-b4e7-a0e33543b74b
md"Read the DF? $(@bind readdf CheckBox(default=true))"

# ╔═╡ 4c9a0060-3a5d-40ec-aefd-c08a88c77cbe
md""" 
select number of events for analysis: $(@bind ndf NumberField(10000:1000000, default=10000))"""

# ╔═╡ 587699d9-93d3-42dd-8bf7-3eb3db3b18e0
md"""
# Detector
"""

# ╔═╡ f77578be-5f37-4223-b47c-bc6c030a4238
md"""
## Data features
"""

# ╔═╡ 53c37edd-2cef-4743-9d9e-21e5518ddfeb
md"""
## Voxelize the detector
- We now divide the detector in boxes around phi and along z
"""

# ╔═╡ c83ce131-901f-451e-923f-18fc0f0386d4
md""" 
select transverse side of the boxes (in mm): $(@bind xyb NumberField(3.0:15.0, default=3.0))"""

# ╔═╡ d77923e2-662d-4bae-871b-4b7504136a94
md"""
- Define voxelization functions in z and phi
- Define inverse voxelization function (give a position in z and phi from index)
"""

# ╔═╡ 6d90338f-3833-4824-bc69-127441e9b94a
md"""
- indexes of zmin, zmax 
- Reverse zmin, zmax
"""

# ╔═╡ d781c59d-aca8-4437-a393-c45bf1bbdbd3
md"""
- indexes of phim (-π), phimax (π)
- Reverse phimin, phimax
"""

# ╔═╡ 7923ef0f-d663-408f-87ed-6606d5d77341
md"""
- Add box index to DF
"""

# ╔═╡ 8c35d719-c918-4b1d-b005-b914b6dfc177
md"""
# Selection

- The condition that trackid has a unique length of 2 reflects that track_id identifies the track number (1 or 2), and also, by definition the hemisphere. Thus a possible sequence is [1,1,1,2], which implies that track1 interacted four times and 2 only one, but the unique series is [1,2] and thus the length is 2 unless there is only one hemisphere, in which case the length of the unique series is 1
"""

# ╔═╡ b55205b1-46cc-4d95-bdad-aed3422f8d66
md"""
## Select events with at least one gamma per hemisphere

- The condition that trackid has a unique length of 2 reflects that track_id identifies the track number (1 or 2), and also, by definition the hemisphere. Thus a possible sequence is [1,1,1,2], which implies that track1 interacted four times and 2 only one, but the unique series is [1,2] and thus the length is 2 unless there is only one hemisphere, in which case the length of the unique series is 1
"""

# ╔═╡ 8f86481d-e9de-4f06-b097-980a5afd5e6b
md"Operate on the DF? $(@bind seltrk CheckBox(default=true))"

# ╔═╡ d2d18193-5d9a-4bae-8174-6aaaac733bbf
md"""
- The next transformation is needed to prepare the calculation of the R barycenter in the case that there are more than one gamma in a given box. 
"""

# ╔═╡ f638568e-4d84-429b-af03-6b57fec70af5
#@transform dfc begin
#    :rbar = :rho .* :edep 
#end

# ╔═╡ 51c29abb-a716-46d9-8ba5-a269c455cff3
md"""
## Group gammas in boxes
"""

# ╔═╡ 2ff47636-91d9-40de-92f5-e9a4363645d4
md"""
- The next combination groups first per event and then groups in terms of iphi and z, that is forms groups of tracks within an event that are in the same box. Then one computes the energy in the box (ebox) as the sum of the energies of the track in that box and the radius of the box (rbox) as the barycenter. Also we keep the number of tracks in the box
"""

# ╔═╡ babc3f53-0308-4b86-a3ce-d0cc4ebf7136
md"Operate on the DF? $(@bind selzphi CheckBox(default=true))"

# ╔═╡ 98d241e1-5695-4f31-b274-8081777b877d
md"""
- The next step is to obtain z and phi (zb, phib) from the indexes (iz, iphi) 
"""

# ╔═╡ c7dc5a77-6f9d-4e44-b8af-3c56e24ffd70
#begin
#	transform!(dfb, [:iz]  => idxz => :zb)
#	transform!(dfb, [:iphi]  => idxph => :phib)
#end

# ╔═╡ 4bd8ca38-e337-41f9-bd70-0a41d9a707b3
md"""
- Histograms show that box position is centered around true values
"""

# ╔═╡ ef2c6236-1033-46f9-b0b1-dfd6e4f3b4a0
md"""
- We can now drop phi, z and rho, keeping the "smeared values" of zb, phib and rbos
"""

# ╔═╡ cbc7a410-87fa-477b-aae5-61348df055d7
md"""
 ## Deal with events in boxes 
"""

# ╔═╡ da30c24f-d28c-4320-be8d-feed2361317c
md"""
- The next step is to treat the events in which there are more than one gamma in the box (they are duplicated and need to be cleaned up). It is only woth to consider nb=1 to 4
"""

# ╔═╡ a80f174a-401b-45b1-bf0a-a7ce050c8ced
md"""
- In the case of nb =1 there is nothing to do. We can histogram the difference between edep and ebox and show that both of them are identical. 
"""

# ╔═╡ e4dd52bc-dfa8-4bc8-b53a-db8f15e5afe4
md"""
- And we can histogram the energy spectrum
"""

# ╔═╡ 7a896f17-cf80-4214-97e6-feb91e3df79d
md"""
- We consider also the cases with 2, 3 and 4 photons
"""

# ╔═╡ 70be139b-e14f-4a67-a33e-52412915688e
md"""
- We add an index to the dataframes
"""

# ╔═╡ 9b496757-c3b2-4c09-bbf2-033e228ed344
md"""
- when nb = 2 we have duplicated rows, so we can drop the even rows
"""

# ╔═╡ c7d29501-b78d-45bd-a572-528c65f05353
md"""
In the case of 3 photons we filter two times the evens

a, a, a, b, b, b, c, c, c, d, d, d, e, e, e

1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15

- filter odd
a, a, b, c, c, d, e, e

1, 2, 3, 4, 5, 6, 7, 8

- filter odd again 

a, b, c, d, e...
"""

# ╔═╡ c3d8cf4c-9d71-4e52-bc63-f6d89f0e9622


# ╔═╡ 28a6e468-ad55-4790-8af5-dd338c65d3af


# ╔═╡ 1b4c91c9-d0dc-4f62-951d-2c431c1d71e8
md"""
- Finally we concat all the dataframes and drop nb and edep
"""

# ╔═╡ b3b132fa-af07-4c6a-aefd-5baa3bd8884d
md"""
## Energy cut
"""

# ╔═╡ cad85aa0-0c94-4d94-974f-ad1f86ecfbe8
md""" 
Select energy cut: $(@bind ecut NumberField(300.0:450.0, default=350.0))"""

# ╔═╡ abb82d8c-30bd-47d7-898f-ff09b8a4362b
md"""
## Impose that there is at least one gamma per hemisphere after energy cut
"""

# ╔═╡ d0f9e26b-2a3c-42a4-a826-7c8694f5d470
md"""
# Functions
"""

# ╔═╡ 3a0af47d-c297-4324-bfe8-a267b7cedb6d
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

# ╔═╡ 33b98f9c-0386-482b-a348-2060aa9da37d
function g_xxz_xxphi(xyb::Float64, dphi::Float64, zmax::Float64)
	function xxz(x::Vector{Float64}) 
		Int.(floor.((x .+ zmax) ./xyb) .+1)
	end
	function xxphi(x::Vector{Float64}) 
		Int.(floor.((x .+ π)./dphi) .+1)
	end
	xxz, xxphi
end

# ╔═╡ 942c2700-5576-4f59-829d-bcf038c8ab7a
function g_idxz_idxph(xyb::Float64, dphi::Float64, zmax::Float64)
	function zfri(iz::Vector{Int64})  
		iz * xyb  .- zmax .- 0.5 * xyb
	end
	function phifri(ip::Vector{Int64}) 
		ip * dphi  .- π .- 0.5 * dphi
	end
	zfri, phifri
end

# ╔═╡ 43c9fd8e-e94a-49d4-8237-d9b8f3a6a02c
function select_gammas1(df, nevt)
	# group by event
	evtgp = groupby(df, :event_id)
end

# ╔═╡ 70914efe-c290-4f43-8c2a-f521281c977e
function getdirs(bdir::AbstractString)
	fdrs = Glob.glob("*", bdir)
	[split(f,"/")[end] for f in fdrs]
end

# ╔═╡ ae89f5dc-e958-496a-91ac-0bd977355563
let
	readdir(cmdir)
	dirs = getdirs(cmdir)
	md""" Select set  : $(@bind sdata Select(dirs))"""
end

# ╔═╡ 9527f976-5e7f-413f-a474-a45320259d5e
begin
	ddir    = joinpath(cmdir, sdata)
	dfiles = getdirs(ddir)
	md""" Select data  : $(@bind xdata Select(dfiles))"""
end

# ╔═╡ 4af9a3ef-e883-4bc3-a2f1-212102e4951b
begin
	xfile    = joinpath(cmdir, sdata, xdata)
end

# ╔═╡ 37a23226-4e3c-4e21-a4ee-a9aef75b2093
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

# ╔═╡ 99b43220-0ce1-41e2-a1f5-447387c51bdf
dfs = df[1:ndf, :]

# ╔═╡ 78648d07-cddc-4276-afa5-6d58b9cbc48e
begin
	rhomin = round(minimum(df.rho), sigdigits=2)
	rhomax = round(maximum(df.rho), sigdigits=2)
	phimin = round(minimum(df.phi), sigdigits=4)
	phimax = round(maximum(df.phi), sigdigits=4)
	zmin = round(minimum(df.z), sigdigits=2)
	zmax = round(maximum(df.z), sigdigits=2)
	dr = rhomax - rhomin
	dz = zmax - zmin
	md"""
	- Detector (all dimensiones in mm )
	
	- (rmin, rmx) = ($rhomin, $rhomax)
	
	- (zmin, zmax) = ($zmin, $zmax)
	
	- (phimn, phimax) = ($phimin, $phimax)

	- Scintillator thickness = $dr

	- Scanner length = $dz
	
	
	"""

end

# ╔═╡ 0eb16588-4eee-4634-885f-ea82e3b40fef
md"""
	- Box size (mm) = ($xyb, $xyb, $dr)
	"""

# ╔═╡ fbe9062c-548e-4b5d-8b98-2325d69ed8ec
begin
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

# ╔═╡ c78a4526-59a3-48f7-9dc3-9fb4ff1e4131
begin
	zindx, phidx = g_xxz_xxphi(xyb, dphi, zmax)
	idxz, idxph = g_idxz_idxph(xyb, dphi, zmax)
end

# ╔═╡ ddd7466e-8c0b-4a22-bf9e-33c97e213dd6
idxz([0,667])

# ╔═╡ abdb8508-bd96-42b9-84fb-a5be97cc18e2
phidx([-π,π])

# ╔═╡ 38d339fe-09ed-488c-8036-945a8a25c577
idxph([1,734])

# ╔═╡ 76e89264-4b4b-4df2-a290-813294b091e8
dfs2 = @chain dfs begin 
    transform(_, [:z]  => zindx => :iz)
	transform(_, [:phi]  => phidx => :iphi)
end

# ╔═╡ b0273b4a-7339-4021-b58c-10cbac518baa
if seltrk
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
end

# ╔═╡ bd766b2d-0c43-4a5a-ab26-765bfdbe7540
dfc = transform(dfc2, [:rho, :edep] => (.*) => :rhodep)

# ╔═╡ 2f78c866-afcd-480f-9bb1-b644a5ee8530
begin
	rf1 = round(size(dfc)[1]/size(dfs)[1], sigdigits=2)
md"""
- reduction factor imposing two gammas in DF = $(rf1)
"""
end

# ╔═╡ 3839fe4a-dc59-4725-bf5f-64157319d267
if selzphi
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
end

# ╔═╡ 5ff0ae52-26f2-4c10-9c5a-ee2ff39ca60f
dfb = @chain dfbx begin 
    transform(_, [:iz]  => idxz => :zb)
	transform(_, [:iphi]  => idxph => :phib)
end

# ╔═╡ 57c38187-18ad-4598-9aab-dbb42a44a2c7
begin
	phi_range = range(-0.01, 0.01, length=20)
	histogram(dfb.phi - dfb.phib, label="phi - phi box", bins=phi_range, color=:gray)
end

# ╔═╡ 70b6b320-7a1b-46f9-ad67-846c282f5fbd
begin
	z_range = range(-5, 5.0, length=20)
	histogram(dfb.z - dfb.zb, label="z - z box", bins=z_range, color=:gray)
end

# ╔═╡ feb2a48f-8f0c-4785-907c-a8e1fa7876ac
begin
	r_range = range(-2.0, 2.0, length=20)
	histogram(dfb.rho - dfb.rbox, label="r - r box", bins=r_range, color=:gray)
end

# ╔═╡ 7012f684-5bfd-44fd-9630-40a74a27d325
dfx =@chain dfb begin
		@select  $(Not(:z)) 
		@select  $(Not(:phi)) 
		@select  $(Not(:rho))
	end

# ╔═╡ d085bbcc-7ca6-48b2-aa4e-7502044db28c
begin
	nb_range = range(0, 10, length=20)
	histogram(dfx.nb, label="nb", bins=nb_range, color=:gray)
end

# ╔═╡ c7ddda3a-ec97-4c87-9f04-6346879485e3
dfy = @subset dfx :nb .< 5

# ╔═╡ d3654b89-6af0-43a1-a515-b2919d9d7c65
dfb1 = @subset dfy :nb .== 1

# ╔═╡ a2148a73-0647-4c61-8de0-2a80447dc97e
dfb2 = @subset dfy :nb .== 2

# ╔═╡ 3fa015c9-e21f-453d-94c0-e0f5af2312ef
dfb2f = filter(r -> isodd(r.Index), dfb2)

# ╔═╡ 1d68a5c7-b87f-4ee4-bb40-37ef6d19718a
dfb3 = @subset dfy :nb .== 3

# ╔═╡ 7be15005-cf88-4c92-a967-426d87a3d0d5
dfb3x = filter(r -> isodd(r.Index), dfb3)

# ╔═╡ 90addf2e-a7c0-4324-8b36-db16e86a5a49
dfb3x[!,"Index"] = 1:size(dfb3x)[1]

# ╔═╡ 4bb732f5-5e09-4b3d-9ff8-fd9d611d6b77
dfb3f = filter(r -> iseven(r.Index), dfb3x)

# ╔═╡ 6366a088-d40f-4897-8887-c6dd8ee333bd
dfb4 = @subset dfy :nb .== 4

# ╔═╡ 81e5e3d3-97f3-4165-ac07-230c2909c99e
begin
	dfb2[!,"Index"] = 1:size(dfb2)[1]
	dfb3[!,"Index"] = 1:size(dfb3)[1]
	dfb4[!,"Index"] = 1:size(dfb4)[1]
end

# ╔═╡ 5a911953-ca2b-422a-8792-6c950d7d2028
begin
	dfb4x = filter(r -> isodd(r.Index), dfb4)
	dfb4x[!,"Index"] = 1:size(dfb4x)[1]
	dfb4f = filter(r -> iseven(r.Index), dfb4x)
end

# ╔═╡ 98437a6a-c685-4600-b5a0-9b29be3b462c
begin
	dfb2fx = @select dfb2f  $(Not(:Index))
	dfb3fx = @select dfb3f  $(Not(:Index))
	dfb4fx = @select dfb4f  $(Not(:Index))
end

# ╔═╡ ec7440d3-0011-4f9f-82b7-9303771d9037
dff = reduce(vcat, [dfb1,dfb2fx,dfb3fx,dfb4fx]) 

# ╔═╡ d0b571ce-330c-4aea-8091-4c558eec3dd6
dffx =@chain dff begin
	@select  $(Not(:nb))
	@select!  $(Not(:edep))
end

# ╔═╡ 3e850fa3-b9f8-4e8d-971d-42ea820b824b
dfec = @subset dffx :ebox .> ecut

# ╔═╡ 4d765fbc-139a-4abb-994a-29a2e6e0d31d
begin
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
end

# ╔═╡ f22e8a10-75db-4512-b923-cf65e6988d93
ngdf = combine(groupby(dfec2, :event_id)) do sdf
	combine(groupby(sdf, :track_id)) do tdf
		ng = nrow(tdf)
	end
end

# ╔═╡ 802e8e17-6ffd-422b-975f-21d52d7287af
begin
	ng_range = range(0,10., length=20)
	histogram(ngdf.x1, label="ng", bins=ng_range, color=:gray)
end

# ╔═╡ cde0a154-9a93-4b3a-bcae-b6fff44c0e5d
begin
	e_range = range(-5.0, 5.0, length=20)
	histogram(dfb.edep - dfb.ebox, label="E - E box", bins=e_range, color=:gray)
end

# ╔═╡ 700eeea5-9632-4e77-87e4-8c20f9dec2d6
begin
	es_range = range(0.0, 550., length=50)
	histogram(dfb.ebox, label="Ebox", bins=es_range, color=:gray)
end

# ╔═╡ 477aecc3-b21e-454c-bb1f-f34c828455bf
histogram(dfb2f.ebox, label="Ebox", bins=es_range, color=:gray)

# ╔═╡ 30389e7d-aa0f-45cd-8424-7f7662518091
histogram(dfb3f.ebox, label="Ebox", bins=es_range, color=:gray)

# ╔═╡ e45cd68b-defc-4478-a311-b88dd5713c68
histogram(dfb4f.ebox, label="Ebox", bins=es_range, color=:gray)

# ╔═╡ 59226b94-f112-4840-be05-8c72493423ba
histogram(dffx.ebox, label="Ebox", bins=es_range, legend=:topleft, color=:gray)

# ╔═╡ 6ef4055d-1eb8-4f6a-9775-49a9a9bf365f
histogram(dfec.ebox, label="Ebox", bins=es_range, legend=:topleft, color=:gray)

# ╔═╡ 5602db59-3161-436b-bfcc-64aa7b9fb05e
histogram(dfec2.ebox, label="Ebox", bins=es_range, legend=:topleft, color=:gray)

# ╔═╡ beb6e5ca-67c9-406e-81b1-770dceff69ab
zindx([zmin,zmax])

# ╔═╡ cd3cc0bd-b5f8-402a-9ffd-9d4e9b032d06
function select_gammas(df, nevt)
	# group by event
	evtgp = groupby(df, :event_id)
	

	EVTID = zeros(Int64, 0)
	IH    = zeros(Int64, 0)
	IPHI  = zeros(Int64, 0)
	IZ    = zeros(Int64, 0)
	RB    = zeros(Float64, 0)
	EB    = zeros(Float64, 0)
	
	for ievt::Int64 in 1:length(evtgp)
		dfe = evtgp[ievt]  
		length(unique(evtgp[3].track_id))
		
		# group by hemisphere
		ihgp = groupby(dfe, :trak_id) 
		nh = length(ihgp)
		#println("nh =",nh)
		if nh != 2
			#println("not two hemispheres, ignore event")
			continue
		end

		for ih = 1:2
			dfh = ihgp[ih]
			if size(dfh)[1] == 1
				push!(EVTID, dfh[1,"event_id"])
				push!(IH,    dfh[1,"ih"])
				push!(IPHI,  dfh[1,"iphi"])
				push!(IZ,    dfh[1,"iz"])
				push!(RB,    dfh[1,"rho"])
				push!(EB,    dfh[1,"deposited"])
			else
				# group by iphi, iz
				boxgp = groupby(dfh, [:iphi, :iz])
				nbx = size(boxgp)[1]
				for ib in 1:nbx
					dfb = boxgp[ib]
					if size(dfb)[1] == 1
						push!(EVTID, dfb[1,"event_id"])
						push!(IH,    dfb[1,"ih"])
						push!(IPHI,  dfb[1,"iphi"])
						push!(IZ,    dfb[1,"iz"])
						push!(RB,    dfb[1,"rho"])
						push!(EB,    dfb[1,"deposited"])
					else
						dfbx = transform(dfb, [:deposited]  => sum => :ebox)
						dfgs = @chain dfbx begin
        					@rtransform  begin
           						:rbar = :rho * :deposited / :ebox
       						end
							@combine :rbox = sum(:rbar)
						end
						push!(EVTID, dfbx[1,"event_id"])
						push!(IH,    dfbx[1,"ih"])
						push!(IPHI,  dfbx[1,"iphi"])
						push!(IZ,    dfbx[1,"iz"])
						push!(RB,    dfgs[1,"rbox"])
						push!(EB,    dfbx[1,"deposited"])
					end	
				end
			end
		end	
	end

# ╔═╡ Cell order:
# ╠═04b446d6-f34f-11ed-2565-0b15d65b6781
# ╠═d1d74a85-8c66-41c2-9bf0-1d2dc18945de
# ╠═349825ff-7ffe-4fa1-ba26-a772041f0323
# ╠═6c59aeae-7990-4b43-8378-0de210a3291a
# ╠═871bd8bf-8e4b-40fb-a9c7-fdeb47589c5a
# ╠═ae89f5dc-e958-496a-91ac-0bd977355563
# ╠═9527f976-5e7f-413f-a474-a45320259d5e
# ╠═4af9a3ef-e883-4bc3-a2f1-212102e4951b
# ╠═5caa7c21-82e5-420b-b4e7-a0e33543b74b
# ╠═37a23226-4e3c-4e21-a4ee-a9aef75b2093
# ╠═4c9a0060-3a5d-40ec-aefd-c08a88c77cbe
# ╠═99b43220-0ce1-41e2-a1f5-447387c51bdf
# ╠═587699d9-93d3-42dd-8bf7-3eb3db3b18e0
# ╠═f77578be-5f37-4223-b47c-bc6c030a4238
# ╠═78648d07-cddc-4276-afa5-6d58b9cbc48e
# ╠═53c37edd-2cef-4743-9d9e-21e5518ddfeb
# ╠═c83ce131-901f-451e-923f-18fc0f0386d4
# ╠═0eb16588-4eee-4634-885f-ea82e3b40fef
# ╠═fbe9062c-548e-4b5d-8b98-2325d69ed8ec
# ╠═d77923e2-662d-4bae-871b-4b7504136a94
# ╠═c78a4526-59a3-48f7-9dc3-9fb4ff1e4131
# ╠═6d90338f-3833-4824-bc69-127441e9b94a
# ╠═beb6e5ca-67c9-406e-81b1-770dceff69ab
# ╠═ddd7466e-8c0b-4a22-bf9e-33c97e213dd6
# ╠═d781c59d-aca8-4437-a393-c45bf1bbdbd3
# ╠═abdb8508-bd96-42b9-84fb-a5be97cc18e2
# ╠═38d339fe-09ed-488c-8036-945a8a25c577
# ╠═7923ef0f-d663-408f-87ed-6606d5d77341
# ╠═76e89264-4b4b-4df2-a290-813294b091e8
# ╠═8c35d719-c918-4b1d-b005-b914b6dfc177
# ╠═b55205b1-46cc-4d95-bdad-aed3422f8d66
# ╠═8f86481d-e9de-4f06-b097-980a5afd5e6b
# ╠═b0273b4a-7339-4021-b58c-10cbac518baa
# ╠═d2d18193-5d9a-4bae-8174-6aaaac733bbf
# ╠═bd766b2d-0c43-4a5a-ab26-765bfdbe7540
# ╠═f638568e-4d84-429b-af03-6b57fec70af5
# ╠═2f78c866-afcd-480f-9bb1-b644a5ee8530
# ╠═51c29abb-a716-46d9-8ba5-a269c455cff3
# ╠═2ff47636-91d9-40de-92f5-e9a4363645d4
# ╠═babc3f53-0308-4b86-a3ce-d0cc4ebf7136
# ╠═3839fe4a-dc59-4725-bf5f-64157319d267
# ╠═98d241e1-5695-4f31-b274-8081777b877d
# ╠═c7dc5a77-6f9d-4e44-b8af-3c56e24ffd70
# ╠═5ff0ae52-26f2-4c10-9c5a-ee2ff39ca60f
# ╠═4bd8ca38-e337-41f9-bd70-0a41d9a707b3
# ╠═57c38187-18ad-4598-9aab-dbb42a44a2c7
# ╠═70b6b320-7a1b-46f9-ad67-846c282f5fbd
# ╠═feb2a48f-8f0c-4785-907c-a8e1fa7876ac
# ╠═ef2c6236-1033-46f9-b0b1-dfd6e4f3b4a0
# ╠═7012f684-5bfd-44fd-9630-40a74a27d325
# ╠═d085bbcc-7ca6-48b2-aa4e-7502044db28c
# ╠═cbc7a410-87fa-477b-aae5-61348df055d7
# ╠═da30c24f-d28c-4320-be8d-feed2361317c
# ╠═c7ddda3a-ec97-4c87-9f04-6346879485e3
# ╠═a80f174a-401b-45b1-bf0a-a7ce050c8ced
# ╠═d3654b89-6af0-43a1-a515-b2919d9d7c65
# ╠═cde0a154-9a93-4b3a-bcae-b6fff44c0e5d
# ╠═e4dd52bc-dfa8-4bc8-b53a-db8f15e5afe4
# ╠═700eeea5-9632-4e77-87e4-8c20f9dec2d6
# ╠═7a896f17-cf80-4214-97e6-feb91e3df79d
# ╠═a2148a73-0647-4c61-8de0-2a80447dc97e
# ╠═1d68a5c7-b87f-4ee4-bb40-37ef6d19718a
# ╠═6366a088-d40f-4897-8887-c6dd8ee333bd
# ╠═70be139b-e14f-4a67-a33e-52412915688e
# ╠═81e5e3d3-97f3-4165-ac07-230c2909c99e
# ╠═9b496757-c3b2-4c09-bbf2-033e228ed344
# ╠═3fa015c9-e21f-453d-94c0-e0f5af2312ef
# ╠═477aecc3-b21e-454c-bb1f-f34c828455bf
# ╠═c7d29501-b78d-45bd-a572-528c65f05353
# ╠═7be15005-cf88-4c92-a967-426d87a3d0d5
# ╠═90addf2e-a7c0-4324-8b36-db16e86a5a49
# ╠═c3d8cf4c-9d71-4e52-bc63-f6d89f0e9622
# ╠═4bb732f5-5e09-4b3d-9ff8-fd9d611d6b77
# ╠═30389e7d-aa0f-45cd-8424-7f7662518091
# ╠═5a911953-ca2b-422a-8792-6c950d7d2028
# ╠═e45cd68b-defc-4478-a311-b88dd5713c68
# ╠═98437a6a-c685-4600-b5a0-9b29be3b462c
# ╠═28a6e468-ad55-4790-8af5-dd338c65d3af
# ╠═1b4c91c9-d0dc-4f62-951d-2c431c1d71e8
# ╠═ec7440d3-0011-4f9f-82b7-9303771d9037
# ╠═d0b571ce-330c-4aea-8091-4c558eec3dd6
# ╠═59226b94-f112-4840-be05-8c72493423ba
# ╠═b3b132fa-af07-4c6a-aefd-5baa3bd8884d
# ╠═cad85aa0-0c94-4d94-974f-ad1f86ecfbe8
# ╠═3e850fa3-b9f8-4e8d-971d-42ea820b824b
# ╠═6ef4055d-1eb8-4f6a-9775-49a9a9bf365f
# ╠═abb82d8c-30bd-47d7-898f-ff09b8a4362b
# ╠═4d765fbc-139a-4abb-994a-29a2e6e0d31d
# ╠═5602db59-3161-436b-bfcc-64aa7b9fb05e
# ╠═f22e8a10-75db-4512-b923-cf65e6988d93
# ╠═802e8e17-6ffd-422b-975f-21d52d7287af
# ╠═d0f9e26b-2a3c-42a4-a826-7c8694f5d470
# ╠═3a0af47d-c297-4324-bfe8-a267b7cedb6d
# ╠═33b98f9c-0386-482b-a348-2060aa9da37d
# ╠═942c2700-5576-4f59-829d-bcf038c8ab7a
# ╠═43c9fd8e-e94a-49d4-8237-d9b8f3a6a02c
# ╠═70914efe-c290-4f43-8c2a-f521281c977e
# ╠═cd3cc0bd-b5f8-402a-9ffd-9d4e9b032d06
