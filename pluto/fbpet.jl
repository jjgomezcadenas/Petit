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

# ╔═╡ 805955a4-f018-4e68-a68a-1e6b64a99378
using Pkg; Pkg.activate(ENV["JPetit"])

# ╔═╡ 53dee408-58f8-11ed-2b44-db5e37f3a91e
begin
	using PlutoUI
	using CSV
	using DataFrames
	#using Images
	#using ImageBinarization
	#using Colors
	using Plots
	using Printf
	using Interpolations
	using QuadGK
	using Markdown
	using InteractiveUtils
	using LsqFit
	using Statistics
	using Chain
	using StatsBase
	using Distributions
	using Unitful 
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


# ╔═╡ a69771ca-f756-4f84-8bb4-1c044fa4db1f
using ATools

# ╔═╡ f75fa3dd-d3dc-4f5d-877b-cc9d8f57d95c
PlutoUI.TableOfContents(title="PETIT analysis", indent=true)

# ╔═╡ 921fb297-2c3d-4f45-847c-b2f009052e9f
md"""
# Analysis
"""

# ╔═╡ 918dc9ec-f027-4d68-8290-96f7f56a416a
md"""
## Read the data
"""

# ╔═╡ 3274825d-61a0-4856-9340-97420d237c27
begin
	dfiles ="*.h5"
	xpath = "/Users/jjgomezcadenas/Data"
	xfiles = Glob.glob(dfiles, xpath)
end

# ╔═╡ aebd4a68-087f-4ea7-8865-72254727a550
md"""
## Load a DF with the MC data
"""

# ╔═╡ 96fb0461-3e42-4030-9394-79ed63c9f423
md"""
Interactions: 0 = Compton, 1 = Photo, 2 = Rayl
"""

# ╔═╡ 442ff64b-9569-4718-b4b2-d9e811479e1a
md"""
## Prepare dataframe
"""

# ╔═╡ 3cb04bd1-9ad0-4f82-861f-44d9f01dad15
#begin
#	dfv0 = filter(:volume_id => x -> x == 0, df)
#	df2 = filter(:process_id => n->n!=2, dfv0)
#	select!(df2, Not(:volume_id))
#	df2
#end

# ╔═╡ 87a35667-62a2-429c-a244-31f964931647
md"""
- Events which interact in the detector (volume_id = 0)
- Eliminate Rayleigh 
"""

# ╔═╡ 06f1185b-e43d-48dc-97b3-f7041324f483
md"""
### Compute cylindrical coordinates
"""

# ╔═╡ 992fca6f-ada2-4be5-8217-95d4323d45fd
md"""
### Detector
"""

# ╔═╡ bd496019-51e2-46ee-96e7-048bf53dc547
md""" 
select transverse side of the box (in mm): $(@bind xyb NumberField(3.0:15.0, default=3.0))"""

# ╔═╡ a4336af9-dc53-4f36-87cf-7a5f3d757907
#mbox = zeros(Int8, 2, 3)

# ╔═╡ a314bbe4-72fe-4b01-8084-5b7f99cd6e6c
xyb

# ╔═╡ 32f728f2-7617-431e-963d-0c2ba88fdfd6
function g_xxz_xxphi(xyb::Float64, dphi::Float64, zmax::Float64)
	function xxz(x::Vector{Float64}) 
		Int.(floor.((x .+ zmax) ./xyb))
	end
	function xxphi(x::Vector{Float64}) 
		Int.(floor.((x .+ π)./dphi))
	end
	xxz, xxphi
end

# ╔═╡ d8a2d1fd-8b20-42b7-9893-bfd7e34d623b
function g_idxz_idxph(xyb::Float64, dphi::Float64, zmax::Float64)
	function zfri(iz::Vector{Int64})  
		iz * xyb  .- zmax 
	end
	function phifri(ip::Vector{Int64}) 
		ip * dphi  .- π
	end
	zfri, phifri
end

# ╔═╡ 3e0582a0-efdf-4f58-a250-c14d5b26090d
md"""
### Event selection
"""

# ╔═╡ f848b71e-341a-4443-aeb9-b4319a850e24
md"""
### First subset: pure photoelectric interaction, no losses.

- This corresponds to events that have pre_KE = 511.0 and deposited = 511.0. In these, the incoming gamma didn't loose energy before interacting with the scintillator and deposits all its energy in a single step via photoelectric interaction. 
"""

# ╔═╡ 72161bfb-2e27-40e3-979d-c3ef9d802c46
#dfpreke = filter(:pre_KE => x -> x == 511.0, df2)

# ╔═╡ d29f2a93-0439-4e8c-bf85-7a619e1be725
#@rsubset df2 :pre_KE == 511.0

# ╔═╡ 9b2e95aa-c85c-469b-96d8-e669b236a17b
#dfphot = filter(:deposited => x -> x == 511.0, dfpreke)

# ╔═╡ 0928e121-3ea5-4ebd-9d12-63102a794793
md"""
- **groupby** *event_id* and **combine** in terms of *nrow*. The nrow column gives the number of gammas interacting in the scintillator. If only one, the second gamma has escaped detection and the event cannot be used to make a LOR and must be filtered away
"""

# ╔═╡ 16a272f5-cd53-4df2-8ec6-cb9a4768665a
y = DataFrame(id=rand('a':'d', 100), v=rand(100))

# ╔═╡ a00aac30-bf3f-49be-b3ba-c7165c894e3f
#grouped_df  = groupby(dfphot, :event_id)

# ╔═╡ cb9cddb0-668a-47ad-9f73-70eeb37c90d0
typeof(gpphot)

# ╔═╡ a783d217-7861-4532-8774-f7d2e0f3e86b
#combine(gpphot, y -> size(y)[1] < 2 ? DataFrame() : y)

# ╔═╡ 81543561-f9eb-4dc9-bbb3-2f498c66d1a7
gpphot

# ╔═╡ 9f5563d2-64aa-402a-b64e-431e4b2466a1
#begin
#grouped_df  = groupby(x, :id)
#
#end

# ╔═╡ 0c56f110-dcf5-4507-9459-00aad5f35b97
combine(groupby(x, :id)) do sdf
    n = size(sdf)[1]
    n < 25 ? DataFrame() : DataFrame(n=n) # drop groups with low number of rows
end

# ╔═╡ 294832c6-7ffa-43e7-8f46-5cbf4833d784
combine(groupby(x, :id)) do sdf
    n = nrow(sdf)
    n < 25 ? DataFrame() : DataFrame(n=n) # drop groups with low number of rows
end

# ╔═╡ 82fdad6f-0c38-483e-9247-a51ccfe4a31b
DataFrame()

# ╔═╡ 7a0bca86-f120-49ee-b2ed-7daee1508636
md"""
- histogram shows the proportion of 1 vs 2 photons
"""

# ╔═╡ f8402e1a-1b4b-4573-a591-ddda211cc0a4
begin
	hp2nrow, pp2nrow = hist1d(ceventidgp.nrow, "number of gammas", 10, 0, 10)
	plot(pp2nrow)
end

# ╔═╡ 3c29ebab-2be4-4c5e-9fbe-a038f6dde3cc
md"""
- **filter** events with nrow =1
"""

# ╔═╡ 9f3825d7-081d-4fb3-863d-ff1f845b66ac
begin
	n1ceventidgp = filter(:nrow => x -> x== 1, ceventidgp)
	l1p = length(n1ceventidgp.nrow)
	md"""
	Number of photoelectric events with only 1 gamma =$(l1p)
	"""
end

# ╔═╡ 2376caa4-3b68-4b52-aab4-d5fa786b525a
md"""
- **filter** events with nrow =2
"""

# ╔═╡ e5c1feba-b371-4b9f-9585-406e484ccdc2
begin
	n2ceventidgp = filter(:nrow => x -> x== 2, ceventidgp)
	l2p = length(n2ceventidgp.nrow)
	md"""
	Number of photoelectric events with 2 gamma =$(l2p)
	"""
end


# ╔═╡ 6da39c1e-3885-4bca-923f-2d47cf7b48d7
n2evts = n2ceventidgp.event_id;

# ╔═╡ 2c90ed9c-3128-4535-8aaf-04a9660a3607
md"""
- Uncomment to select events in which both gammas interact exactly by photoelectric
(very slow)
"""

# ╔═╡ 86051cb1-fe25-4029-bb24-2610cd7fe646
#dfphot[in.(dfphot[!,"event_id"], (n2evts,)), :]

# ╔═╡ a5797d58-845d-4358-80c0-f4d0a93db650
#dfn2phot = dfphot[in(n2evts).(dfphot.event_id), :]

# ╔═╡ 0aa8c55c-c931-43e9-bf7a-f9576efac54e
md"""
### Extended subset: at least two gammas with sum of post_KE above cut

- This corresponds to events that may be used for real reconstruction. 
"""

# ╔═╡ 965c8192-1e24-4c9b-bc71-df8c9c7a26ed
md"""
#### Selection algorithm
- Groupby events
- Separate each event in two hemispheres
- Add the energy of the photons in a box
- Compute R as the baricenter of the Ri (where i is gammas in the box)
- In each box: (E -> Esum, R - > Barycenter Ri, z->zbox, phi->phibox)
"""

# ╔═╡ 8975b62c-9ec9-4a9c-a9ec-5a4d4ea51196
md"""
- Select and event
"""

# ╔═╡ 9d82fd3d-4e90-43ed-99fc-35930c05b66b
function assign_hemisphere(df)
	x1 = df.x[1]
	y1 = df.y[1]
	z1 = df.z[1]
	nrm1 = sqrt(x1^2 + y1^2 + z1^2)
	nrmi =[sqrt(df.x[i]^2 + df.y[i]^2 + df.z[i]^2) for i in 1:length(df.x)]
	pdot = [df.x[i] * x1 + df.y[i] * y1 + df.z[i] * z1 for i in 1:length(df.x)]
	pdotn = pdot./(nrm1*nrmi)
	[p > 0 ? 1 : -1 for p in pdotn]
	
end

# ╔═╡ 88343dd9-d306-4a53-b28d-60086a8d70a9
ievt = 3

# ╔═╡ 7d539e78-5756-4dad-8491-18bbfde24915


# ╔═╡ 3f29923c-3d21-44c5-940c-13328277b4d3
md"""
- Group by hemisphere
"""

# ╔═╡ f78bdd04-c4f6-480a-960a-9219abd6fb13
md"""
- Group based on phi, z
"""

# ╔═╡ d64d7c14-b6a1-4219-aaef-2f182d3ff0a4
md"""
- For groups with more than one photon, compute rbar y esum
"""

# ╔═╡ e68c51e4-e78a-4994-b19d-ec2048437054
aaa = zeros(Int64, 0)

# ╔═╡ 7b58081f-a2e0-477c-a313-f89453dd81a6
push!(aaa,1)

# ╔═╡ 8fadf7ab-9a2c-49ea-89e3-d64aea1baf81
function select_gammas(df, nevt)
	# group by event
	evtgp = groupby(df, :event_id)
	nevtt = size(evtgp)

	EVTID = zeros(Int64, 0)
	IH    = zeros(Int64, 0)
	IPHI  = zeros(Int64, 0)
	IZ    = zeros(Int64, 0)
	RB    = zeros(Float64, 0)
	EB    = zeros(Float64, 0)
	for ievt in 1:nevt
		dfe = evtgp[ievt]  
		dfe[!, "ih"] = assign_hemisphere(dfe) 

		# group by hemisphere
		ihgp = groupby(dfe, :ih) 
		nh = size(ihgp)[1]
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
	dctg = Dict("event_id" => EVTID, 
		"ebox" => EB,
		"rbox" => RB,
		"iphi" => IPHI,
	    "iz"   => IZ,
	    "ih"   => IH)
	
	DataFrame(dctg)
end
	

# ╔═╡ 5b8c66ec-a038-42d8-be20-91281cb52ee4
function

# ╔═╡ 42550f71-0d65-4569-8332-41774a979e67
md"""
## Control plots
"""

# ╔═╡ 7fa8480c-ac4e-4dec-be26-52b49a343b9f
#pid0df.post_KE

# ╔═╡ 99eb99db-56f7-4884-bb1a-ca362fee1ba0
md""" Plane 2: select rmax histogram: $(@bind rmax NumberField(400.0:600.0, default=500.0))"""

# ╔═╡ f43ddffc-c0da-4140-9ddd-5668708129a4
md""" Plane 2: select width of peak for fit: $(@bind wpeak NumberField(5.0:20.0, default=10.0))"""

# ╔═╡ 5fc2c3b0-3cdf-4bdd-9f83-408a68a78f97


# ╔═╡ 8efc9eb5-603e-4799-92f0-46a30de30dc3


# ╔═╡ c0816fb9-00c7-41c4-b2e4-da4dc50fdbd1
md"""
# Functions
"""

# ╔═╡ 88aab6f6-07a8-42c4-b39f-33ed7bb0a75b
"""
Read h5 files with PETIT data, return a dataframe

"""
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

# ╔═╡ 082a3d2a-c915-4c98-9221-d4610ae112bc
begin
	df = get_dfs("/Users/jjgomezcadenas/Data/MC-000.h5")
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
	select!(df, Not(:parent_id))
	select!(df, Not(:moved))
	df[!, "deposited"] = df.pre_KE - df.post_KE
	df
end

# ╔═╡ 7d52ddc9-3c15-4358-ad70-1dbba22a30dc
names(df)

# ╔═╡ 55a2ca09-b190-4214-b943-c18e87aea87a
begin
	df2 = @chain df begin
		@rsubset :volume_id == 0 
		@rsubset :process_id !=2   
	end
	select!(df2, Not(:volume_id))
end

# ╔═╡ 2f1351aa-b954-4494-bb83-06038cb29549
dfc = @chain df2 begin
    transform(_, [:x, :y] => ((x,y) -> atan.(y,x)) => :phi)
	transform(_, [:x, :y] => ((x,y) -> sqrt.(x.^2+y.^2)) => :rho)
end

# ╔═╡ cd03a379-5d2b-4633-881c-77e9f10d2c31
dfl = length(unique(dfc.event_id))

# ╔═╡ 44c411e6-dda6-45be-980e-0b807847bd4e
md"""
- number of events in DF = $dfl
"""

# ╔═╡ f92670e3-6ac2-4f78-a115-2c124f3141a3
begin
	rhomin = round(minimum(dfc.rho), sigdigits=2)
	rhomax = round(maximum(dfc.rho), sigdigits=2)
	phimin = round(minimum(dfc.phi), sigdigits=4)
	phimax = round(maximum(dfc.phi), sigdigits=4)
	zmin = round(minimum(dfc.z), sigdigits=2)
	zmax = round(maximum(dfc.z), sigdigits=2)
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

# ╔═╡ df9e2fad-886e-40a7-a0d6-832b81b84661
begin
	#xyb = 9.0
	zb = dr
	md"""
	Box size (mm) = ($xyb, $xyb, $dr)
	"""
end

# ╔═╡ bd38bc17-7fcf-4767-b1e3-852b606fa939
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

# ╔═╡ aaed8c95-8062-4993-9507-581c1a9b788d
zindx, phidx = g_xxz_xxphi(xyb, dphi, zmax)

# ╔═╡ bbe4047d-a5da-4548-b0a4-61bee9bcc82d
zindx([-1000.0,1000.0])

# ╔═╡ 11150995-4664-418a-877c-e56199ba9087
phidx([-π,π])

# ╔═╡ be59ff78-a14b-4451-bf40-c188c36fe6e7
idxz, idxph = g_idxz_idxph(xyb, dphi, zmax)

# ╔═╡ ffcb8e5d-bf77-43a5-a399-2c5fbf349b7e
idxz([1,0])

# ╔═╡ a3f1ba28-54a3-4cc5-8191-8d2ebb2e68f2
dfc

# ╔═╡ f69968ff-587d-489d-bb4d-348a2ee1e26d
dfcp = transform(dfc, [:z]  => zindx => :iz)

# ╔═╡ 95346821-547a-4743-a302-4dc2bc0b9149
histogram(dfcp.iz)

# ╔═╡ b03afad1-4d61-48c1-850d-737889cf5f9c
dfcz = transform(dfcp, [:phi]  => phidx => :iphi)

# ╔═╡ f4ed7867-9668-4060-89f4-578edaba7347
histogram(dfcz.iphi)

# ╔═╡ e555c4ec-dac1-4d8c-a18c-daf801a13b77
maximum(dfcz.iphi)

# ╔═╡ 125ba2a4-a345-4146-b4c6-2a5075526842
histogram(idxz(dfcz.iz))

# ╔═╡ ade2e2a1-a86c-4147-8db3-32a8fa14890f
idxz(dfcz.iz)

# ╔═╡ ef37af9e-5fde-4c71-9735-c85f39d5df4c
typeof(dfcz.iz)

# ╔═╡ 7607b7e3-5e49-4ac9-b51d-392e66463bad
histogram(idxph(dfcz.iz))

# ╔═╡ 4919adc9-8c8c-4324-b5de-abdd92e64172
evtgp = groupby(dfcz, :event_id)

# ╔═╡ 03526651-e04b-43db-874f-b4a80757bbf2
size(evtgp)

# ╔═╡ 83a12404-352a-4ab1-96da-160e666c4cec
dfe = evtgp[ievt]

# ╔═╡ 64ab2bd2-da2e-4aae-b271-ad4463e640e7
dfe[!, "ih"] = assign_hemisphere(dfe)

# ╔═╡ e4bbaf3c-43bc-423b-ae5c-c72c62a95926
ihgp = groupby(dfe, :ih)

# ╔═╡ e369d39c-f283-40df-971b-b81cc9429e95
ihgp[1]

# ╔═╡ aa4d9b4c-499c-45c8-8812-a2a3339bcbf1
size(ihgp[1])

# ╔═╡ 18f041fe-c89c-45ac-b25a-3e8b5ccc0e8e
ihgp[2]

# ╔═╡ c4202fbe-7f45-4b7a-a627-b5749072386b
boxgp = groupby(ihgp[1], [:iphi, :iz])

# ╔═╡ 978be772-d6e2-4009-bfe2-e9ce7ee12248
boxgp[1]

# ╔═╡ a9924d52-9674-4cd9-a0c3-26c6115a819b
boxgp[2]

# ╔═╡ 72ec35ac-992e-4cfb-b7bd-74239ae860ee
boxgp[3]

# ╔═╡ 823f2c63-a3dd-44ec-b0dd-e3a8e60f147c
dfbx = transform(boxgp[2], [:deposited]  => sum => :ebox)

# ╔═╡ ee9884d4-0256-43de-b20e-0d4a4fa2c7f1
dfgs = @chain dfbx begin
        @rtransform  begin
           :rbar = :rho * :deposited / :ebox
       end
	@combine :rbox = sum(:rbar)
end


# ╔═╡ dbda5eb8-3183-4ab9-8906-a7a3524f8fa5
begin
	dctg = Dict("event_id" => dfbx[1,"event_id"], 
		"ebox" => dfbx[1,"ebox"],
		"rbox" => dfgs[1,"rbox"],
		"iphi" => dfbx[1,"iphi"],
	    "iz"   => dfbx[1,"iz"],
	    "ih"   => dfbx[1,"ih"])
	
	dfg = DataFrame(dctg)
end

# ╔═╡ fce159d5-e6e7-4ca4-9f6b-5b5d9be6e085
dfgg = @time select_gammas(dfcz, 2000)

# ╔═╡ 9beee121-1392-48ae-9ce4-cd881d0d4958
histogram(dfgg.iphi)

# ╔═╡ 8b4aff6a-8c6e-42dc-9da7-f1215a2c6d92
histogram(dfgg.ebox)

# ╔═╡ 991a69af-5163-4721-9dc5-53965dba6508
histogram(dfgg.rbox)

# ╔═╡ 954ce38b-f89f-411b-9621-6e6319c99123
begin
evtdfg = groupby(dfgg, :event_id)
evtdfgnr = combine(evtdfg, nrow)
end

# ╔═╡ 362d7ec1-a33a-480a-a3cd-353d2d882e96
histogram(evtdfgnr.nrow)

# ╔═╡ e4b5d37b-6090-408f-8087-5e90b7de4e1a
dfphot = @rsubset(df2, :pre_KE == 511.0, :deposited == 511.0)

# ╔═╡ a1bfb274-a912-4d86-af1f-c0a0a8bfed24
dfphotl = length(unique(dfphot.event_id))

# ╔═╡ df043f80-9873-4959-a6ca-f88ab90f4336
md"""
- Select events in which pre_KE = 511 keV (thus no losses before scintillator), and deposited=511 keV (thus all energy deposited by photoelectric) reduction factor =$(dfphotl / dfl)
"""

# ╔═╡ 7c993118-9c81-4fa3-a439-e8d898bf107d
md"""
 - Fraction of photoelectric events with 2 gamma interacting wrt all photoelectric =$(l2p/dfphotl)
 -  Fraction of photoelectric events with 2 gamma interacting wrt all events =$(l2p/dfl)
"""

# ╔═╡ 7fe29d2e-4fcf-4884-b379-086107a96010
begin
grouped_df = groupby(dfphot, :event_id)
#grouped_df  = groupby(y, :id)
combine(grouped_df, x -> size(x)[1] < 2 ? DataFrame() : x)
end

# ╔═╡ 6cdd8029-d3ed-461e-aa7d-1087ff72eba2
combine(grouped_df, x -> size(x)[1] == 1 ? DataFrame() : x)

# ╔═╡ a9b68546-cee8-4f6b-bd26-6fd6f1ca3bcd
typeof(grouped_df)

# ╔═╡ 16a70b43-3f82-432e-ad53-6b2a6985669d
combine(grouped_df, x -> size(x)[1] < 2 ? DataFrame() : x)

# ╔═╡ eb0e11bd-5418-4923-a094-53222aee6ce7
combine(groupby(dfphot, :event_id)) do sdf
    n = size(sdf)[1]
	DataFrame()
end

# ╔═╡ 0f316ae1-e6be-446f-81d9-56c1a64fe988
@chain dfphot begin
    groupby(:event_id)
    transform(:deposited => mean)
end

# ╔═╡ 52f6a0e2-a5f6-452e-ba45-14958b9ee3a4
begin
	hpke, ppke = hist1d(df.post_KE, "post_KE", 50, 100., 600.0)
	hprke, pprke = hist1d(df.pre_KE, "pre_KE", 50, 100., 600.0)
	plot(ppke, pprke)
end

# ╔═╡ b7793747-d3ee-49eb-a88a-68f730358931
begin
	hx, px = hist1d(df.x, "x", 50, -400., 400.0)
	hy, py = hist1d(df.y, "y", 50, -400., 400.0)
	hz, pz = hist1d(df.z, "z", 50, -400., 400.0)
	plot(px, py, pz)
end

# ╔═╡ 519399b4-8f8c-4fcb-bfe3-4f7d7d67e761
begin
	hxy, pxy = hist2d(df.x,df.y, 100,
                "X", "Y", -500., 500.0, -500., 500.0; title="X-Y")
	plot(pxy)
end

# ╔═╡ 9093a05f-ecfd-497b-95b1-7530b5a49541
begin
	hxz, pxz = hist2d(df.x,df.z, 100,
                "X", "z", -500., 500.0, -500., 500.0; title="X-z")
	plot(pxz)
end

# ╔═╡ ee8468ab-2987-497c-83ac-720902561490
"""
Read files and output a concatenated dataframes

"""
function concdf(xfiles::Vector{String})
	df = get_dfs(xfiles[1])
	for file in xfiles[2:end]
		dfx = get_dfs(file)
		df = vcat(df, dfx)
	end
	df
end

# ╔═╡ c0ba6116-ae04-4ab5-a5c8-08d623444a3e
"""
Interpolate a function to vectors (x, y)

"""
function v2f(x::Vector{Float64}, y::Vector{Float64}, overflow::Float64=0.0)
	function gfpdf_(fi, xmin::Real, xmax::Real, bkgnd::Real=0.0)
		function fn(x)
		   	if x < xmin || x > xmax
				return bkgnd
			else
				return fi(x) 
			end
		end
		return fn
	end
    li = LinearInterpolation(x, y)
	gfpdf_(li, x[1], x[end], overflow)
end

# ╔═╡ 2d18834e-70f2-47f0-9238-026d1d6fe462
function rfromft(fg::ATools.FGauss)
	R = 2.355 * fg.std[1]/fg.mu[1]
	δR = R * fg.δstd[1]/fg.std[1]  # mu is negligible
	(var=R, delta=δR)
end

# ╔═╡ Cell order:
# ╠═805955a4-f018-4e68-a68a-1e6b64a99378
# ╠═53dee408-58f8-11ed-2b44-db5e37f3a91e
# ╠═a69771ca-f756-4f84-8bb4-1c044fa4db1f
# ╠═f75fa3dd-d3dc-4f5d-877b-cc9d8f57d95c
# ╠═921fb297-2c3d-4f45-847c-b2f009052e9f
# ╠═918dc9ec-f027-4d68-8290-96f7f56a416a
# ╠═3274825d-61a0-4856-9340-97420d237c27
# ╠═aebd4a68-087f-4ea7-8865-72254727a550
# ╠═082a3d2a-c915-4c98-9221-d4610ae112bc
# ╠═7d52ddc9-3c15-4358-ad70-1dbba22a30dc
# ╠═96fb0461-3e42-4030-9394-79ed63c9f423
# ╠═442ff64b-9569-4718-b4b2-d9e811479e1a
# ╠═3cb04bd1-9ad0-4f82-861f-44d9f01dad15
# ╠═87a35667-62a2-429c-a244-31f964931647
# ╠═55a2ca09-b190-4214-b943-c18e87aea87a
# ╠═06f1185b-e43d-48dc-97b3-f7041324f483
# ╠═2f1351aa-b954-4494-bb83-06038cb29549
# ╠═cd03a379-5d2b-4633-881c-77e9f10d2c31
# ╠═44c411e6-dda6-45be-980e-0b807847bd4e
# ╠═992fca6f-ada2-4be5-8217-95d4323d45fd
# ╠═f92670e3-6ac2-4f78-a115-2c124f3141a3
# ╠═bd496019-51e2-46ee-96e7-048bf53dc547
# ╠═df9e2fad-886e-40a7-a0d6-832b81b84661
# ╠═bd38bc17-7fcf-4767-b1e3-852b606fa939
# ╠═a4336af9-dc53-4f36-87cf-7a5f3d757907
# ╠═a314bbe4-72fe-4b01-8084-5b7f99cd6e6c
# ╠═32f728f2-7617-431e-963d-0c2ba88fdfd6
# ╠═d8a2d1fd-8b20-42b7-9893-bfd7e34d623b
# ╠═aaed8c95-8062-4993-9507-581c1a9b788d
# ╠═be59ff78-a14b-4451-bf40-c188c36fe6e7
# ╠═bbe4047d-a5da-4548-b0a4-61bee9bcc82d
# ╠═11150995-4664-418a-877c-e56199ba9087
# ╠═a3f1ba28-54a3-4cc5-8191-8d2ebb2e68f2
# ╠═f69968ff-587d-489d-bb4d-348a2ee1e26d
# ╠═95346821-547a-4743-a302-4dc2bc0b9149
# ╠═b03afad1-4d61-48c1-850d-737889cf5f9c
# ╠═f4ed7867-9668-4060-89f4-578edaba7347
# ╠═e555c4ec-dac1-4d8c-a18c-daf801a13b77
# ╠═125ba2a4-a345-4146-b4c6-2a5075526842
# ╠═ffcb8e5d-bf77-43a5-a399-2c5fbf349b7e
# ╠═ade2e2a1-a86c-4147-8db3-32a8fa14890f
# ╠═ef37af9e-5fde-4c71-9735-c85f39d5df4c
# ╠═7607b7e3-5e49-4ac9-b51d-392e66463bad
# ╠═3e0582a0-efdf-4f58-a250-c14d5b26090d
# ╠═f848b71e-341a-4443-aeb9-b4319a850e24
# ╠═72161bfb-2e27-40e3-979d-c3ef9d802c46
# ╠═d29f2a93-0439-4e8c-bf85-7a619e1be725
# ╠═9b2e95aa-c85c-469b-96d8-e669b236a17b
# ╠═e4b5d37b-6090-408f-8087-5e90b7de4e1a
# ╠═a1bfb274-a912-4d86-af1f-c0a0a8bfed24
# ╠═df043f80-9873-4959-a6ca-f88ab90f4336
# ╠═0928e121-3ea5-4ebd-9d12-63102a794793
# ╠═16a272f5-cd53-4df2-8ec6-cb9a4768665a
# ╠═7fe29d2e-4fcf-4884-b379-086107a96010
# ╠═a00aac30-bf3f-49be-b3ba-c7165c894e3f
# ╠═cb9cddb0-668a-47ad-9f73-70eeb37c90d0
# ╠═a783d217-7861-4532-8774-f7d2e0f3e86b
# ╠═6cdd8029-d3ed-461e-aa7d-1087ff72eba2
# ╠═81543561-f9eb-4dc9-bbb3-2f498c66d1a7
# ╠═9f5563d2-64aa-402a-b64e-431e4b2466a1
# ╠═a9b68546-cee8-4f6b-bd26-6fd6f1ca3bcd
# ╠═16a70b43-3f82-432e-ad53-6b2a6985669d
# ╠═eb0e11bd-5418-4923-a094-53222aee6ce7
# ╠═0c56f110-dcf5-4507-9459-00aad5f35b97
# ╠═294832c6-7ffa-43e7-8f46-5cbf4833d784
# ╠═82fdad6f-0c38-483e-9247-a51ccfe4a31b
# ╠═0f316ae1-e6be-446f-81d9-56c1a64fe988
# ╠═7a0bca86-f120-49ee-b2ed-7daee1508636
# ╠═f8402e1a-1b4b-4573-a591-ddda211cc0a4
# ╠═3c29ebab-2be4-4c5e-9fbe-a038f6dde3cc
# ╠═9f3825d7-081d-4fb3-863d-ff1f845b66ac
# ╠═2376caa4-3b68-4b52-aab4-d5fa786b525a
# ╠═e5c1feba-b371-4b9f-9585-406e484ccdc2
# ╠═7c993118-9c81-4fa3-a439-e8d898bf107d
# ╠═6da39c1e-3885-4bca-923f-2d47cf7b48d7
# ╠═2c90ed9c-3128-4535-8aaf-04a9660a3607
# ╠═86051cb1-fe25-4029-bb24-2610cd7fe646
# ╠═a5797d58-845d-4358-80c0-f4d0a93db650
# ╠═0aa8c55c-c931-43e9-bf7a-f9576efac54e
# ╠═965c8192-1e24-4c9b-bc71-df8c9c7a26ed
# ╠═4919adc9-8c8c-4324-b5de-abdd92e64172
# ╠═03526651-e04b-43db-874f-b4a80757bbf2
# ╠═8975b62c-9ec9-4a9c-a9ec-5a4d4ea51196
# ╠═9d82fd3d-4e90-43ed-99fc-35930c05b66b
# ╠═88343dd9-d306-4a53-b28d-60086a8d70a9
# ╠═83a12404-352a-4ab1-96da-160e666c4cec
# ╠═64ab2bd2-da2e-4aae-b271-ad4463e640e7
# ╠═7d539e78-5756-4dad-8491-18bbfde24915
# ╠═3f29923c-3d21-44c5-940c-13328277b4d3
# ╠═e4bbaf3c-43bc-423b-ae5c-c72c62a95926
# ╠═e369d39c-f283-40df-971b-b81cc9429e95
# ╠═aa4d9b4c-499c-45c8-8812-a2a3339bcbf1
# ╠═18f041fe-c89c-45ac-b25a-3e8b5ccc0e8e
# ╠═f78bdd04-c4f6-480a-960a-9219abd6fb13
# ╠═c4202fbe-7f45-4b7a-a627-b5749072386b
# ╠═978be772-d6e2-4009-bfe2-e9ce7ee12248
# ╠═a9924d52-9674-4cd9-a0c3-26c6115a819b
# ╠═72ec35ac-992e-4cfb-b7bd-74239ae860ee
# ╠═d64d7c14-b6a1-4219-aaef-2f182d3ff0a4
# ╠═823f2c63-a3dd-44ec-b0dd-e3a8e60f147c
# ╠═ee9884d4-0256-43de-b20e-0d4a4fa2c7f1
# ╠═dbda5eb8-3183-4ab9-8906-a7a3524f8fa5
# ╠═e68c51e4-e78a-4994-b19d-ec2048437054
# ╠═7b58081f-a2e0-477c-a313-f89453dd81a6
# ╠═8fadf7ab-9a2c-49ea-89e3-d64aea1baf81
# ╠═fce159d5-e6e7-4ca4-9f6b-5b5d9be6e085
# ╠═9beee121-1392-48ae-9ce4-cd881d0d4958
# ╠═8b4aff6a-8c6e-42dc-9da7-f1215a2c6d92
# ╠═991a69af-5163-4721-9dc5-53965dba6508
# ╠═954ce38b-f89f-411b-9621-6e6319c99123
# ╠═362d7ec1-a33a-480a-a3cd-353d2d882e96
# ╠═5b8c66ec-a038-42d8-be20-91281cb52ee4
# ╠═42550f71-0d65-4569-8332-41774a979e67
# ╠═52f6a0e2-a5f6-452e-ba45-14958b9ee3a4
# ╠═b7793747-d3ee-49eb-a88a-68f730358931
# ╠═519399b4-8f8c-4fcb-bfe3-4f7d7d67e761
# ╠═9093a05f-ecfd-497b-95b1-7530b5a49541
# ╠═7fa8480c-ac4e-4dec-be26-52b49a343b9f
# ╠═99eb99db-56f7-4884-bb1a-ca362fee1ba0
# ╠═f43ddffc-c0da-4140-9ddd-5668708129a4
# ╠═5fc2c3b0-3cdf-4bdd-9f83-408a68a78f97
# ╠═8efc9eb5-603e-4799-92f0-46a30de30dc3
# ╠═c0816fb9-00c7-41c4-b2e4-da4dc50fdbd1
# ╠═88aab6f6-07a8-42c4-b39f-33ed7bb0a75b
# ╠═ee8468ab-2987-497c-83ac-720902561490
# ╠═c0ba6116-ae04-4ab5-a5c8-08d623444a3e
# ╠═2d18834e-70f2-47f0-9238-026d1d6fe462
