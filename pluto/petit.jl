### A Pluto.jl notebook ###
# v0.19.14

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
	using StatsBase
	using Distributions
	using Unitful 
	#using UnitfulEquivalences 
	#using PhysicalConstants
	using Peaks
	#using FFTW
	#using DSP
	#using Clustering
	using HDF5
	using Roots
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
	xpath = "/Users/jjgomezcadenas/LaserLab/PETALO/petitr1/"
	xfiles = Glob.glob(dfiles, xpath)
end

# ╔═╡ 42550f71-0d65-4569-8332-41774a979e67
md"""
## Control plots
"""

# ╔═╡ 59a1dffb-6902-4cb0-8650-a23498794f49
md"""
## Single SiPM analysis
"""

# ╔═╡ bd496019-51e2-46ee-96e7-048bf53dc547
md""" Plane 2: select rmin histogram: $(@bind rmin NumberField(50.0:650.0, default=300.0))"""

# ╔═╡ 99eb99db-56f7-4884-bb1a-ca362fee1ba0
md""" Plane 2: select rmax histogram: $(@bind rmax NumberField(400.0:600.0, default=500.0))"""

# ╔═╡ a88457ad-0f1c-4ca8-babf-5eaf3f040de7
rbins = Int(rmax - rmin)

# ╔═╡ f43ddffc-c0da-4140-9ddd-5668708129a4
md""" Plane 2: select width of peak for fit: $(@bind wpeak NumberField(5.0:20.0, default=10.0))"""

# ╔═╡ 1b9a7002-8ad6-42aa-a691-620923daedf0
md"""
### ML unbinned fit 
"""

# ╔═╡ 7c4bf607-3d3c-4f70-9ffb-fe601ac92487
md"""
### Binned fit 
"""

# ╔═╡ 08c2ea5c-137a-42db-9043-135cc4b3231e
md"""
## Multi SiPM analysis
"""

# ╔═╡ 44eb878c-7c78-42e9-87b1-32d919ae2c3f
bad_sipms=[116,119,120, 123, 129, 130, 133, 134, 139, 140, 149, 150, 159, 160, 169, 170, 171, 179, 180, 181]

# ╔═╡ 61fb8aab-598c-4dd7-9b45-4ade34436d05
all_sipms = collect(111:188)

# ╔═╡ e451863b-1c4f-49de-9c6c-7d2f41388d37
nsipms = setdiff(all_sipms, bad_sipms) #sipms in all not in bad

# ╔═╡ 7ac98d50-0edd-40c9-94bf-eb4fb8cdf706
md"""
### Charge
"""

# ╔═╡ 60101f34-d408-4b6e-9813-33edafcdabcf
md"""
### Unbinned resolution 
"""

# ╔═╡ 5fc2c3b0-3cdf-4bdd-9f83-408a68a78f97


# ╔═╡ 5b05ebab-e416-4550-8dcb-a8a3fcf41f2f
#begin
	#hrfw, prfw = hist1d(rml, "R (FWHM)", 20, 0.0, 0.05)
	#plot(prfw)
#end

# ╔═╡ 120a9a67-45dc-410f-85cf-2e66f679b276
#for sp in nsipms[1]
#	println("sp = ", sp)
#	rmx =rfmle(sp)
#	@printf "R = %7.4f" rmx
#	println("")
#end

# ╔═╡ 1004c6bb-07bf-4334-a0e8-feeb0d94e295
md"""
### Binned
"""

# ╔═╡ 8efc9eb5-603e-4799-92f0-46a30de30dc3


# ╔═╡ a49f93d5-5106-48ac-8107-4f579374b85d
md"""
## Resolution intrinsic
"""

# ╔═╡ 26db5f86-93ac-440a-b6d5-cce090f4866c
RFWHM(σm, σp) = sqrt(σm^2 - σp^2)

# ╔═╡ d3a19c59-602f-4e26-9a19-a0b8ebb2584b
RFWHM(0.03, 0.025) 

# ╔═╡ 2f125bdc-c551-49d3-959f-cc539be6d69c
sqrt(0.03^2 - 0.0255^2)

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
	dset =fid["analysis"]["table"]
	dda = read(dset)
	ddic = Dict()
    for (i, key) in enumerate(keys(dda[1]))
        ddic[key] = [ddi[i] for ddi in dda]
    end
	close(fid)
    DataFrame(ddic)
end

# ╔═╡ 71f427b9-4a44-492c-a826-1bb403fcf7aa
df2 = get_dfs("/Users/jjgomezcadenas/LaserLab/PETALO/petitr1/channel_analysis_he_12334_10_1.h5")

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

# ╔═╡ 082a3d2a-c915-4c98-9221-d4610ae112bc
df = concdf(xfiles)

# ╔═╡ 52f6a0e2-a5f6-452e-ba45-14958b9ee3a4
begin
	hmc1, pmc1 = hist1d(df.max_charge0, "max charge plane 1", 50, 100., 500.0)
	hmsc1, pmsc1 = hist1d(df.sum_charge0, "max charge plane 1", 50, 100., 500.0)
	plot(pmc1, pmsc1)
end

# ╔═╡ 7690127a-b6ab-4661-aea3-6be93c0f8760
begin
	hmc2, pmc2 = hist1d(df.max_charge2, "max charge plane 2", 50, 100., 500.0)
	hmsc2, pmsc2 = hist1d(df.sum_charge2, "max charge plane 2", 50, 100., 500.0)
	plot(pmc2, pmsc2)
end

# ╔═╡ 519399b4-8f8c-4fcb-bfe3-4f7d7d67e761
begin
	hmc12, pmc12 = hist2d(df.max_charge0,df.max_charge2, 25,
                "Q1", "Q2", 100., 500.0, 100., 500.0; title="Q1-Q2")
	hmsc12, pmsc12 = hist2d(df.sum_charge0,df.sum_charge2, 25,
                "Q1", "Q2", 100., 500.0, 100., 500.0; title="Q1-Q2")
	plot(pmc12, pmsc12)
end

# ╔═╡ 58670f80-1bf2-4af1-a588-42d52f0394fb
begin
	snsp1 = sort(unique(df.sns_id0))
	snsp2 = sort(unique(df.sns_id2))
	gnspm2 = groupby(df, :sns_id2)
	md""" 
	- Plane1 SiPms betwwen $(snsp1[1]) and $(snsp1[end])
	- Plane2 SiPms betwwen $(snsp2[1]) and $(snsp2[end])
	"""
end

# ╔═╡ 4271b7ae-62ad-42e9-a298-b7d9cb07bb86
length(snsp2)

# ╔═╡ f2476981-326f-4f7e-b796-90a9f357ee4a
md""" Plane 2: select sipm index: $(@bind nspm NumberField(snsp2[1]:snsp2[end], default=snsp2[1]))"""


# ╔═╡ de099068-fbe4-428a-89d7-6b6119c52fb9
begin
	spmdf = gnspm2[(nspm,)]
	hspmc2, pspmc2 = hist1d(Float64.(spmdf.max_charge2), "max charge plane 2", 250, 100., 600.0)
	#hspsc2, pspsc2 = hist1d(Float64.(spmdf.sum_charge2), "sum charge plane 2", 50, 100., 600.0)
	plot(pspmc2)
end

# ╔═╡ 768d7845-a601-46e4-a86c-09d98d35f715
begin
	h2spmc2, p2spmc2 = hist1d(Float64.(spmdf.max_charge2), "max charge plane 2", rbins, rmin, rmax)
	plot(p2spmc2)
end

# ╔═╡ 3dbe8495-cee9-4d1a-8a9a-4003db5e9149
"""
Compute charge 
"""
function qsipms(nspmx)
	spmdf = gnspm2[(nspmx,)]
	Float64.(spmdf.max_charge2)
end


# ╔═╡ a578676e-3daa-46f0-90dd-1699b62d4605
qs = map(x->qsipms(x), nsipms)

# ╔═╡ 82cd3877-56bb-47ae-b0fd-0f2062cd4386
qsm = map(x->maximum(x), qs)

# ╔═╡ 85760432-a6ed-4823-8e00-4835571ed340
begin
	hqs, pqs = hist1d(qsm, "Q (LSB)", 25, 400.0, 600.0; datap=false)
	plot(pqs)
end

# ╔═╡ 10184c1e-601c-4763-98b7-73a42f64e7e5
length(qs)

# ╔═╡ 0ad899a0-5848-4f27-ada1-fc9aa5e7339b
"""
Finds the positin of the photopeak as the peak with larger prominence in histogram
"""
function find_photopeak(h)
	pks, vals = findmaxima(h.weights)
	peaks, proms = peakproms(pks, h.weights)
	mxprom = findmax(proms)
	ppeak = pks[mxprom[2]]
	centers(h2spmc2)[ppeak]
end

# ╔═╡ 7bea5385-9bc5-4eb5-86ef-d3030b9e2159
xpeak = find_photopeak(h2spmc2)

# ╔═╡ 23c2e8d7-18de-4494-94e6-4acc2aff3127
begin
	xmin = xpeak - wpeak 
	xmax = xpeak + wpeak +2 # asymmetric on the right 
end

# ╔═╡ db275e86-b644-4982-8666-6eda78ed5072
md""" Plane 2: select xmin fit: $(@bind sxmin NumberField(100.0:600.0, default=xmin))"""

# ╔═╡ b8cb9917-bb46-43a8-95a1-a05e6c2d0d7f
md""" Plane 2: select max fit: $(@bind sxmax NumberField(100.0:600.0, default=xmax))"""

# ╔═╡ ca6b160a-bbea-4ae0-8ab0-de31c387d439
ftbins = Int(sxmax - sxmin)

# ╔═╡ 65bb374d-aa0a-4c97-aab7-1be14120d137
md"""
- maximum of distribution = $(xpeak)
- semidiwth for fit = $(wpeak)

fit limits: 
- xmin = $(xmin)
- xmax = $(xmax)
- number of bins = $(ftbins)
"""

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

# ╔═╡ 54874b4d-2faf-47d5-b68f-17fb532de95c
"""
Performs a Maximum Likelihood fit to a Normal distribution

Returns the fit object and the loglikelihood

"""
function norm_mle_fit(x, xmin, xmax)
	fdset = filter(x->x>xmin && x<xmax, x)
	fmle = fit_mle(Normal, fdset)
	fmle, -loglikelihood(fmle, fdset)
end

# ╔═╡ e16ea9bc-7aa6-47a6-9a86-8d0f656dc002
ftmle = norm_mle_fit(Float64.(spmdf.max_charge2), sxmin, sxmax)

# ╔═╡ 41ec942b-4c04-4b43-8475-d4a14fae02c2
"""
Vary the paramters of the fit around the best value. Return a named
tuple (lft object) with information to compute parabolic fit
"""
function norm_mle_var(x::Vector{Float64}, xmin::Float64, xmax::Float64; var="σ",
                      sr=2.0, ss=100.0)
	
	fdset = filter(x->x>xmin && x<xmax, x)
	fmle = fit_mle(Normal, fdset)
	mll = -loglikelihood(fmle, fdset)
	
	if var=="μ"
		s0 = fmle.μ
	else
		s0 = fmle.σ
	end
	smin = s0 - sr
	smax = s0 + sr
	sx = (smax-smin) / ss
	srr = collect(smin:sx:smax)
	if var=="σ"
		dists = map(x->Normal(fmle.μ, x), srr)
	else
		dists = map(x->Normal(x, fmle.σ), srr)
	end
	
	mlls = map(x->-loglikelihood(x, fdset), dists)
	(fbest = fmle, llbest = mll, vars = srr, logs = mlls)
end

# ╔═╡ 328768e0-c381-4060-b547-282ad0ae43dc
lft = norm_mle_var(Float64.(spmdf.max_charge2), sxmin, sxmax, var="σ", sr=0.15)

# ╔═╡ b72281dd-e4b1-4c0f-ab1d-8d9df19aa837
begin
	plot(lft.vars, lft.logs, lw=2, label="")
	scatter!([lft.fbest.σ], [lft.llbest], label="best fit", legend=true)
	hline!([lft.llbest + 0.5], label="")
	xlabel!("σ")
	ylabel!("-log(L)")
end

# ╔═╡ fe759a1a-4ba2-4b55-ad56-d9859e6f8d3a
lft2 = norm_mle_var(Float64.(spmdf.max_charge2), sxmin, sxmax, var="μ", sr=0.25)

# ╔═╡ 1ff412ee-1613-4d2e-9d22-7205ad7ffdd4
begin
	plot(lft2.vars, lft2.logs, lw=2, label="")
	scatter!([lft2.fbest.μ], [lft.llbest], label="best fit", legend=true)
	hline!([lft.llbest + 0.5], label="")
	xlabel!("μ")
	ylabel!("-log(L)")
end

# ╔═╡ 0fadbb90-5815-44e7-b3bf-8f624a087cbd
"""
Take an lft obect and return best value and error for sigma 
"""
function fsigma(lft)
	fll = v2f(lft.vars, lft.logs)
	fzero(x) = fll(x) - lft.llbest - 0.5
	sigma = find_zero(fzero , (lft.fbest.σ,lft.vars[end]), Bisection()) - find_zero(fzero , (lft.vars[1],lft.fbest.σ), Bisection())
	(var=lft.fbest.σ, delta=sigma)
end

# ╔═╡ 25f704f9-f280-4f97-aeba-2ad3856c5c4e
lftsigma = fsigma(lft)

# ╔═╡ 15cfdbbb-436e-4e19-b036-a7527521e094
lftsigma.delta/lftsigma.var

# ╔═╡ 5f18aef7-c557-458a-a14f-15f9986924b4
"""
Take an lft obect and return best value and error for mu 
"""
function fmu(lft)
	fll = v2f(lft.vars, lft.logs)
	fzero(x) = fll(x) - lft.llbest - 0.5
	sigma = find_zero(fzero , (lft.fbest.μ,lft.vars[end]), Bisection()) - find_zero(fzero , (lft.vars[1],lft.fbest.μ), Bisection())
	(var=lft.fbest.μ, delta=sigma)
end

# ╔═╡ 6f34b6ab-ba5a-4e5b-aee8-f53b78c4e260
lftmu = fmu(lft2)

# ╔═╡ 42e41d5a-bbf4-4423-ab3d-0e8afc326990
fg, pg = fitg1(Float64.(spmdf.max_charge2),"E (LSB)", rbins, rmin, rmax;
	           xgmin=sxmin, xgmax=sxmax, fbins=ftbins, norm=true,
               fm=lftmu.var, flex_mean=true)

# ╔═╡ 026d5f3e-0e02-4692-8cbd-3af88a61f603
plot(pg, legend=false)

# ╔═╡ af825b53-e0a3-4296-ab5b-a467800c18ee
"""
Takes a fit object and returns relative resolution FWHM (2.355 * sigma/mu)
"""
function rfromft(lftmu, lftsigma)
	R = 2.355 * lftsigma.var/lftmu.var
	δR = R * lftsigma.delta/lftsigma.var  # mu is negligible
	(var=R, delta=δR)
end

# ╔═╡ 2d18834e-70f2-47f0-9238-026d1d6fe462
function rfromft(fg::ATools.FGauss)
	R = 2.355 * fg.std[1]/fg.mu[1]
	δR = R * fg.δstd[1]/fg.std[1]  # mu is negligible
	(var=R, delta=δR)
end

# ╔═╡ c70f9dde-2cb3-4057-b928-b02d1d3b46cf
rfft = rfromft(lftmu, lftsigma)

# ╔═╡ b408ba82-5b66-4200-8300-dbf50d67caf1
md"""
Fit result:
- μ = $(round(lftmu.var, sigdigits=3)) +- $(round(lftmu.delta, sigdigits=2))
- σ = $(round(lftsigma.var, sigdigits=3)) +- $(round(lftsigma.delta, sigdigits=2))
- R = $(round(rfft.var, sigdigits=3)) +- $(round(rfft.delta,sigdigits=2))
"""

# ╔═╡ 62806e24-1db6-4f57-8a46-c0a7875e431c
R, dR = rfromft(fg)

# ╔═╡ 540fa65f-0d85-457c-a6f2-ad4bc1799408
"""
    gausg(μ::Real, σ::Real, C::Real)::Function
Return a Gaussian function with parameters μ, σ and normalisation C.
"""
function gausg(μ::Real, σ::Real, C::Real)::Function
	function gausx(x::Real)
		return C * pdf(Normal(μ, σ,), x)
	end
	return gausx
end

# ╔═╡ 38e444f9-1bc1-43ac-972b-75e072ca8cac
"""
Histogram vector x, fit a Normal distribution and plot
"""
function histfit(x::Vector{Float64}, nbins::Int64, xmin::Float64, xmax::Float64,
                  ftbins::Int64, ftmin::Float64, ftmax::Float64, norm=false)

	fdset = filter(x->x>ftmin && x<ftmax, x)
	fmle = fit_mle(Normal, fdset)
	
	@. gauss1(x, p)   = p[1] * pdf(Normal(fmle.μ, fmle.σ), x)

	#@debug "fit min. likelihood"  fmle
	
	hf   = hist1d(x, ftbins, ftmin, ftmax, norm)
	h, p = hist1d(x, "X", nbins, xmin, xmax, norm=norm, legend=true)
	c = centers(hf)
	w = eltype(c).(hf.weights)
	

	# fit parameters lb, ub, po are lower, upper bounds and pars
    lb = [         0.0]
    ub = [100 * sum(w)]
    p0 = [      sum(w)]
	fq = curve_fit(gauss1, c, w, p0, lower=lb, upper=ub)
	CC = coef(fq)[1]
	
	@info " CC, μ, σ, R" CC, fmle.μ, fmle.σ, 100.0*2.355 * fmle.σ/ fmle.μ
	gx = gausg(fmle.μ, fmle.σ, CC)
	lbl = @sprintf "R =%4.5f  FWHM" 2.355 * fmle.σ/ fmle.μ
    p   = plot!(p, c, gx.(c), lw=2, label=lbl, legend=true, fmt = :png)
	xlabel!("E (LSB)")
	(μ=fmle.μ, σ=fmle.σ, C = CC, gx= gx, plt=p)

end

# ╔═╡ cac298c3-c090-4884-82d2-20cb2bc53e64
"""
	Compute the chi2 of a binned fit to x, with nbins between xmin and xmax
"""
function fitchi2(x::Vector{Float64}, nbins::Int64, xmin::Float64, xmax::Float64, gx; 
                 imin::Int64=1)
	h = hist1d(x, nbins, xmin, xmax)
	c = centers(h)
	O = h.weights
	E = gx.(c)
	num = [(O[i] - E[i])^2 for i in imin:length(c)]
	#num = (O .- E) .* (O .- E)
	denom = E[imin:end]
	ndf = length(c) - imin
	chi2 = sum(num ./ denom) /ndf
	(C=c, E = E, O = O, chi2=chi2)
end

# ╔═╡ d63eee83-0284-4ec6-95ba-f2f93d5f9ffc
"""
	Compute the chi2 of a binned fit to x, with nbins between xmin and xmax
"""
function fitchi2(fg::ATools.FGauss; imin::Int64=1)
	h = fg.h
	c = centers(h)
	O = h.weights
	E = fg.Y
	num = [(O[i] - E[i])^2 for i in imin:length(c)]
	denom = E[imin:end]
	ndf = length(c) - imin
	chi2 = sum(num ./ denom) /ndf
	(C=c, E = E, O = O, chi2=chi2)
end

# ╔═╡ 2de40ccb-d795-449a-9939-c0f23cbbe1bc
fch2= fitchi2(fg, imin=4)

# ╔═╡ 8d002405-6bb7-450d-8368-1026df3190b9
begin
	scatter(fch2.C, fch2.E,  yerr=sqrt.(fch2.E), label="fit")
	scatter!(fch2.C, fch2.O,  yerr=sqrt.(fch2.O), label="data")
	plot!(fch2.C, fg.g[1].(fch2.C), label="gauss")
end

# ╔═╡ b17ef99c-221a-47f7-9080-9e04c722281e
md"""
Binned Fit results:
- μ = $(round(fg.mu[1], sigdigits=2))
- σ = $(round(fg.std[1], sigdigits=2))
- chi2 = $(round(fch2.chi2, sigdigits=2))
- R = $(round(R, sigdigits=3)) +- $(round(dR, sigdigits=3)) FWHM
"""

# ╔═╡ 9294ca23-4a6a-465e-a840-5edc04dbe2eb
"""
Compute R (FWHM) using ML
"""
function rfmle(nspmx, nbin=200, bmin=300.0, bmax=500.0; ftm="unbin", wpeak=10)
	#println("nsipm =", nspmx)
	spmdf = gnspm2[(nspmx,)]
	h2spmc2, _ = hist1d(Float64.(spmdf.max_charge2), "max charge plane 2", nbin, bmin, bmax)
	
	xpeak = find_photopeak(h2spmc2)
	xmin = ceil(xpeak - wpeak)
	xmax = ceil(xpeak + wpeak +2)  # asymmetric on the right 
	ftbins = Int(xmax - xmin)
	
	#println("xmin =", xmin, " xmax = ", xmax, " xpeak = ", xpeak, " bins =", ftbins)

	if ftm == "unbin"
		ftmle = norm_mle_fit(Float64.(spmdf.max_charge2), xmin, xmax)
		lft   = norm_mle_var(Float64.(spmdf.max_charge2), xmin, xmax, var="σ", sr=0.15)
		lftsigma = fsigma(lft)
		lft2 = norm_mle_var(Float64.(spmdf.max_charge2), xmin, xmax, var="μ", sr=0.25)
		lftmu = fmu(lft2)
		rfft = rfromft(lftmu, lftsigma)
		return rfft.var
	else

		fg, pg = fitg1(Float64.(spmdf.max_charge2),"E (LSB)", nbin, bmin, bmax;
		               xgmin=xmin, xgmax=xmax, fbins=ftbins, flex_mean=true)
		R, dR = rfromft(fg)
		fch2= fitchi2(fg, imin=4)
		return R
	end
end

# ╔═╡ 95afa3bc-7cbc-4f93-bc31-5829a96037f2
rml = map(x->rfmle(x, ftm="unbin"), nsipms)

# ╔═╡ 5fbae887-ff4f-4aa0-823e-f1284cc5d0b2
begin
	srq = scatter(qsm, rml, legend=false)
	xlabel!("Q (LSB)")
	ylabel!("σ (FWHM)")
end

# ╔═╡ 41005dd9-3303-4994-b78f-13580fa90a86
rmlb = map(x->rfmle(x, ftm="bin"), nsipms)

# ╔═╡ b5a5a985-1111-43f9-8cce-5cd68531c735
begin
	srml = scatter(nsipms, rml, label="R unbin")
	srmlb = scatter(nsipms, rmlb, label="R bin")
	plot(size=(600,600),srml, srmlb)
end

# ╔═╡ 12ce6070-a101-410d-a018-0df3633bef0d
begin
	hrfw, prfw = hist1d(rml, "R (FWHM)", 40, 0.0, 0.08; datap=false)
	hrfwb, prfwb = hist1d(rmlb, "R (FWHM)", 40, 0.0, 0.08; datap=false)
	plot( prfw, prfwb)
end

# ╔═╡ d63c1850-a05c-414c-b132-460381bb8767
plot( prfw, legend=false)

# ╔═╡ 8c87548d-6bbb-4130-adf9-3ef4c2ec91ae
md"""
- Unbinned:  R = $(round(mean(rml), sigdigits=3)) +- $(round(std(rml), sigdigits=3))
- Binned: R = $(round(mean(rmlb), sigdigits=3)) +- $(round(std(rmlb), sigdigits=3))
"""

# ╔═╡ 3d29bc36-530c-4999-b9fe-58bf46a6c712
begin
	rml2 = map(x->rfmle(x, ftm="unbin", wpeak=11), nsipms)
	rmlb2 = map(x->rfmle(x, ftm="bin", wpeak=11), nsipms)
	hrfw2, prfw2 = hist1d(rml2, "R (FWHM)", 40, 0.0, 0.08)
	hrfwb2, prfwb2 = hist1d(rmlb2, "R (FWHM)", 40, 0.0, 0.08)
	plot(prfw2, prfwb2)
end

# ╔═╡ 632a6398-8213-45c1-8c29-fa419e1c1739
md"""
- Unbinned:  R = $(round(mean(rml2), sigdigits=3)) +- $(round(std(rml2), sigdigits=3))
- Binned: R = $(round(mean(rmlb2), sigdigits=3)) +- $(round(std(rmlb2), sigdigits=3))
"""

# ╔═╡ adfc0a99-656c-4b3c-aec3-e0a2ceb4cbec
begin
	rml3 = map(x->rfmle(x, ftm="unbin", wpeak=9), nsipms)
	rmlb3 = map(x->rfmle(x, ftm="bin", wpeak=9), nsipms)
	hrfw3, prfw3 = hist1d(rml3, "R (FWHM)", 40, 0.0, 0.08, datap=false)
	hrfwb3, prfwb3 = hist1d(rmlb3, "R (FWHM)", 40, 0.0, 0.08)
	plot(prfw3, prfwb3)
end

# ╔═╡ 36d29dd1-7df7-457f-97a4-70db9e707c3b
md"""
- Unbinned:  R = $(round(mean(rml3), sigdigits=3)) +- $(round(std(rml3), sigdigits=3))
- Binned: R = $(round(mean(rmlb3), sigdigits=3)) +- $(round(std(rmlb3), sigdigits=3))
"""

# ╔═╡ da19c1d9-0682-47de-b69c-50a9938fe897
plot(prfw3, legend=false)

# ╔═╡ e1ce8dd8-e66e-444e-bf04-a8d1378b6d71
fg

# ╔═╡ Cell order:
# ╠═805955a4-f018-4e68-a68a-1e6b64a99378
# ╠═53dee408-58f8-11ed-2b44-db5e37f3a91e
# ╠═a69771ca-f756-4f84-8bb4-1c044fa4db1f
# ╠═f75fa3dd-d3dc-4f5d-877b-cc9d8f57d95c
# ╠═921fb297-2c3d-4f45-847c-b2f009052e9f
# ╠═918dc9ec-f027-4d68-8290-96f7f56a416a
# ╠═3274825d-61a0-4856-9340-97420d237c27
# ╠═082a3d2a-c915-4c98-9221-d4610ae112bc
# ╠═71f427b9-4a44-492c-a826-1bb403fcf7aa
# ╠═42550f71-0d65-4569-8332-41774a979e67
# ╠═52f6a0e2-a5f6-452e-ba45-14958b9ee3a4
# ╠═7690127a-b6ab-4661-aea3-6be93c0f8760
# ╠═519399b4-8f8c-4fcb-bfe3-4f7d7d67e761
# ╠═59a1dffb-6902-4cb0-8650-a23498794f49
# ╠═4271b7ae-62ad-42e9-a298-b7d9cb07bb86
# ╠═58670f80-1bf2-4af1-a588-42d52f0394fb
# ╠═f2476981-326f-4f7e-b796-90a9f357ee4a
# ╠═de099068-fbe4-428a-89d7-6b6119c52fb9
# ╠═bd496019-51e2-46ee-96e7-048bf53dc547
# ╠═99eb99db-56f7-4884-bb1a-ca362fee1ba0
# ╠═a88457ad-0f1c-4ca8-babf-5eaf3f040de7
# ╠═768d7845-a601-46e4-a86c-09d98d35f715
# ╠═7bea5385-9bc5-4eb5-86ef-d3030b9e2159
# ╠═f43ddffc-c0da-4140-9ddd-5668708129a4
# ╠═23c2e8d7-18de-4494-94e6-4acc2aff3127
# ╠═65bb374d-aa0a-4c97-aab7-1be14120d137
# ╠═db275e86-b644-4982-8666-6eda78ed5072
# ╠═b8cb9917-bb46-43a8-95a1-a05e6c2d0d7f
# ╠═ca6b160a-bbea-4ae0-8ab0-de31c387d439
# ╠═1b9a7002-8ad6-42aa-a691-620923daedf0
# ╠═e16ea9bc-7aa6-47a6-9a86-8d0f656dc002
# ╠═328768e0-c381-4060-b547-282ad0ae43dc
# ╠═b72281dd-e4b1-4c0f-ab1d-8d9df19aa837
# ╠═25f704f9-f280-4f97-aeba-2ad3856c5c4e
# ╠═15cfdbbb-436e-4e19-b036-a7527521e094
# ╠═fe759a1a-4ba2-4b55-ad56-d9859e6f8d3a
# ╠═1ff412ee-1613-4d2e-9d22-7205ad7ffdd4
# ╠═6f34b6ab-ba5a-4e5b-aee8-f53b78c4e260
# ╠═c70f9dde-2cb3-4057-b928-b02d1d3b46cf
# ╠═b408ba82-5b66-4200-8300-dbf50d67caf1
# ╠═7c4bf607-3d3c-4f70-9ffb-fe601ac92487
# ╠═42e41d5a-bbf4-4423-ab3d-0e8afc326990
# ╠═026d5f3e-0e02-4692-8cbd-3af88a61f603
# ╠═62806e24-1db6-4f57-8a46-c0a7875e431c
# ╠═2de40ccb-d795-449a-9939-c0f23cbbe1bc
# ╠═8d002405-6bb7-450d-8368-1026df3190b9
# ╠═b17ef99c-221a-47f7-9080-9e04c722281e
# ╠═08c2ea5c-137a-42db-9043-135cc4b3231e
# ╠═44eb878c-7c78-42e9-87b1-32d919ae2c3f
# ╠═61fb8aab-598c-4dd7-9b45-4ade34436d05
# ╠═e451863b-1c4f-49de-9c6c-7d2f41388d37
# ╠═7ac98d50-0edd-40c9-94bf-eb4fb8cdf706
# ╠═a578676e-3daa-46f0-90dd-1699b62d4605
# ╠═82cd3877-56bb-47ae-b0fd-0f2062cd4386
# ╠═85760432-a6ed-4823-8e00-4835571ed340
# ╠═60101f34-d408-4b6e-9813-33edafcdabcf
# ╠═95afa3bc-7cbc-4f93-bc31-5829a96037f2
# ╠═10184c1e-601c-4763-98b7-73a42f64e7e5
# ╠═5fc2c3b0-3cdf-4bdd-9f83-408a68a78f97
# ╠═5fbae887-ff4f-4aa0-823e-f1284cc5d0b2
# ╠═5b05ebab-e416-4550-8dcb-a8a3fcf41f2f
# ╠═120a9a67-45dc-410f-85cf-2e66f679b276
# ╠═1004c6bb-07bf-4334-a0e8-feeb0d94e295
# ╠═41005dd9-3303-4994-b78f-13580fa90a86
# ╠═8efc9eb5-603e-4799-92f0-46a30de30dc3
# ╠═b5a5a985-1111-43f9-8cce-5cd68531c735
# ╠═12ce6070-a101-410d-a018-0df3633bef0d
# ╠═d63c1850-a05c-414c-b132-460381bb8767
# ╠═8c87548d-6bbb-4130-adf9-3ef4c2ec91ae
# ╠═3d29bc36-530c-4999-b9fe-58bf46a6c712
# ╠═632a6398-8213-45c1-8c29-fa419e1c1739
# ╠═adfc0a99-656c-4b3c-aec3-e0a2ceb4cbec
# ╠═36d29dd1-7df7-457f-97a4-70db9e707c3b
# ╠═da19c1d9-0682-47de-b69c-50a9938fe897
# ╠═a49f93d5-5106-48ac-8107-4f579374b85d
# ╠═26db5f86-93ac-440a-b6d5-cce090f4866c
# ╠═d3a19c59-602f-4e26-9a19-a0b8ebb2584b
# ╠═2f125bdc-c551-49d3-959f-cc539be6d69c
# ╠═c0816fb9-00c7-41c4-b2e4-da4dc50fdbd1
# ╠═3dbe8495-cee9-4d1a-8a9a-4003db5e9149
# ╠═9294ca23-4a6a-465e-a840-5edc04dbe2eb
# ╠═88aab6f6-07a8-42c4-b39f-33ed7bb0a75b
# ╠═ee8468ab-2987-497c-83ac-720902561490
# ╠═0ad899a0-5848-4f27-ada1-fc9aa5e7339b
# ╠═c0ba6116-ae04-4ab5-a5c8-08d623444a3e
# ╠═54874b4d-2faf-47d5-b68f-17fb532de95c
# ╠═41ec942b-4c04-4b43-8475-d4a14fae02c2
# ╠═0fadbb90-5815-44e7-b3bf-8f624a087cbd
# ╠═5f18aef7-c557-458a-a14f-15f9986924b4
# ╠═af825b53-e0a3-4296-ab5b-a467800c18ee
# ╠═2d18834e-70f2-47f0-9238-026d1d6fe462
# ╠═540fa65f-0d85-457c-a6f2-ad4bc1799408
# ╠═38e444f9-1bc1-43ac-972b-75e072ca8cac
# ╠═cac298c3-c090-4884-82d2-20cb2bc53e64
# ╠═d63eee83-0284-4ec6-95ba-f2f93d5f9ffc
# ╠═e1ce8dd8-e66e-444e-bf04-a8d1378b6d71
