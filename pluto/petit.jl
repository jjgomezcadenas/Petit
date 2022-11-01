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
	#using Peaks
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

# ╔═╡ 3274825d-61a0-4856-9340-97420d237c27
begin
	dfiles ="*.h5"
	xpath = "/Users/jjgomezcadenas/LaserLab/PETALO/petitr1/"
	xfiles = Glob.glob(dfiles, xpath)
end

# ╔═╡ a2baefc2-236c-4b98-a9d1-1506feacaee0


# ╔═╡ bd496019-51e2-46ee-96e7-048bf53dc547
md""" Plane 2: select rmin histogram: $(@bind rmin NumberField(300.0:375.0, default=330.0))"""

# ╔═╡ 99eb99db-56f7-4884-bb1a-ca362fee1ba0
md""" Plane 2: select rmax histogram: $(@bind rmax NumberField(400.0:450.0, default=430.0))"""

# ╔═╡ a88457ad-0f1c-4ca8-babf-5eaf3f040de7
rbins = Int(rmax - rmin)

# ╔═╡ db275e86-b644-4982-8666-6eda78ed5072
md""" Plane 2: select xmin fit: $(@bind sxmin NumberField(350.0:450.0, default=375.0))"""

# ╔═╡ b8cb9917-bb46-43a8-95a1-a05e6c2d0d7f
md""" Plane 2: select max fit: $(@bind sxmax NumberField(350.0:450.0, default=400.0))"""

# ╔═╡ ca6b160a-bbea-4ae0-8ab0-de31c387d439
ftbins = Int(sxmax - sxmin)

# ╔═╡ a71cbe76-0a65-4cf9-be46-0cc695b55204
md"""
test
	-a
 	-b
"""

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

# ╔═╡ f2476981-326f-4f7e-b796-90a9f357ee4a
md""" Plane 2: select sipm index: $(@bind nspm NumberField(snsp2[1]:snsp2[end], default=snsp2[1]))"""


# ╔═╡ de099068-fbe4-428a-89d7-6b6119c52fb9
begin
	spmdf = gnspm2[(nspm,)]
	hspmc2, pspmc2 = hist1d(Float64.(spmdf.max_charge2), "max charge plane 2", 50, 100., 600.0)
	hspsc2, pspsc2 = hist1d(Float64.(spmdf.sum_charge2), "sum charge plane 2", 50, 100., 600.0)
	plot(pspmc2, pspsc2)
end

# ╔═╡ 768d7845-a601-46e4-a86c-09d98d35f715
begin
	h2spmc2, p2spmc2 = hist1d(Float64.(spmdf.max_charge2), "max charge plane 2", 100, rmin, rmax)
	plot(p2spmc2)
end

# ╔═╡ 67acec5f-a69e-4af7-98e8-8c8fcdec77b1
argmax(h2spmc2.weights)

# ╔═╡ 0020aaac-a114-4cfb-bdba-955b4f0fdabd
centers(h2spmc2)[argmax(h2spmc2.weights)]

# ╔═╡ 65bb374d-aa0a-4c97-aab7-1be14120d137
md"""
fit limits: 
- xmin = $(centers(h2spmc2)[argmax(h2spmc2.weights)] -12)
- xmax = $(centers(h2spmc2)[argmax(h2spmc2.weights)] +12)
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

# ╔═╡ af825b53-e0a3-4296-ab5b-a467800c18ee
"""
Takes a fit object and returns relative resolution FWHM (2.355 * sigma/mu)
"""
function rfromft(lftmu, lftsigma)
	R = 2.355 * lftsigma.var/lftmu.var
	δR = R * lftsigma.delta/lftsigma.var  # mu is negligible
	(var=R, delta=δR)
end

# ╔═╡ c70f9dde-2cb3-4057-b928-b02d1d3b46cf
rfft = rfromft(lftmu, lftsigma)

# ╔═╡ b408ba82-5b66-4200-8300-dbf50d67caf1
md"""
Fit result:
- μ = $(round(lftmu.var, sigdigits=3)) +- $(round(lftmu.delta, sigdigits=2))
- σ = $(round(lftsigma.var, sigdigits=2)) +- $(round(lftsigma.delta, sigdigits=2))
- R = $(round(rfft.var, sigdigits=2)) +- $(round(rfft.delta,sigdigits=2))
"""

# ╔═╡ 2d18834e-70f2-47f0-9238-026d1d6fe462
#@. gauss1(x, p)   = p[1] * pdf(Normal(p[2], p[3]), x)

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
function plotfit(x::Vector{Float64}, nbins::Int64, xmin::Float64, xmax::Float64,
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
	p

end

# ╔═╡ c7879c6a-c905-4cd0-9661-2cb924a09702
plotfit(Float64.(spmdf.max_charge2), rbins, rmin, rmax, 
               ftbins, sxmin, sxmax)

# ╔═╡ cac298c3-c090-4884-82d2-20cb2bc53e64
"""
	Compute the chi2 of a binned fit to x, with nbins between xmin and xmax
"""
function fitchi2(x, fg, nbins, xmin, xmax)
	h = hist1d(x, nbins, xmin+0.5, xmax+0.5)
	e = edges(h)
	ec = [(e[i+1] + e[i])/2.0 for i in 1:length(e)-1]
	w = h.weights
	O = w
	ff = fg.g[1]
	E = ff.(ec)
	length(O), length(E)
	num = (O .- E) .* (O .- E)
	denom = E
	ndf = length(e) 
	sum(num ./ denom) /ndf
end

# ╔═╡ Cell order:
# ╠═805955a4-f018-4e68-a68a-1e6b64a99378
# ╠═53dee408-58f8-11ed-2b44-db5e37f3a91e
# ╠═a69771ca-f756-4f84-8bb4-1c044fa4db1f
# ╠═f75fa3dd-d3dc-4f5d-877b-cc9d8f57d95c
# ╠═921fb297-2c3d-4f45-847c-b2f009052e9f
# ╠═3274825d-61a0-4856-9340-97420d237c27
# ╠═082a3d2a-c915-4c98-9221-d4610ae112bc
# ╠═71f427b9-4a44-492c-a826-1bb403fcf7aa
# ╠═52f6a0e2-a5f6-452e-ba45-14958b9ee3a4
# ╠═7690127a-b6ab-4661-aea3-6be93c0f8760
# ╠═519399b4-8f8c-4fcb-bfe3-4f7d7d67e761
# ╠═58670f80-1bf2-4af1-a588-42d52f0394fb
# ╠═f2476981-326f-4f7e-b796-90a9f357ee4a
# ╠═a2baefc2-236c-4b98-a9d1-1506feacaee0
# ╠═de099068-fbe4-428a-89d7-6b6119c52fb9
# ╠═bd496019-51e2-46ee-96e7-048bf53dc547
# ╠═99eb99db-56f7-4884-bb1a-ca362fee1ba0
# ╠═a88457ad-0f1c-4ca8-babf-5eaf3f040de7
# ╠═768d7845-a601-46e4-a86c-09d98d35f715
# ╠═67acec5f-a69e-4af7-98e8-8c8fcdec77b1
# ╠═0020aaac-a114-4cfb-bdba-955b4f0fdabd
# ╠═65bb374d-aa0a-4c97-aab7-1be14120d137
# ╠═db275e86-b644-4982-8666-6eda78ed5072
# ╠═b8cb9917-bb46-43a8-95a1-a05e6c2d0d7f
# ╠═ca6b160a-bbea-4ae0-8ab0-de31c387d439
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
# ╠═a71cbe76-0a65-4cf9-be46-0cc695b55204
# ╠═c7879c6a-c905-4cd0-9661-2cb924a09702
# ╠═c0816fb9-00c7-41c4-b2e4-da4dc50fdbd1
# ╠═88aab6f6-07a8-42c4-b39f-33ed7bb0a75b
# ╠═ee8468ab-2987-497c-83ac-720902561490
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
