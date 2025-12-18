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
	using BayesHistogram
	import Glob
end


# ╔═╡ a69771ca-f756-4f84-8bb4-1c044fa4db1f
using ATools

# ╔═╡ f75fa3dd-d3dc-4f5d-877b-cc9d8f57d95c
PlutoUI.TableOfContents(title="Toy MC", indent=true)

# ╔═╡ 921fb297-2c3d-4f45-847c-b2f009052e9f
md"""
# Generate a gaussian
"""

# ╔═╡ 4350f088-c0a8-4a6e-9246-a574a97e106d
begin
	gaus = Normal(5000.0, 80.0)
	ng=rand(gaus,1000000)
end

# ╔═╡ 52f6a0e2-a5f6-452e-ba45-14958b9ee3a4
begin
	hmc1, pmc1 = hist1d(ng, "events", 100, 4500.0, 5500.0)
	plot(pmc1)
end

# ╔═╡ bd496019-51e2-46ee-96e7-048bf53dc547
md""" Plane 2: select min value to fit: $(@bind rmin NumberField(4500.0:5550.0, default=4600.0))"""

# ╔═╡ 99eb99db-56f7-4884-bb1a-ca362fee1ba0
md""" Plane 2: select max value to fit: $(@bind rmax NumberField(4500.0:5500.0, default=5400.0))"""

# ╔═╡ 7aab66b9-1796-41e0-a800-d188ac8bc944
tgaus = truncated(gaus; lower=rmin, upper=rmax)

# ╔═╡ 1b9a7002-8ad6-42aa-a691-620923daedf0
md"""
### ML unbinned fit 
"""

# ╔═╡ 24443046-a148-445c-becd-8c58117ff9b6
fdset = filter(x->x>rmin && x<rmax, ng)

# ╔═╡ 360e8386-2b85-4f8a-913f-2c3dea59a4d2
begin
	gaus2 = Normal(5000.0, 520.0)
	ng2=rand(gaus,length(fdset))
end

# ╔═╡ 0bf4b0a3-90bd-4152-8e9e-d218451a7435
begin
	hmcx1, pmcx1 = hist1d(fdset, "events", 100, 4500.0, 5500.0)
	hmcx2, pmcx2 = hist1d(ng2, "events2", 100, 4500.0, 5500.0)
	plot!(pmcx1,pmcx2)
	
end

# ╔═╡ 8029c54f-c45d-4c1a-9815-7ccf8654ec24
length(ng2)

# ╔═╡ 7c4bf607-3d3c-4f70-9ffb-fe601ac92487
md"""
### Binned fit 
"""

# ╔═╡ c0816fb9-00c7-41c4-b2e4-da4dc50fdbd1
md"""
# Functions
"""

# ╔═╡ 9659ec6d-6eee-46df-b19a-940db87bb7f3
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
	#fmle = fit_mle(Normal, fdset)
	#fmle = fit(Normal, fdset)
	tgaus = truncated(Normal; lower=xmin, upper=xmax)
	fmle = fit_mle(tgaus, fdset)
	fmle, -loglikelihood(fmle, fdset)
end

# ╔═╡ e16ea9bc-7aa6-47a6-9a86-8d0f656dc002
ftmle = norm_mle_fit(ng, rmin, rmax)

# ╔═╡ 41ec942b-4c04-4b43-8475-d4a14fae02c2
"""
Vary the paramters of the fit around the best value. Return a named
tuple (lft object) with information to compute parabolic fit
"""
function norm_mle_var(x::Vector{Float64}, xmin::Float64, xmax::Float64; var="σ",
                      sr=2.0, ss=100.0)
	
	fdset = filter(x->x>xmin && x<xmax, x)
	tgaus = Truncated(Normal, xmin, xmax)
	fmle = fit_mle(tgaus, fdset)
	#fmle = fit_mle(Normal, fdset)
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
lft = norm_mle_var(ng, rmin, rmax, var="σ", sr=0.15)

# ╔═╡ b72281dd-e4b1-4c0f-ab1d-8d9df19aa837
begin
	plot(lft.vars, lft.logs, lw=2, label="")
	scatter!([lft.fbest.σ], [lft.llbest], label="best fit", legend=true)
	hline!([lft.llbest + 0.5], label="")
	xlabel!("σ")
	ylabel!("-log(L)")
end

# ╔═╡ fe759a1a-4ba2-4b55-ad56-d9859e6f8d3a
lft2 = norm_mle_var(ng, rmin, rmax, var="μ", sr=0.15)

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

# ╔═╡ 284c283a-48e5-4551-a60b-5fc1dcf0842f
lftsigma.var * length(ng) /(length(ng2))

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
fg, pg = fitg1(ng,"E (LSB)", 100, rmin, rmax;
	           xgmin=rmin, xgmax=rmax, fbins=100, norm=true,
               fm=lftmu.var, flex_mean=true)

# ╔═╡ 026d5f3e-0e02-4692-8cbd-3af88a61f603
plot(pg, legend=false)

# ╔═╡ e30a4fcf-b7ba-4c7a-b9ce-c455c41033d6
fg2, pg2 = fitg1(ng2,"E (LSB)", 100, rmin, rmax;
	           xgmin=rmin, xgmax=rmax, fbins=100, norm=true,
               fm=lftmu.var, flex_mean=true)

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

# ╔═╡ Cell order:
# ╠═805955a4-f018-4e68-a68a-1e6b64a99378
# ╠═53dee408-58f8-11ed-2b44-db5e37f3a91e
# ╠═a69771ca-f756-4f84-8bb4-1c044fa4db1f
# ╠═f75fa3dd-d3dc-4f5d-877b-cc9d8f57d95c
# ╠═921fb297-2c3d-4f45-847c-b2f009052e9f
# ╠═4350f088-c0a8-4a6e-9246-a574a97e106d
# ╠═7aab66b9-1796-41e0-a800-d188ac8bc944
# ╠═52f6a0e2-a5f6-452e-ba45-14958b9ee3a4
# ╠═bd496019-51e2-46ee-96e7-048bf53dc547
# ╠═99eb99db-56f7-4884-bb1a-ca362fee1ba0
# ╠═1b9a7002-8ad6-42aa-a691-620923daedf0
# ╠═e16ea9bc-7aa6-47a6-9a86-8d0f656dc002
# ╠═24443046-a148-445c-becd-8c58117ff9b6
# ╠═360e8386-2b85-4f8a-913f-2c3dea59a4d2
# ╠═0bf4b0a3-90bd-4152-8e9e-d218451a7435
# ╠═328768e0-c381-4060-b547-282ad0ae43dc
# ╠═b72281dd-e4b1-4c0f-ab1d-8d9df19aa837
# ╠═25f704f9-f280-4f97-aeba-2ad3856c5c4e
# ╠═284c283a-48e5-4551-a60b-5fc1dcf0842f
# ╠═8029c54f-c45d-4c1a-9815-7ccf8654ec24
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
# ╠═e30a4fcf-b7ba-4c7a-b9ce-c455c41033d6
# ╠═c0816fb9-00c7-41c4-b2e4-da4dc50fdbd1
# ╠═9659ec6d-6eee-46df-b19a-940db87bb7f3
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
