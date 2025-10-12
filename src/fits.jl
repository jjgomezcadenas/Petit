using LinearAlgebra
using GLM
using LsqFit
using Distributions
using StatsBase
using Printf
using Plots

# Constants for fit bounds
const FIT_SIGMA_MULTIPLIER_LOWER = 100.0
const FIT_SIGMA_MULTIPLIER_UPPER = 100.0
const FIT_WEIGHT_MULTIPLIER = 100.0

"""
    lfit(ndf::DataFrame)

Perform a linear fit on DataFrame with columns x_mean and y_mean.

# Arguments
- `ndf::DataFrame`: DataFrame with columns x_mean and y_mean

# Returns
- Tuple of (fitted function, predictions, coefficients)
"""
function lfit(ndf::DataFrame)
    lr = lm(@formula(y_mean ~ x_mean), ndf)
    c = coef(lr)
    return x -> c[1] + c[2]*x, predict(lr), c
end

"""
    RFit

Structure to hold polynomial fit results.

# Fields
- `fitpar::Vector{<:Number}`: Fitted parameters
- `fitstd::Vector{<:Number}`: Standard errors of parameters
- `ci::Vector{Tuple{Number, Number}}`: Confidence intervals
- `g::Function`: Fitted function
"""
struct RFit
	fitpar::Vector{<:Number}
	fitstd::Vector{<:Number}
	ci::Vector{Tuple{Number, Number}}
	g::Function
end


"""
    gpol1(ct::Vector{<:Real})::Function
return a degree 1 polynomial function with parameters ct.
"""
function gpol1(ct::Vector{<:Real})::Function
    function f1(z::Real)
        return ct[1] + ct[2] * z
    end
    return f1
end


"""
    gpol2(ct::Vector{<:Real})::Function
return a degree 2 polynomial function with parameters ct.
"""
function gpol2(ct::Vector{<:Real})::Function
    function f2(z::Real)
        return ct[1] + ct[2] * z + ct[3] * z^2
    end
    return f2
end


"""
    gpol2(ct::Vector{<:Real})::Function
return a degree 3 polynomial function with parameters ct.
"""
function gpol3(ct::Vector{<:Real})::Function
    function f3(z::Real)
        return ct[1] + ct[2] * z + ct[3] * z^2 + ct[4] * z^3
    end
    return f3
end


"""
    func1dfit(ffit::Function, x::Vector{<:Real},
              y::Vector{<:Real}, p0::Vector{<:Real},
              lb::Vector{Float64}, ub::Vector{Float64})

Fit a function to the data x, y with start prediction p0
and return coefficients and errors.
"""
function func1dfit(ffit::Function, x::Vector{<:Real},
                   y::Vector{<:Real}, p0::Vector{<:Real},
                   lb::Vector{Float64}, ub::Vector{Float64}; ci::Float64=1.0-0.66)
    fq = curve_fit(ffit, x, y, p0, lower=lb, upper=ub)
    cfq = coef(fq)
    #@info "coef(fq)" cfq
    sfq = stderror(fq)
    #@info "std(fq)" sfq
    #@info "margin_of_error (90%)" margin_error(fq, ci) # default 1 sigma (1-0.66)
    #@info " confidence_interval (90%)" confidence_interval(fq,ci)
    return cfq, sfq, confidence_interval(fq,ci)
end


"""
    func1dfit(ffit::Function, x::Vector{<:Real}, y::Vector{<:Real},
              yerr::Vector{<:Real}, p0::Vector{<:Real},
              lb::Vector{<:Real}, ub::Vector{<:Real})
Fit function ffit to data x,y taking into account the errors
on the y (yerr). Assumes weights of 1/sigma^2 as in standard least squares
and returns fit result.
"""
function func1dfit(ffit::Function, x::Vector{<:Real}, y::Vector{<:Real},
	               yerr::Vector{<:Real}, p0::Vector{<:Real},
                   lb::Vector{<:Real}, ub::Vector{<:Real}; ci::Float64=1.0-0.66)
	fq = curve_fit(ffit, x, y, yerr.^-2, p0, lower=lb, upper=ub)
    return fq
end


"""
    fit_pol1(x::Vector{<:Real}, y::Vector{<:Real}, ci=0.1)
Fit a 1st degree polynomial to the data x, y and return
a RFit object with the result.

# Arguments
- `x::Vector{<:Real}`: X data
- `y::Vector{<:Real}`: Y data
- `ci::Real=0.1`: Confidence interval level (must be between 0 and 1)

# Returns
- `RFit`: Fitted polynomial with parameters, errors, and confidence intervals
"""
function fit_pol1(x::Vector{<:Real}, y::Vector{<:Real}, ci=0.1)
    # Validate inputs
    length(x) != length(y) && throw(ArgumentError("x and y must have same length"))
    length(x) < 2 && throw(ArgumentError("Need at least 2 points for linear fit"))
    !(0 < ci < 1) && throw(ArgumentError("Confidence interval must be between 0 and 1, got $ci"))

    @. pol(x, p) = p[1] + p[2] * x
    p0 = [1.0, 1.0]
    fq = curve_fit(pol, x, y, p0)

    RFit(coef(fq), stderror(fq), confidence_interval(fq, ci),
         gpol1(coef(fq)))
end


"""
    fit_pol2(x::Vector{<:Real}, y::Vector{<:Real}, ci=0.1)
Fit a 2nd degree polynomial to the data x, y and return
a RFit object with the result.

# Arguments
- `x::Vector{<:Real}`: X data
- `y::Vector{<:Real}`: Y data
- `ci::Real=0.1`: Confidence interval level (must be between 0 and 1)

# Returns
- `RFit`: Fitted polynomial with parameters, errors, and confidence intervals
"""
function fit_pol2(x::Vector{<:Real}, y::Vector{<:Real}, ci=0.1)
    # Validate inputs
    length(x) != length(y) && throw(ArgumentError("x and y must have same length"))
    length(x) < 3 && throw(ArgumentError("Need at least 3 points for quadratic fit"))
    !(0 < ci < 1) && throw(ArgumentError("Confidence interval must be between 0 and 1, got $ci"))

    @. pol(x, p) = p[1] + p[2] * x + p[3] * x^2
    p0 = [1.0, 1.0, 1.0]
    fq = curve_fit(pol, x, y, p0)
    RFit(coef(fq), stderror(fq), confidence_interval(fq, ci),
         gpol2(coef(fq)))
end


"""
    fit_pol3(x::Vector{<:Real}, y::Vector{<:Real}, ci=0.1)
Fit a 3rd degree polynomial to the data x, y and return
a RFit object with the result.

# Arguments
- `x::Vector{<:Real}`: X data
- `y::Vector{<:Real}`: Y data
- `ci::Real=0.1`: Confidence interval level (must be between 0 and 1)

# Returns
- `RFit`: Fitted polynomial with parameters, errors, and confidence intervals
"""
function fit_pol3(x::Vector{<:Real}, y::Vector{<:Real}, ci=0.1)
    # Validate inputs
    length(x) != length(y) && throw(ArgumentError("x and y must have same length"))
    length(x) < 4 && throw(ArgumentError("Need at least 4 points for cubic fit"))
    !(0 < ci < 1) && throw(ArgumentError("Confidence interval must be between 0 and 1, got $ci"))

    @. pol(x, p) = p[1] + p[2] * x + p[3] * x^2 + p[4] * x^3
    p0 = [1.0, 1.0, 1.0, 1.0]
    fq = curve_fit(pol, x, y, p0)
    RFit(coef(fq), stderror(fq), confidence_interval(fq, ci),
         gpol3(coef(fq)))
end

"""
    FGauss

Structure to hold Gaussian fit results.

# Fields
- `mu::Vector{<:Real}`: Fitted mean(s)
- `std::Vector{<:Real}`: Fitted standard deviation(s)
- `C::Vector{<:Real}`: Normalization constant(s)
- `δmu::Vector{<:Real}`: Errors on mean(s)
- `δstd::Vector{<:Real}`: Errors on std(s)
- `δC::Vector{<:Real}`: Errors on normalization(s)
- `cimu::Vector{Tuple{Float64, Float64}}`: Confidence intervals for mean(s)
- `cistd::Vector{Tuple{Float64, Float64}}`: Confidence intervals for std(s)
- `ciC::Vector{Tuple{Float64, Float64}}`: Confidence intervals for C(s)
- `h`: Histogram (can be Histogram or Histo1d)
- `X::Vector{<:Real}`: X values
- `Y::Vector{<:Real}`: Y values (fitted)
- `g::Vector{Function}`: Fitted function(s)
"""
struct FGauss
	mu ::Vector{<:Real}
	std::Vector{<:Real}
	C  ::Vector{<:Real}
    δmu ::Vector{<:Real}
	δstd::Vector{<:Real}
	δC  ::Vector{<:Real}
    cimu ::Vector{Tuple{Float64, Float64}}
	cistd::Vector{Tuple{Float64, Float64}}
	ciC  ::Vector{Tuple{Float64, Float64}}
	h  ::Any  # Can be Histogram or Histo1d
	X  ::Vector{<:Real}
	Y  ::Vector{<:Real}
	g  ::Vector{Function}
end


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


"""
    gausg2(μ1::Real, σ1::Real, C1::Real, μ2::Real, σ2::Real, C2::Real)
Return a function for the sum of two Gaussians with parameters
    μ1, σ1 and μ2, σ2 and normalisations C1 and C2.
"""
function gausg2(μ1::Real, σ1::Real, C1::Real, μ2::Real, σ2::Real, C2::Real)::Function
	function gausx(x::Real)
		return C1 * pdf(Normal(μ1, σ1,), x) + C2 * pdf(Normal(μ2, σ2,), x)
	end
	return gausx
end

@. gauss1(x, p)   = p[1] * pdf(Normal(p[2], p[3]), x)
@. gauss2(x, p)   = p[1] * pdf(Normal(p[2], p[3]), x) + p[4] * pdf(Normal(p[5], p[6]), x)
@. gauss2cm(x, p) = p[1] * pdf(Normal(p[2], p[3]), x) + p[4] * pdf(Normal(p[2], p[5]), x)


"""
    gaussfm(mu::Real)::Function
Return a Gaussian Fitting function with a fixed mean mu.
"""
function gaussfm(mu::Real)::Function
	function gauss(x::Vector{<:Real}, p::Vector{<:Real})
		return @. p[1]* pdf(Normal(mu, p[2]), x)
	end
	return gauss
end


"""
    gauss2fm(mu::Real)::Function
Return a function of the sum of two Gaussians with a fixed common mean mu.
"""
function gauss2fm(mu::Real)::Function
	function gauss2(x::Vector{<:Real}, p::Vector{<:Real})
		return @. p[1]* pdf(Normal(mu, p[2]), x) + p[3]* pdf(Normal(mu, p[4]), x)
	end
	return gauss2
end


"""
	fit_gauss(h::Histo1d)

Fit a normal distribution to a Histo1d and return FGauss object.
"""
function fit_gauss(h::Histo1d)
	c = h.centers
	w = eltype(c).(h.weights)
	@debug "histo"  w c
	mu, sigma = mean_and_std(c, Weights(h.weights); corrected = false)
	@debug "mu, std" mu, sigma

	# fit parameters lb, ub, po are lower, upper bounds and pars
    lb = [         0.0, mu - FIT_SIGMA_MULTIPLIER_LOWER * sigma, sigma / FIT_SIGMA_MULTIPLIER_LOWER]
    ub = [FIT_WEIGHT_MULTIPLIER * sum(w), mu + FIT_SIGMA_MULTIPLIER_UPPER * sigma, FIT_SIGMA_MULTIPLIER_UPPER * sigma]
    p0 = [      sum(w), mu                ,         sigma]
	(CC, μ, σ), (δCC, δμ, δσ),(ciCC, ciμ, ciσ) = func1dfit(gauss1, c, w, p0, lb, ub)
	gx = gausg(μ, σ, CC)
	return FGauss([μ], [σ], [CC], [δμ], [δσ], [δCC], [ciμ], [ciσ], [ciCC], h, c, gx.(c), [gx])

end


"""
	fit_gauss(x::Vector{Float64}, xmin::Float64, xmax::Float64;
              bins::Integer=50, norm=false)

Fit a normal distribution to data x,y and return FGauss object.
"""
fit_gauss(x::Vector{Float64}, xmin::Float64, xmax::Float64;
	      bins::Integer=50, norm=false) = fit_gauss(hist1d(x; nbins=bins, xlim=(xmin, xmax), norm=norm))


"""
	fit_gauss_fm(x::Vector{Float64}, xmin::Float64, xmax::Float64, bins=50, fm=0.0)

Fit a gaussian with a fixed mean fm and return an FGauss object.
"""
function fit_gauss_fm(x::Vector{Float64}, xmin::Float64, xmax::Float64;
	                  bins=50, norm=false, fm=0.0)

	# fit the unbinned distribution
	xx =  in_range(x, xmin, xmax)
	σ  = std(xx)
	@debug "gfit_gauss_fm: σ = $σ"

	# bin distribution
    h = hist1d(xx; nbins=bins, xlim=(xmin, xmax), norm=norm)
    c = centers(h)
    w = eltype(c).(h.weights)
    @debug "histo w and c"  w c

	# fit parameters lb, ub, po are lower, upper bounds and pars
    lb = [            0.0, σ / 10.0]
    ub = [1.0e+6 * sum(w), σ * 10.0]
    p0 = [         sum(w), σ       ]

	g1 = gaussfm(fm)
	(CC, sigma), (δCC, δσ),(ciCC, ciσ) = func1dfit(g1, c, w, p0, lb, ub)

    
	@debug "CC,  sigma"  CC  sigma
	gx = gausg(fm, sigma, CC)

	return FGauss([fm], [sigma], [CC], [0.0], [δσ], [δCC], 
                  [(0.0,0.0)], [ciσ], [ciCC], h, c, gx.(c), [gx])
end


"""
    plot_fit_gauss(x::Vector{Float64}, xs::String, ys::String,
                   bins::Integer, xmin::Float64, xmax::Float64;
                   xgmin::Float64, xgmax::Float64, gbins::Integer=50)
Fit a Gaussian to the histogram of data x between xmin and xmax with bins
and plot the results.
"""
function plot_fit_gauss(x::Vector{Float64}, xs::String, ys::String,
                        bins::Integer, xmin::Float64, xmax::Float64;
                        xgmin::Float64, xgmax::Float64, gbins::Integer=50)

    # Create histogram with proper keyword arguments
    h = hist1d(x; nbins=bins, xlim=(xmin, xmax), norm=false)
    fg   = fit_gauss(x, xgmin, xgmax, bins=gbins, norm=false)
    gx   = fg.g[1]
    X    = centers(h)
    Y    = h.weights
    σY   = sqrt.(Y)
    lbl  = @sprintf " μ=%5.1f, σ =%5.1f " fg.mu[1] fg.std[1]
    lbl  = string("gaussian fit:\n", lbl)
    p    = scatter(X, Y, yerror=σY,fmt=:png, shape=:circle, color=:black, label="data", legend=true)
    p    = plot!(p, X, gx.(X), lw=2, label=lbl, legend=true, fmt = :png)
    xlabel!(xs)
    ylabel!(ys)  # Fixed: was xlabel!(ys)
    return fg, p
end


"""
    fitg1(x::Vector{<:Real}, xs::String, nbins::Int64, xmin::Real, xmax::Real;
	      xgmin::Real, xgmax::Real, fbins::Int64=100, norm::Bool=true,
          fm::Real=0.0, flex_mean=false)

Histograms data x and returns a plot and a Gaussian fit to the histogram.
Option flex_mean allows the mean to float or be fixed to fm.
"""
function fitg1(x::Vector{<:Real}, xs::String, nbins::Int64, xmin::Real, xmax::Real;
	           xgmin::Real, xgmax::Real, fbins::Int64=100, norm::Bool=true,
               fm::Real=0.0, flex_mean=false)

	h, p = step_hist(x, nbins, xs; norm=norm, legend=true)
    if flex_mean
        fg = fit_gauss(x, xgmin, xgmax, bins=fbins, norm=norm)
    else
        fg = fit_gauss_fm(x, xgmin, xgmax, bins=fbins, norm=norm, fm=fm)
    end
	gx  = fg.g[1]
	X   = centers(h)
	lbl = @sprintf "σ =%4.1f " fg.std[1]
    p   = plot!(p, X, gx.(X), lw=2, label=lbl, legend=true, fmt = :png)
	xlabel!(xs)
	return fg, p
end


"""
	gfit_gauss2_cmean(y, xmin, xmax, bins, sigmas, cs, cmean=0.0)

Fit a double gaussian (with sigmas -->[sigma1, sigma2] cs -->[c1, c2] )
and a fixed mean, cmean, to data.
"""
function gfit_gauss2_cmean(y::Vector{Float64}, xmin::Float64, xmax::Float64,
	                       bins::Integer, sigmas::Vector{Float64}, cs::Vector{Float64},
						   norm=false, cmean=0.0)

    x = in_range(y, xmin, xmax)
    h = hist1d(x; nbins=bins, xlim=(xmin, xmax), norm=norm)
    c = centers(h)
    w = h.weights
    @debug "histo centers and weights in full region"  w c

    g2 = gauss2fm(cmean)
    # fit parameters lb, ub, po are lower, upper bounds and pars

    lb = [cs[1]/100.0, sigmas[1]/5.0, cs[2]/100.0, sigmas[2]/5.0]
    ub = [cs[1]*100.0, sigmas[1]*5.0, cs[2]*100.0, sigmas[2]*5.0]
    p0_bounds = [cs[1], sigmas[1], cs[2], sigmas[2]]

    @debug "pars" p0 lb ub
    # fit double gaussian
    fq = curve_fit(g2, c, w, p0_bounds, lower=lb, upper=ub)
    C1, sigma1, C2,  sigma2   = coef(fq)
    @debug "C1 sigma1 C2 sigma2" C1 sigma1 C2  sigma2

    #
    gx  = gausg2(cmean, sigma1, C1, cmean, sigma2, C2)
    gx1 = gausg( cmean, sigma1, C1)
    gx2 = gausg( cmean, sigma2, C2)

	return FGauss([cmean, cmean], [sigma1, sigma2], [C1, C2],
	               [0.0, 0.0], [0.0, 0.0], [0.0, 0.0],
	               [(0.0,0.0), (0.0,0.0)], [(0.0,0.0), (0.0,0.0)], [(0.0,0.0), (0.0,0.0)],
			       h, c, gx.(c), [gx, gx1, gx2])
end

"""
	fit_2gauss_cmean(data, gp, g1p, g2p, cm)

Fit two gaussian with common mean JJ Review!!!
"""
function fit_2gauss_cmean(data, gp, g1p, g2p, cm, norm=false)
    gf1 = fit_gauss_fm(data, g1p.xmin,g1p.xmax, bins=g1p.nbin, norm=norm, fm=cm)
    @debug gf1
    gf2 = fit_gauss_fm(data, g2p.xmin,g2p.xmax,bins=g2p.nbin, norm= norm, fm=cm)
    @debug gf2
    gf = gfit_gauss2_cmean(data, gp.xmin,gp.xmax,gp.nbin,
	                       [gf1.std[1], gf2.std[1]], [gf1.C[1], gf2.C[1]], norm)
    @debug gf
    return gf
end


"""
    fitg2(x, xs, xmin, xmax, xg1min, xg1max, xg2min, xg2max, xgmin, xgmax; bins=100)
Fits 2 gaussians with common mean (0 by default) to vector x. JJ Review!!
"""
function fitg2(x::Vector{<:Real}, xs::String, bins::Int64, xmin::Real, xmax::Real;
	           xg1min::Real, xg1max::Real, xg2min::Real, xg2max::Real, xgmin::Real, xgmax::Real, cm=0.0,
         	   g1bins::Int64=100, g2bins::Int64=100, gbins::Int64=100, norm::Bool=true)

	_, p = step_hist(x, bins, xs; norm=norm, legend=true)

    g1p = (xmin = xg1min, xmax = xg1max, nbin=g1bins)
    g2p = (xmin = xg2min, xmax = xg2max, nbin=g2bins)
    gp  = (xmin = xgmin , xmax = xgmax , nbin=gbins )

    fg  = fit_2gauss_cmean(x, gp, g1p, g2p, cm, norm)
	gx  = fg.g[1]
	gx1 = fg.g[2]
	gx2 = fg.g[3]

	lbl = @sprintf "σt =%4.1f mm, σ =%4.1f mm" fg.std[1] fg.std[2]
	st  = @sprintf "σt =%4.1f mm " fg.std[1]
	sf  = @sprintf "σ  =%4.1f mm" fg.std[2]
    p = plot!(p, fg.X, fg.Y, label=lbl, lw=2, fmt = :png)
    p = plot!(p, fg.X, gx1.(fg.X), label=st, lw=1, fmt = :png)
    p = plot!(p, fg.X, gx2.(fg.X), label=sf, lw=1, fmt = :png)
    #xlabel!(xs)
    return fg, p
end


"""
	fit_profile(x1, x2, tx1, ty1, fit="pol1", bins=25)
    fit_profile(df1, c1, c2, tx1, ty1, fit="pol1", bins=25)
Create and fit a profile with pol1 or poli2 functions.
Return fit parameters, fit function and plot
"""
fit_profile(df1::DataFrame, c1::String, c2::String,
            tx1::String, ty1::String, fit="pol1", bins=25;
            ybin_width::Real=0.1, ymin::Real=352.4, ymax::Real=392.4,
            min_proportion::Real=0.0) =
	fit_profile(df1[!,c1], df1[!,c2], tx1, ty1, fit, bins,
                ybin_width=ybin_width, ymin=ymin, ymax=ymax, min_proportion=min_proportion)


function fit_profile(x1::Vector{<:Real}, x2::Vector{<:Real},
	                 tx1::String, ty1::String, fit="pol1", bins=25;
                     ybin_width::Real=0.1, ymin::Real=minimum(x2),
                     ymax::Real=maximum(x2), min_proportion::Real=0.0)

    pdf1, _ = p1df(x1,x2, bins, ybin_width=ybin_width, ymin=ymin, ymax=ymax, min_proportion=min_proportion)

    if fit == "pol1"
        fr = fit_pol1(pdf1.x_mean, pdf1.y_mean)
    elseif fit == "pol2"
        fr = fit_pol2(pdf1.x_mean, pdf1.y_mean)
	elseif fit == "pol3"
		fr = fit_pol3(pdf1.x_mean, pdf1.y_mean)
	else
		println("option not implemented")
		return nothing
    end

    p1 = scatter(pdf1.x_mean, pdf1.y_mean, yerror=pdf1.y_std, fmt=:png,
	             shape=:circle, color=:black, legend=false)
    p1 = plot!(p1, pdf1.x_mean, fr.g.(pdf1.x_mean), fmt=:png)
    xlabel!(tx1)
    ylabel!(ty1)

   return fr, p1
end
