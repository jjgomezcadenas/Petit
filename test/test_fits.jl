"""
Comprehensive tests for fits.jl functions
"""

using Test
using DataFrames
using Distributions
using StatsBase
using Random

# Load the Petit module
push!(LOAD_PATH, dirname(dirname(@__FILE__)))
using Petit

# Import specific functions we need for testing
import Petit: fit_pol1, fit_pol2, fit_pol3, gpol1, gpol2, gpol3
import Petit: fit_gauss, fit_gauss_fm, gausg, gausg2
import Petit: lfit, RFit, func1dfit

@testset "Polynomial Functions" begin
    @testset "gpol1" begin
        # Test linear polynomial
        f = gpol1([1.0, 2.0])
        @test f(0.0) ≈ 1.0
        @test f(1.0) ≈ 3.0
        @test f(2.0) ≈ 5.0
    end

    @testset "gpol2" begin
        # Test quadratic polynomial
        f = gpol2([1.0, 2.0, 3.0])
        @test f(0.0) ≈ 1.0
        @test f(1.0) ≈ 6.0  # 1 + 2*1 + 3*1^2
        @test f(2.0) ≈ 17.0  # 1 + 2*2 + 3*4
    end

    @testset "gpol3" begin
        # Test cubic polynomial
        f = gpol3([1.0, 2.0, 3.0, 4.0])
        @test f(0.0) ≈ 1.0
        @test f(1.0) ≈ 10.0  # 1 + 2*1 + 3*1 + 4*1
        @test f(2.0) ≈ 49.0  # 1 + 2*2 + 3*4 + 4*8
    end
end

@testset "Polynomial Fitting" begin
    @testset "fit_pol1" begin
        # Generate linear data with small noise
        Random.seed!(42)
        x = collect(0.0:0.1:10.0)
        y_true = 2.0 .+ 3.0 .* x
        y = y_true .+ 0.1 .* randn(length(x))

        result = fit_pol1(x, y)

        @test length(result.fitpar) == 2
        @test result.fitpar[1] ≈ 2.0 atol=0.5  # intercept
        @test result.fitpar[2] ≈ 3.0 atol=0.5  # slope
        @test length(result.fitstd) == 2
        @test length(result.ci) == 2
        @test isa(result.g, Function)

        # Test the fitted function
        @test result.g(0.0) ≈ result.fitpar[1]
        @test result.g(1.0) ≈ result.fitpar[1] + result.fitpar[2]
    end

    @testset "fit_pol2" begin
        # Generate quadratic data
        Random.seed!(42)
        x = collect(0.0:0.1:10.0)
        y_true = 1.0 .+ 2.0 .* x .+ 0.5 .* x.^2
        y = y_true .+ 0.1 .* randn(length(x))

        result = fit_pol2(x, y)

        @test length(result.fitpar) == 3
        @test result.fitpar[1] ≈ 1.0 atol=0.5
        @test result.fitpar[2] ≈ 2.0 atol=0.5
        @test result.fitpar[3] ≈ 0.5 atol=0.5
    end

    @testset "fit_pol3" begin
        # Generate cubic data
        Random.seed!(42)
        x = collect(0.0:0.1:10.0)
        y_true = 1.0 .+ 2.0 .* x .+ 0.5 .* x.^2 .+ 0.1 .* x.^3
        y = y_true .+ 0.1 .* randn(length(x))

        result = fit_pol3(x, y)

        @test length(result.fitpar) == 4
        @test result.fitpar[1] ≈ 1.0 atol=0.5
        @test result.fitpar[2] ≈ 2.0 atol=0.5
        @test result.fitpar[3] ≈ 0.5 atol=0.5
        @test result.fitpar[4] ≈ 0.1 atol=0.5
    end
end

@testset "Gaussian Functions" begin
    @testset "gausg" begin
        # Test single Gaussian
        μ, σ, C = 5.0, 1.0, 100.0
        g = gausg(μ, σ, C)

        # Test at mean (maximum of Gaussian)
        @test g(μ) ≈ C * pdf(Normal(μ, σ), μ)

        # Test symmetry
        @test g(μ + 1.0) ≈ g(μ - 1.0) atol=1e-10

        # Test that it's a valid function
        @test g(μ + σ) > 0.0
        @test g(μ - σ) > 0.0
    end

    @testset "gausg2" begin
        # Test double Gaussian
        μ1, σ1, C1 = 3.0, 0.5, 50.0
        μ2, σ2, C2 = 7.0, 1.0, 100.0
        g = gausg2(μ1, σ1, C1, μ2, σ2, C2)

        # Test that it's sum of two Gaussians
        g1 = gausg(μ1, σ1, C1)
        g2 = gausg(μ2, σ2, C2)

        test_x = 5.0
        @test g(test_x) ≈ g1(test_x) + g2(test_x)
    end
end

@testset "Gaussian Fitting" begin
    @testset "fit_gauss - basic" begin
        # Generate Gaussian data
        Random.seed!(42)
        μ_true, σ_true = 10.0, 2.0
        n_samples = 10000
        data = μ_true .+ σ_true .* randn(n_samples)

        # Fit Gaussian
        xmin, xmax = 0.0, 20.0
        result = fit_gauss(data, xmin, xmax, bins=100)

        # Check results
        @test length(result.mu) == 1
        @test length(result.std) == 1
        @test length(result.C) == 1
        @test result.mu[1] ≈ μ_true atol=0.5
        @test result.std[1] ≈ σ_true atol=0.5
        @test result.C[1] > 0.0

        # Check that errors are reasonable
        @test result.δmu[1] > 0.0
        @test result.δstd[1] > 0.0
        @test result.δC[1] > 0.0

        # Check confidence intervals
        @test length(result.cimu) == 1
        @test length(result.cistd) == 1
        @test length(result.ciC) == 1
    end

    @testset "fit_gauss_fm - fixed mean" begin
        # Generate Gaussian data with known mean
        Random.seed!(42)
        μ_true, σ_true = 0.0, 2.0
        n_samples = 10000
        data = μ_true .+ σ_true .* randn(n_samples)

        # Fit Gaussian with fixed mean
        xmin, xmax = -10.0, 10.0
        result = fit_gauss_fm(data, xmin, xmax, bins=100, fm=μ_true)

        # Check results
        @test result.mu[1] ≈ μ_true  # Should be exactly the fixed mean
        @test result.std[1] ≈ σ_true atol=0.5
        @test result.δmu[1] == 0.0  # Error should be 0 for fixed parameter
    end
end

@testset "Linear Fit" begin
    @testset "lfit" begin
        # Create test DataFrame
        x_data = [1.0, 2.0, 3.0, 4.0, 5.0]
        y_data = [2.0, 4.0, 6.0, 8.0, 10.0]  # Perfect linear relationship y = 2x

        df = DataFrame(x_mean = x_data, y_mean = y_data)

        # Perform fit
        f, pred, coef = lfit(df)

        # Check coefficients
        @test coef[1] ≈ 0.0 atol=1e-10  # intercept
        @test coef[2] ≈ 2.0 atol=1e-10  # slope

        # Check fitted function
        @test f(1.0) ≈ 2.0 atol=1e-10
        @test f(5.0) ≈ 10.0 atol=1e-10

        # Check predictions
        @test length(pred) == length(x_data)
    end
end

@testset "Function Fitting - func1dfit" begin
    @testset "func1dfit without errors" begin
        # Test simple exponential fit
        Random.seed!(42)
        x = collect(0.0:0.1:5.0)
        # Exponential decay: y = A * exp(-λ * x)
        A_true, λ_true = 10.0, 0.5
        y_true = A_true .* exp.(-λ_true .* x)
        y = y_true .+ 0.1 .* randn(length(x))

        # Define fit function
        exp_model(x, p) = @. p[1] * exp(-p[2] * x)

        p0 = [10.0, 0.5]
        lb = [0.0, 0.0]
        ub = [100.0, 10.0]

        params, errors, ci = func1dfit(exp_model, x, y, p0, lb, ub)

        @test length(params) == 2
        @test params[1] ≈ A_true atol=1.0
        @test params[2] ≈ λ_true atol=0.2
        @test length(errors) == 2
        @test length(ci) == 2
    end
end

@testset "Edge Cases and Error Handling" begin
    @testset "Empty data" begin
        # Test with empty vectors
        @test_throws Exception fit_pol1(Float64[], Float64[])
    end

    @testset "Single point" begin
        # Cannot fit with single point
        @test_throws Exception fit_pol1([1.0], [2.0])
    end

    @testset "Insufficient points for polynomial degree" begin
        # Need at least 3 points for pol2
        x = [1.0, 2.0]
        y = [1.0, 2.0]
        @test_throws Exception fit_pol2(x, y)
    end

    @testset "Invalid confidence interval" begin
        # CI should be between 0 and 1
        x = collect(0.0:0.1:10.0)
        y = 2.0 .+ 3.0 .* x

        @test_throws Exception fit_pol1(x, y, -0.1)  # negative CI
        @test_throws Exception fit_pol1(x, y, 1.5)   # CI > 1
    end
end

@testset "RFit Structure" begin
    @testset "RFit construction" begin
        # Test RFit structure
        params = [1.0, 2.0]
        stds = [0.1, 0.2]
        cis = [(0.8, 1.2), (1.6, 2.4)]
        func = x -> 1.0 + 2.0*x

        rfit = RFit(params, stds, cis, func)

        @test rfit.fitpar == params
        @test rfit.fitstd == stds
        @test rfit.ci == cis
        @test rfit.g(0.0) ≈ 1.0
        @test rfit.g(1.0) ≈ 3.0
    end
end

println("\n" * "="^60)
println("All fits.jl tests completed!")
println("="^60)
