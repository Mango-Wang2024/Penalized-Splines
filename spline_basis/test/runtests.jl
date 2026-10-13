using Test, LinearAlgebra, Random, Statistics
include(joinpath(@__DIR__,"..","src","NaturalCubic.jl"))
using .NaturalCubic

@testset "Natural cubic construction" begin
    t = [-2.,-1.4,-0.3,0.2,1.2,2.]
    F = second_derivative_map(t)
    x = collect(range(-2.7,2.8;length=101))
    B = evaluate_basis(t,F,x)
    @test evaluate_basis(t,F,t) ≈ Matrix{Float64}(I,6,6) atol=1e-12
    @test B*ones(6) ≈ ones(length(x)) atol=1e-12
    @test B*t ≈ x atol=1e-12
    @test F*ones(6) ≈ zeros(6) atol=1e-12
    @test F*t ≈ zeros(6) atol=1e-12
    @test F[[1,end],:] == zeros(2,6)
    @test evaluate_basis(t,F,[-3.,3.];derivative=2) == zeros(2,6)
    for j in 2:5
        left=evaluate_basis(t,F,[t[j]-1e-9];derivative=1)
        right=evaluate_basis(t,F,[t[j]+1e-9];derivative=1)
        @test maximum(abs,left-right) < 1e-7
    end
end

@testset "Penalty and centering" begin
    t = [-2.,-1.4,-0.3,0.2,1.2,2.]
    x = collect(range(-2,2;length=101))
    c = build_spline(x,t)
    @test c.S_raw ≈ c.S_raw'
    @test minimum(eigvals(Symmetric(c.S_raw))) > -1e-10
    @test rank(c.S_raw;rtol=1e-10) == 4
    @test c.S_raw*ones(6) ≈ zeros(6) atol=1e-10
    @test c.S_raw*t ≈ zeros(6) atol=1e-10
    @test c.H'*c.H ≈ Matrix{Float64}(I,5,5) atol=1e-12
    @test maximum(abs,vec(mean(c.B;dims=1))) < 1e-12
    @test rank(c.S;rtol=1e-10) == 4
    # Two-point Gauss quadrature exactly integrates squared linear curvature.
    # This independently checks the penalty's numerical value.
    z = randn(MersenneTwister(42),6)
    integral = 0.0
    for j in 1:5
        h=t[j+1]-t[j]
        mid=(t[j+1]+t[j])/2
        points=[mid-h/(2sqrt(3)),mid+h/(2sqrt(3))]
        curvature=evaluate_basis(t,c.F,points;derivative=2)*z
        integral += h/2*sum(abs2,curvature)
    end
    @test dot(z,c.S_raw*z) ≈ integral rtol=1e-12
    alpha=randn(MersenneTwister(7),5)
    @test c.B*alpha ≈ c.B_raw*(c.H*alpha) atol=1e-12
    @test dot(alpha,c.S*alpha) ≈ dot(c.H*alpha,c.S_raw*(c.H*alpha))/c.penalty_scale rtol=1e-12
end

@testset "Default knots and invalid inputs" begin
    @test default_knots([3.,0.,1.,2.,2.],3) == [0.,1.5,3.]
    @test_throws ArgumentError second_derivative_map([0.,1.])
    @test_throws ArgumentError second_derivative_map([0.,1.,1.])
    @test_throws ArgumentError second_derivative_map([0.,NaN,2.])
    @test_throws ArgumentError default_knots([0.,1.],3)
    @test_throws ArgumentError build_spline([0.,0.,1.],[0.,0.5,1.])
end
