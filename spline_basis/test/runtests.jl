using Test, LinearAlgebra, Random, Statistics
include(joinpath(@__DIR__,"..","src","BSpline.jl"))
using .BSplineBasis

@testset "Cubic B-spline values and derivatives" begin
    # Independent polynomial formulas for uniform cubic B-splines on [0,1].
    t = collect(-3.:4.)
    x = [0.2,0.5,0.8]
    B = evaluate_basis(t,x)
    expected = hcat((1 .- x).^3,3 .* x.^3-6 .* x.^2 .+4,
                    -3 .* x.^3+3 .* x.^2+3 .* x .+1,x.^3)/6
    first = hcat(-(1 .- x).^2/2,(9 .* x.^2-12 .* x)/6,
                 (-9 .* x.^2+6 .* x .+3)/6,x.^2/2)
    second = hcat(1 .- x,3 .* x .-2,-3 .* x .+1,x)
    @test B ≈ expected atol=1e-13
    @test evaluate_basis(t,x;derivative=1) ≈ first atol=1e-13
    @test evaluate_basis(t,x;derivative=2) ≈ second atol=1e-13
    @test B*ones(4) ≈ ones(3) atol=1e-13
    @test minimum(B)>=0
    @test all(sum(row .>0)<=4 for row in eachrow(B))
    greville = [mean(t[j+1:j+3]) for j in 1:4]
    @test B*greville ≈ x atol=1e-13
    @test evaluate_basis(t,x;derivative=1)*ones(4) ≈ zeros(3) atol=1e-13
    @test evaluate_basis(t,x;derivative=2)*ones(4) ≈ zeros(3) atol=1e-13
    @test evaluate_basis(t,x;derivative=1)*greville ≈ ones(3) atol=1e-13
    @test evaluate_basis(t,x;derivative=2)*greville ≈ zeros(3) atol=1e-13
    # The package adapter must retain endpoint values, duplicates and row order.
    xe = [1.0,0.2,0.0,0.2,0.8]
    expected_edges = hcat((1 .- xe).^3,3 .* xe.^3-6 .* xe.^2 .+4,
                         -3 .* xe.^3+3 .* xe.^2+3 .* xe .+1,xe.^3)/6
    @test evaluate_basis(t,xe) ≈ expected_edges atol=1e-13
    @test evaluate_basis(t,xe;derivative=1)*greville ≈ ones(length(xe)) atol=1e-13
    @test size(evaluate_basis(t,Float64[])) == (0,4)
end

@testset "Curvature penalty and centering" begin
    breaks = [-2.,-0.9,0.1,0.7,2.]
    t = extend_knots(breaks)
    k = length(t)-4
    x = collect(range(-2,2;length=101))
    c = build_spline(x,t)
    greville = [mean(t[j+1:j+3]) for j in 1:k]
    @test c.S_raw ≈ c.S_raw'
    @test minimum(eigvals(Symmetric(c.S_raw))) > -1e-10
    @test rank(c.S_raw;rtol=1e-10)==k-2
    @test c.S_raw*ones(k) ≈ zeros(k) atol=1e-10
    @test c.S_raw*greville ≈ zeros(k) atol=1e-10
    @test c.H'*c.H ≈ Matrix{Float64}(I,k-1,k-1) atol=1e-12
    @test maximum(abs,vec(mean(c.B;dims=1)))<1e-12
    @test rank(c.S;rtol=1e-10)==k-2
    # Simpson integration independently checks the exact quadratic curvature
    # products integrated by BSplineKit's Gauss quadrature.
    coeff = randn(MersenneTwister(42),k)
    integral = 0.0
    for j in 1:length(breaks)-1
        lo,hi = breaks[j],breaks[j+1]
        d = evaluate_basis(t,[lo,(lo+hi)/2,hi];derivative=2)*coeff
        integral += (hi-lo)/6*(d[1]^2+4d[2]^2+d[3]^2)
    end
    @test dot(coeff,c.S_raw*coeff) ≈ integral rtol=1e-12
    alpha = randn(MersenneTwister(7),k-1)
    @test c.B*alpha ≈ c.B_raw*(c.H*alpha) atol=1e-12
    @test dot(alpha,c.S*alpha) ≈ dot(c.H*alpha,c.S_raw*(c.H*alpha))/c.penalty_scale rtol=1e-12
    @test maximum(abs,evaluate_basis(t,[breaks[1],breaks[end]];derivative=2))>0.1
    for value in breaks[2:end-1]
        @test evaluate_basis(t,[value-1e-9];derivative=2) ≈ evaluate_basis(t,[value+1e-9];derivative=2) atol=1e-6
    end
end

@testset "Linear prediction tails" begin
    t = default_knots([-2.,2.],8)
    lo,hi = t[4],t[end-3]
    for (endpoint,point) in [(lo,lo-0.7),(hi,hi+0.7)]
        value = evaluate_basis(t,[endpoint])
        slope = evaluate_basis(t,[endpoint];derivative=1)
        @test evaluate_basis(t,[point]) ≈ value+(point-endpoint)*slope atol=1e-12
        @test evaluate_basis(t,[point];derivative=1) ≈ slope atol=1e-12
        @test evaluate_basis(t,[point];derivative=2)==zeros(1,8)
    end
end

@testset "Knot settings and invalid inputs" begin
    t = default_knots([-2.,2.],8)
    @test length(t)==12
    @test t[4] ≈ -2.004 atol=1e-14
    @test t[end-3] ≈ 2.004 atol=1e-14
    @test default_knots([-2.,0.,0.,2.],8)==t
    @test length(extend_knots([-2.,0.,2.]))==9
    @test_throws ArgumentError evaluate_basis(collect(1.:7.),[4.])
    @test_throws ArgumentError evaluate_basis([0.,1.,2.,3.,3.,5.,6.,7.],[3.])
    @test_throws ArgumentError evaluate_basis([0.,1.,2.,3.,4.,5.,6.,NaN],[3.])
    @test_throws ArgumentError evaluate_basis(t,[NaN])
    @test_throws ArgumentError evaluate_basis(t,[0.];derivative=3)
    @test_throws ArgumentError default_knots([1.,1.],8)
    @test_throws ArgumentError default_knots([-2.,2.],3)
    @test_throws ArgumentError extend_knots([-2.,0.,0.,2.])
    @test_throws ArgumentError build_spline([-2.,0.,2.],t)
    @test_throws ArgumentError build_spline(collect(range(-3,3;length=100)),t)
end
