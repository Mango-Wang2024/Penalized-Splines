using Test, LinearAlgebra, Random
include(joinpath(@__DIR__, "..", "src", "ModelInputs.jl"))
using .ModelInputs

@testset "Independent diagonal example" begin
    x = collect(range(-1,1; length=30))
    B = hcat(ones(30), x, x.^2, x.^3)
    S = Matrix(Diagonal([0.,0.,4.,9.]))
    c = prepare_inputs(B,S; expected_rank=2)
    @test size(c.X0) == (30,2)
    @test size(c.Z) == (30,2)
    @test c.penalty_rank == 2
    @test c.null_rank == 2
    # Expected answers derived directly from this example, without eigendecomposition.
    expected_covariance = B[:,3:4] * Diagonal([1/4,1/9]) * B[:,3:4]'
    @test relative_error(c.Z*c.Z', expected_covariance) < 1e-12
    @test relative_error(space_projector(c.X0),space_projector(B[:,1:2])) < 1e-12
    a = [0.3,-0.4,1.2,-0.8] # random coefficients first, unpenalized coefficients last
    z = c.C*a
    @test relative_error(B*z,hcat(c.Z,c.X0)*a) < 1e-12
    @test dot(z,S*z) ≈ sum(abs2,a[1:2]) atol=1e-12
    @test rank(c.C) == 4
end

@testset "Rotated coefficient coordinates" begin
    rng = MersenneTwister(42)
    Q = Matrix(qr(randn(rng,5,5)).Q)
    B = randn(rng,40,5)
    d = [0.,0.,2.,7.,11.]
    S = Q*Diagonal(d)*Q'
    c = prepare_inputs(B,S)
    expected_covariance = B*Q[:,3:5]*Diagonal(1 ./ d[3:5])*Q[:,3:5]'*B'
    @test c.penalty_rank == 3
    @test relative_error(c.Z*c.Z',expected_covariance) < 1e-11
    @test relative_error(space_projector(c.X0),space_projector(B*Q[:,1:2])) < 1e-11
    @test relative_error(space_projector(B),space_projector(hcat(c.Z,c.X0))) < 1e-11
    for _ in 1:10
        a = randn(rng,5)
        z = c.C*a
        @test relative_error(B*z,hcat(c.Z,c.X0)*a) < 1e-11
        @test isapprox(dot(z,S*z),sum(abs2,a[1:3]); atol=1e-10,rtol=1e-10)
    end
end

@testset "Unpenalized and fully penalized limits" begin
    B = Matrix{Float64}(I,3,3)
    c0 = prepare_inputs(B,zeros(3,3))
    @test size(c0.Z) == (3,0)
    @test c0.null_rank == 3
    @test space_projector(c0.X0) ≈ B
    c1 = prepare_inputs(B,Matrix(Diagonal([1.,4.,9.])))
    @test size(c1.X0) == (3,0)
    @test c1.Z*c1.Z' ≈ Diagonal([1.,1/4,1/9])
end

@testset "Invalid inputs fail explicitly" begin
    B = Matrix{Float64}(I,2,2)
    @test_throws DimensionMismatch prepare_inputs(B,zeros(3,3))
    @test_throws ArgumentError prepare_inputs(B,[1. 1.;0. 1.])
    @test_throws ArgumentError prepare_inputs(B,Matrix(Diagonal([1.,-1.])))
    @test_throws ArgumentError prepare_inputs([NaN 0.;0. 1.],B)
    @test_throws ArgumentError prepare_inputs(B,B;expected_rank=1)
    @test_throws ArgumentError prepare_inputs(B,B;expected_rank=3)
end

println("Workflow preparation mathematics checks passed.")

isdefined(@__MODULE__,:BSplineBasis) || include(joinpath(@__DIR__,"../../spline_basis/src/BSpline.jl"))
using .BSplineBasis
@testset "Additive assembly and retained prediction transformations" begin
    rng=MersenneTwister(77)
    for m in (1,2,3,10)
        xs=[sort(4rand(rng,160).-2) for _ in 1:m]
        # Randomize row pairing independently for every component.
        xs=[x[randperm(rng,length(x))] for x in xs]
        bases=[build_spline(x,default_knots(x,5+mod(j,3))) for (j,x) in enumerate(xs)]
        p=prepare_additive(bases)
        @test size(p.X_fixed)==(160,m+1)
        @test length(p.Z_blocks)==m
        raw=[b.B_raw for b in bases]
        new=prediction_inputs(bases,p,raw)
        @test new.X_fixed ≈ p.X_fixed atol=1e-12
        for j in 1:m
            @test new.Z_blocks[j] ≈ p.Z_blocks[j] atol=1e-12
            @test length(p.fixed_ranges[j])==1
        end
        different=[evaluate_basis(b.knots,[-.5,.5]) for b in bases]
        @test size(prediction_inputs(bases,p,different).X_fixed)==(2,m+1)
    end
    x=collect(range(-1,1;length=40));b=build_spline(x,default_knots(x,5))
    @test_throws ArgumentError prepare_additive([b,b])
    @test_throws ArgumentError prepare_additive([b];names=["a","b"])
end
