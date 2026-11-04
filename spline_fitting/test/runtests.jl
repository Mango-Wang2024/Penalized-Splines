using Test, LinearAlgebra, Random, MixedModels
include(joinpath(@__DIR__,"../src/SplineFit.jl"));using .SplineFit

# Independent original-coordinate Gaussian calculations, without adapter rotations.
function marginal_check(y,X,Zs,theta)
    n,p=size(X)
    V=cholesky(Symmetric(Matrix{Float64}(I,n,n)+sum(theta[j]^2*(Zs[j]*Zs[j]') for j in eachindex(Zs))))
    G=cholesky(Symmetric(X'*(V\X)))
    b=G\(X'*(V\y));res=y-X*b;sigma2=dot(res,V\res)/(n-p)
    tau=[theta[j]^2*Zs[j]'*(V\res) for j in eachindex(Zs)]
    objective=logdet(V)+logdet(G)+(n-p)*(1+log(2pi*sigma2))
    return (;b,tau,sigma2,objective,fitted=X*b+sum(Z*t for (Z,t) in zip(Zs,tau)))
end

@testset "Variable component counts and shared variances" begin
    for count in (1,2,3,10)
        rng=MersenneTwister(100+count);n=160
        X=hcat(ones(n),randn(rng,n))
        Zs=[randn(rng,n,1+mod(j,3)) for j in 1:count]
        y=X*[.6,-.3]+sum(Z*randn(rng,size(Z,2))*.4 for Z in Zs)+.2randn(rng,n)
        fit=fit_reml(y,X,Zs)
        m=fit.prepared.model
        @test length(m.θ)==count
        @test length(fit.lambda)==count
        @test length(fit.tau)==count
        @test length(m.reterms)==count
        @test all(t.native isa MixedModels.ReMat for t in m.reterms)
        @test all(all(==(t.λ[1,1]),diag(t.native.λ)) for t in m.reterms)
        independent=marginal_check(y,X,Zs,fit.theta)
        @test fit.b ≈ independent.b rtol=1e-9
        @test fit.sigma2 ≈ independent.sigma2 rtol=1e-9
        @test fit.objective ≈ independent.objective rtol=1e-9
        @test fit.fitted ≈ independent.fitted rtol=1e-9
        for j in 1:count
            @test fit.tau[j] ≈ independent.tau[j] rtol=1e-9
            @test fit.lambda[j] ≈ fit.sigma2/fit.r2[j]
        end
        # Test arbitrary trial variances, not only the optimized point.
        trial=collect(range(.2,1.4;length=count==1 ? 2 : count))[1:count]
        model_theta=trial[invperm(fit.prepared.model_order)]
        obj=MixedModels.objective!(m,model_theta)
        @test obj ≈ marginal_check(y,X,Zs,trial).objective rtol=1e-9
        if count==3
            prepared=fit.prepared
            T=restore_observations(prepared,Matrix{Float64}(I,n,n))'
            @test T*T' ≈ I atol=1e-12
            @test restore_observations(prepared,T*y) ≈ y atol=1e-12
            @test restore_observations(prepared,T*X) ≈ X atol=1e-12
            for j in 1:count
                @test prepared.Zr_blocks[j] ≈ T*Zs[j] atol=1e-12
            end
            combined=hcat(prepared.Zr_blocks...);q=size(combined,2)
            @test istriu(combined[1:q,:])
        end
    end
end

@testset "Basis invariance, boundary and invalid inputs" begin
    x=collect(range(-1,1;length=80));X=hcat(ones(80),x)
    Z=hcat(x.^2,x.^3,x.^4);y=1 .+2x.+sin.(3x).+.1cos.(83x)
    fit=fit_reml(y,X,[Z])
    Q=Matrix(qr([1.0 2 0;0 1 3;2 0 1]).Q)
    rotated=fit_reml(y,X,[Z*Q])
    @test fit.lambda ≈ rotated.lambda rtol=1e-5
    @test fit.fitted ≈ rotated.fitted rtol=1e-7
    @test fit.tau[1] ≈ Q*rotated.tau[1] rtol=1e-6
    W=hcat(X,Z);noise=cos.(31x)-W*(W\cos.(31x))
    boundary=fit_reml(X*[1.,2.]+.1noise,X,[Z])
    @test iszero(boundary.theta[1])
    @test isinf(boundary.lambda[1])
    @test boundary.fitted ≈ X*[1.,2.] atol=1e-10
    @test_throws DimensionMismatch build_model(y[2:end],X,[Z])
    @test_throws DimensionMismatch build_model(y,X,[Z[2:end,:]])
    @test_throws ArgumentError build_model(fill(NaN,80),X,[Z])
    @test_throws ArgumentError build_model(y,hcat(x,x),[Z])
    @test_throws ArgumentError build_model(y,X,[hcat(x.^2,x.^2)])
    @test_throws ArgumentError build_model(y,X,Matrix{Float64}[])
    @test_throws DimensionMismatch predict_response(fit,X,[])
end
