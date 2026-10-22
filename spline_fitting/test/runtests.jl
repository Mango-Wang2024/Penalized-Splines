using Test, LinearAlgebra, MixedModels
include(joinpath(@__DIR__,"..","src","SplineFit.jl"))
using .SplineFit

# Independent Gaussian marginal calculation in the original observation coordinates.
# This evaluates the fitted solution; it does not fix lambda before fitting.
function marginal_check(y,X,Z,theta)
    n,p=size(X)
    V=cholesky(Symmetric(Matrix{Float64}(I,n,n)+theta^2*(Z*Z')))
    XtViX=cholesky(Symmetric(X'*(V\X)))
    beta=XtViX\(X'*(V\y))
    residual=y-X*beta
    rss=dot(residual,V\residual)
    sigma2=rss/(n-p)
    tau=theta^2*Z'*(V\residual)
    objective=logdet(V)+logdet(XtViX)+(n-p)*(1+log(2pi*sigma2))
    return (;beta,tau,sigma2,objective,fitted=X*beta+Z*tau)
end

@testset "Native MixedModels adapter and automatic REML" begin
    x=collect(range(-1,1;length=80))
    X=hcat(ones(80),x)
    Z=hcat(x.^2,x.^3,x.^4)
    y=1 .+ 2x .+ sin.(3x) .+ 0.1cos.(83x)
    prepared=build_model(y,X,Z)
    m=prepared.model
    @test only(m.reterms) isa MixedModels.ReMat{Float64,1}
    @test length(m.θ)==1
    @test length(only(m.reterms).levels)==size(Z,2)
    # Explicit rotation is constructed only in this small independent test.
    T=restore_observations(prepared,Matrix{Float64}(I,80,80))'
    @test T*T' ≈ I atol=1e-12
    @test restore_observations(prepared,T*y) ≈ y atol=1e-12
    @test restore_observations(prepared,T*X) ≈ X atol=1e-12
    packed=Matrix(only(m.reterms))
    @test packed ≈ T*Z*prepared.V atol=1e-12
    @test packed*packed' ≈ T*(Z*Z')*T' atol=1e-12
    @test isapprox(abs(det(T)),1;atol=1e-12)
    fit=fit_reml(y,X,Z)
    independent=marginal_check(y,X,Z,fit.theta)
    @test fit.status in (:SUCCESS,:STOPVAL_REACHED,:FTOL_REACHED,:XTOL_REACHED)
    @test fit.lambda ≈ fit.sigma2/fit.r2
    @test fit.lambda ≈ inv(fit.theta^2)
    @test fit.beta ≈ independent.beta rtol=1e-9
    @test fit.tau ≈ independent.tau rtol=1e-9
    @test fit.sigma2 ≈ independent.sigma2 rtol=1e-9
    @test fit.objective ≈ independent.objective rtol=1e-9
    @test fit.fitted ≈ independent.fitted rtol=1e-9
    # Rotating equally penalized directions must not change the fitted model.
    Q=Matrix(qr([1.0 2 0; 0 1 3; 2 0 1]).Q)
    rotated=fit_reml(y,X,Z*Q)
    @test fit.lambda ≈ rotated.lambda rtol=1e-5
    @test fit.fitted ≈ rotated.fitted rtol=1e-7
    @test fit.sigma2 ≈ rotated.sigma2 rtol=1e-6
    @test fit.tau ≈ Q*rotated.tau rtol=1e-6
    # An exactly linear signal with noise orthogonal to all spline directions
    # has its optimum at zero spline variance (infinite smoothing).
    W=hcat(X,Z)
    noise=cos.(31x)-W*(W\cos.(31x))
    boundary=fit_reml(X*[1.0,2.0]+0.1noise,X,Z)
    @test iszero(boundary.theta)
    @test isinf(boundary.lambda)
    @test iszero(boundary.r2)
    @test boundary.fitted ≈ X*[1.0,2.0] atol=1e-10
    @test_throws DimensionMismatch build_model(y[2:end],X,Z)
    @test_throws DimensionMismatch build_model(y,X,Z[2:end,:])
    @test_throws ArgumentError build_model(fill(NaN,80),X,Z)
    @test_throws ArgumentError build_model(y,hcat(x,x),Z)
    @test_throws ArgumentError build_model(y,X,hcat(x.^2,x.^2))
    @test_throws ArgumentError build_model(y[1:4],X[1:4,:],Z[1:4,:])
end
