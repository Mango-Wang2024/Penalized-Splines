module SplineFit

using LinearAlgebra, MixedModels, StatsModels
export build_model, fit_reml, restore_observations

"""
Package a dense spline design as a native, scalar ReMat with one shared variance.

For Z = Q*[R;0] and R = U*Diagonal(d)*V', rotate observations by
    T = blockdiag(U',I)*Q'
and random coefficients by V'. Then T*Z*V = [Diagonal(d);0].
Both rotations are orthogonal, so the Gaussian model and REML criterion are
unchanged. This scalar design uses genuine ReMat storage, not a replacement
of adjA that would leave its optimized cross-product methods inconsistent.
Only thin QR factors are stored; an n-by-n observation matrix is never built.

Sources: MixedModels 5.9.0 src/linearmixedmodel.jl, src/remat.jl;
Embrace Uncertainty, Chapter 1, Section 1.4.
"""
function build_model(y::AbstractVector, X::AbstractMatrix, Z::AbstractMatrix)
    yf, Xf, Zf = Vector{Float64}(y), Matrix{Float64}(X), Matrix{Float64}(Z)
    n, p = size(Xf)
    size(Zf,1)==n && length(yf)==n || throw(DimensionMismatch("y, X and Z must share observation rows"))
    q = size(Zf,2)
    p>0 && q>0 && n>p+q || throw(ArgumentError("Require fixed columns, penalized columns and n > p+q"))
    all(isfinite,yf) && all(isfinite,Xf) && all(isfinite,Zf) || throw(ArgumentError("Model inputs must be finite"))
    rank(Xf)==p || throw(ArgumentError("Fixed-effect columns must be linearly independent"))

    rotated = rotate_inputs(yf,Xf,Zf)
    terms = build_model_terms(rotated.yr,rotated.Xr,rotated.singular_values)
    model = MixedModels.LinearMixedModel(
        rotated.yr,terms.feterm,[terms.reterm],terms.form,
        Float64[],nothing,false,1000,false)
    length(model.θ)==1 || error("Expected exactly one shared random-effect variance parameter")

    # Keep the singular values (without row padding) for restoring observations.
    return (;model,qrZ=rotated.qrZ,U=rotated.U,V=rotated.V,
             loading=rotated.singular_values,y=yf,X=Xf,Z=Zf)
end

"""Find QR/SVD coordinates and transform the response and fixed effects together."""
function rotate_inputs(yf, Xf, Zf)
    q = size(Zf,2)
    qrZ = qr(Zf)
    parts = svd(Matrix(qrZ.R))
    minimum(parts.S)>1e-10*maximum(parts.S) || throw(ArgumentError("Penalized spline columns must be linearly independent"))
    rotated = qrZ.Q' * hcat(yf,Xf)
    rotated[1:q,:] = parts.U' * rotated[1:q,:]
    yr, Xr = rotated[:,1], rotated[:,2:end]
    return (;yr,Xr,qrZ,U=parts.U,V=parts.V,singular_values=parts.S)
end

"""Package transformed inputs as the formula, FeTerm and scalar ReMat."""
function build_model_terms(yr, Xr, singular_values)
    n, p = size(Xr)
    q = length(singular_values)
    row_loadings = vcat(singular_values,zeros(n-q))
    labels = ["spline_"*lpad(string(i),length(string(q)),'0') for i in 1:q]
    group = vcat(labels,fill(first(labels),n-q))
    names = Symbol.("fixed_" .* string.(1:p))
    tbl = merge((response=yr,loading=row_loadings,spline=group),NamedTuple{Tuple(names)}(Tuple(eachcol(Xr))))
    fixed = foldl(+,term.(names);init=term(0))
    random = (@formula(response ~ (0 + loading | spline))).rhs
    formula = FormulaTerm(term(:response),fixed+random)
    form = MixedModels.schematize(formula,tbl,Dict{Symbol,Any}())
    _, design = StatsModels.modelcols(form,tbl)
    reterm = only(t for t in design if t isa MixedModels.ReMat)
    feterm = MixedModels.FeTerm(Xr,String.(names))
    return (;feterm,reterm,form)
end

"""Return rotated observation values to the original observation order and units."""
function restore_observations(prepared, values::AbstractVecOrMat)
    out = copy(values)
    q = length(prepared.loading)
    if out isa AbstractVector
        out[1:q] = prepared.U * out[1:q]
    else
        out[1:q,:] = prepared.U * out[1:q,:]
    end
    return prepared.qrZ.Q * out
end

"""Estimate the shared variance ratio with REML and default optimizer settings."""
function fit_reml(y::AbstractVector, X::AbstractMatrix, Z::AbstractMatrix)
    prepared = build_model(y,X,Z)
    m = prepared.model
    MixedModels.fit!(m;REML=true,progress=false)
    status = m.optsum.returnvalue
    status in (:SUCCESS,:STOPVAL_REACHED,:FTOL_REACHED,:XTOL_REACHED) ||
        error("MixedModels did not converge: $status")
    theta = abs(only(m.θ)) # r / sigma; not the spline smoothing parameter.
    sigma2 = MixedModels.varest(m)
    r2 = sigma2 * theta^2
    r = sqrt(r2) # Shared standard deviation; tau below is the coefficient vector.
    lambda = iszero(theta) ? Inf : inv(theta^2)
    beta = MixedModels.coef(m)
    tau = prepared.V * vec(only(MixedModels.ranef(m)))
    fitted = restore_observations(prepared,MixedModels.fitted(m))
    reconstructed = prepared.X*beta+prepared.Z*tau
    isapprox(fitted,reconstructed;rtol=1e-10,atol=1e-10) || error("Fitted-value reconstruction failed")
    isfinite(sigma2) && sigma2>0 && all(isfinite,fitted) || error("Invalid fitted model")
    return (;prepared,lambda,sigma2,r2,r,theta,beta,tau,fitted,status,
            objective=MixedModels.objective(m),evaluations=m.optsum.feval)
end

end
