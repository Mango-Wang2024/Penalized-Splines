module SplineConversion

using LinearAlgebra
export convert_spline, relative_error, space_projector

"""Scaled numerical difference; the denominator has a floor of one."""
relative_error(a, b) = norm(a - b) / max(1.0, norm(a), norm(b))

"""
Convert one positive-semidefinite quadratic spline penalty to fixed/random form.

Input B is n-by-p, and S is p-by-p. The returned coefficient map C satisfies
    z = C * [u; alpha]
    B*z = Z*u + X0*alpha
    z'*S*z = u'*u.

This function converts a supplied basis and penalty. It does NOT construct a
spline basis, estimate lambda, or fit a mixed model.
"""
function convert_spline(B::AbstractMatrix, S::AbstractMatrix;
                        expected_rank::Union{Nothing,Int}=nothing,
                        rtol::Real=1e-10, atol::Real=0.0)
    n, p = size(B)
    n > 0 && p > 0 || throw(ArgumentError("B must have positive dimensions"))
    size(S) == (p, p) || throw(DimensionMismatch("S must be p-by-p for the supplied B"))
    rtol >= 0 && atol >= 0 || throw(ArgumentError("Tolerances must be nonnegative"))
    Bf, Sf = Matrix{Float64}(B), Matrix{Float64}(S)
    all(isfinite, Bf) && all(isfinite, Sf) || throw(ArgumentError("B and S must be finite"))
    norm(Sf - Sf') <= atol + rtol * max(norm(Sf), eps(Float64)) ||
        throw(ArgumentError("S is not symmetric within tolerance"))
    E = eigen(Symmetric((Sf + Sf') / 2))
    order = sortperm(E.values; rev=true) # R's symmetric eigen() sorts descending.
    values = E.values[order]
    U = E.vectors[:, order]
    tol = max(Float64(atol), Float64(rtol) * maximum(abs, values))
    minimum(values) >= -tol || throw(ArgumentError("S is not positive semidefinite"))
    r = isnothing(expected_rank) ? count(>(tol), values) : expected_rank
    0 <= r <= p || throw(ArgumentError("Expected rank must be between zero and p"))
    if !isnothing(expected_rank)
        all(>(tol), values[1:r]) || throw(ArgumentError("Expected positive eigenvalues exceed the rank tolerance"))
        all(x -> abs(x) <= tol, values[r+1:end]) ||
            throw(ArgumentError("Expected null space contains nonzero penalty eigenvalues"))
    end
    scales = vcat(1 ./ sqrt.(values[1:r]), ones(p-r))
    C = U * Diagonal(scales)
    transformed = Bf * C
    return (; B=Bf, S=Sf, C, Z=transformed[:, 1:r],
             X0=transformed[:, r+1:end], penalty_rank=r,
             null_rank=p-r, eigenvalues=values, scales, rank_tolerance=tol)
end

"""Orthogonal projector onto the column space, including rank-deficient inputs."""
function space_projector(A::AbstractMatrix; rtol::Real=1e-10)
    n, p = size(A)
    p == 0 && return zeros(n,n)
    F = svd(Matrix{Float64}(A); full=false)
    cutoff = rtol * maximum(F.S; init=0.0)
    q = count(>(cutoff), F.S)
    Q = F.U[:,1:q]
    return Q * Q'
end

end
