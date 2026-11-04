module ModelInputs

using LinearAlgebra
export prepare_inputs, prepare_additive, prediction_inputs, relative_error, space_projector

# Relative cutoff for numerical zero eigenvalues and singular values.
const RANK_RTOL = 1e-10

"""Scaled numerical difference; the denominator has a floor of one."""
relative_error(a, b) = norm(a - b) / max(1.0, norm(a), norm(b))

"""
Prepare fixed/random spline columns inside the combined workflow.

Input B is n-by-p, and S is p-by-p. The returned coefficient map C satisfies
    z = C * [tau; alpha]
    B*z = Z*tau + X0*alpha
    z'*S*z = tau'*tau.
"""
function prepare_inputs(B::AbstractMatrix, S::AbstractMatrix)
    p = size(B, 2)
    size(S) == (p, p) || throw(DimensionMismatch("S must be p-by-p for the supplied B"))
    Bf, Sf = Matrix{Float64}(B), Matrix{Float64}(S)
    all(isfinite, Bf) && all(isfinite, Sf) || throw(ArgumentError("B and S must be finite"))
    norm(Sf - Sf') <= RANK_RTOL * max(norm(Sf), eps(Float64)) ||
        throw(ArgumentError("S is not symmetric within tolerance"))
    E = eigen(Symmetric((Sf + Sf') / 2))
    order = sortperm(E.values; rev=true) # R's symmetric eigen() sorts descending.
    values = E.values[order]
    U = E.vectors[:, order]
    tol = RANK_RTOL * maximum(abs, values)
    minimum(values) >= -tol || throw(ArgumentError("S is not positive semidefinite"))
    r = count(>(tol), values)
    scales = vcat(1 ./ sqrt.(values[1:r]), ones(p-r))
    C = U * Diagonal(scales)
    transformed = Bf * C
    return (; C, Z=transformed[:, 1:r], X0=transformed[:, r+1:end],
             penalty_rank=r, null_rank=p-r)
end

"""Prepare any number of centered spline components and assemble their fixed columns."""
function prepare_additive(bases, X::AbstractMatrix)
    isempty(bases) && throw(ArgumentError("Supply at least one spline component"))
    X=Matrix{Float64}(X)
    components=[prepare_inputs(b.B,b.S) for b in bases]
    offset=size(X,2)
    fixed_ranges=UnitRange{Int}[]
    for c in components
        push!(fixed_ranges,(offset+1):(offset+size(c.X0,2)))
        offset+=size(c.X0,2)
    end
    X_fixed=hcat(X,(c.X0 for c in components)...)
    rank(X_fixed)==size(X_fixed,2) || throw(ArgumentError("Combined fixed columns must be linearly independent"))
    return (;components,X_fixed,Z_blocks=[c.Z for c in components],
             fixed_ranges,ordinary_columns=size(X,2))
end

"""Evaluate training transformations at new points; do not recenter the new rows."""
function prediction_inputs(bases,prepared,raw_matrices,X::AbstractMatrix)
    length(raw_matrices)==length(bases)==length(prepared.components) || throw(DimensionMismatch("Component counts disagree"))
    size(X,2)==prepared.ordinary_columns || throw(DimensionMismatch("New ordinary fixed columns must match training"))
    X=Matrix{Float64}(X)
    Z_blocks=Matrix{Float64}[];unpenalized=Matrix{Float64}[]
    for (basis,component,raw) in zip(bases,prepared.components,raw_matrices)
        transformed=raw*basis.H*component.C
        r=component.penalty_rank
        push!(Z_blocks,transformed[:,1:r]);push!(unpenalized,transformed[:,r+1:end])
    end
    return (;X_fixed=hcat(X,unpenalized...),Z_blocks)
end

"""Orthogonal projector onto the column space, including rank-deficient inputs."""
function space_projector(A::AbstractMatrix)
    n, p = size(A)
    p == 0 && return zeros(n,n)
    F = svd(Matrix{Float64}(A); full=false)
    cutoff = RANK_RTOL * maximum(F.S; init=0.0)
    q = count(>(cutoff), F.S)
    Q = F.U[:,1:q]
    return Q * Q'
end

end
