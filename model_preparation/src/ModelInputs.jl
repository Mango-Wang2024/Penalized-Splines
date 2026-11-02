module ModelInputs

using LinearAlgebra
export prepare_inputs, prepare_additive, prediction_inputs, relative_error, space_projector

"""Scaled numerical difference; the denominator has a floor of one."""
relative_error(a, b) = norm(a - b) / max(1.0, norm(a), norm(b))

"""
Prepare fixed/random spline columns inside the combined workflow.

Input B is n-by-p, and S is p-by-p. The returned coefficient map C satisfies
    z = C * [tau; alpha]
    B*z = Z*tau + X0*alpha
    z'*S*z = tau'*tau.

This helper prepares a supplied basis and penalty for the combined workflow. It does NOT construct a
spline basis, estimate lambda, or fit a mixed model.
"""
function prepare_inputs(B::AbstractMatrix, S::AbstractMatrix;
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

"""Prepare any number of centered spline components and assemble their fixed columns."""
function prepare_additive(bases;names=["x"*string(j) for j in eachindex(bases)],ordinary=nothing)
    isempty(bases) && throw(ArgumentError("Supply at least one spline component"))
    length(names)==length(bases) && length(unique(names))==length(names) || throw(ArgumentError("Use one distinct name per component"))
    n=size(first(bases).B,1)
    X=ordinary===nothing ? ones(n,1) : Matrix{Float64}(ordinary)
    size(X,1)==n || throw(DimensionMismatch("Ordinary fixed columns must share observation rows"))
    components=[prepare_inputs(b.B,b.S) for b in bases]
    all(size(c.B,1)==n for c in components) || throw(DimensionMismatch("Spline components must share observation rows"))
    offset=size(X,2)
    fixed_ranges=UnitRange{Int}[]
    for c in components
        push!(fixed_ranges,(offset+1):(offset+size(c.X0,2)))
        offset+=size(c.X0,2)
    end
    X_fixed=hcat(X,(c.X0 for c in components)...)
    rank(X_fixed)==size(X_fixed,2) || throw(ArgumentError("Combined fixed columns must be linearly independent"))
    return (;names=String.(names),components,X_fixed,Z_blocks=[c.Z for c in components],
             fixed_ranges,ordinary_columns=size(X,2))
end

"""Evaluate training transformations at new points; do not recenter the new rows."""
function prediction_inputs(bases,prepared,raw_matrices;ordinary=nothing)
    length(raw_matrices)==length(bases)==length(prepared.components) || throw(DimensionMismatch("Component counts disagree"))
    n=size(first(raw_matrices),1)
    X=ordinary===nothing ? ones(n,prepared.ordinary_columns) : Matrix{Float64}(ordinary)
    ordinary===nothing && prepared.ordinary_columns!=1 && throw(ArgumentError("Supply new ordinary fixed columns"))
    size(X)==(n,prepared.ordinary_columns) || throw(DimensionMismatch("New ordinary fixed matrix has incorrect dimensions"))
    Z_blocks=Matrix{Float64}[];unpenalized=Matrix{Float64}[]
    for (basis,component,raw) in zip(bases,prepared.components,raw_matrices)
        size(raw)==(n,size(basis.B_raw,2)) || throw(DimensionMismatch("New raw basis has incorrect dimensions"))
        transformed=raw*basis.H*component.C
        r=component.penalty_rank
        push!(Z_blocks,transformed[:,1:r]);push!(unpenalized,transformed[:,r+1:end])
    end
    return (;X_fixed=hcat(X,unpenalized...),Z_blocks)
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
