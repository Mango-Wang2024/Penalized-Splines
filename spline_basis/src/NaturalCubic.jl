module NaturalCubic

using LinearAlgebra, Statistics
export default_knots, second_derivative_map, evaluate_basis, build_spline

function checked_knots(knots)
    t = Float64.(collect(knots))
    length(t) >= 3 || throw(ArgumentError("Use at least three knots"))
    all(isfinite,t) || throw(ArgumentError("Knots must be finite"))
    all(>(0),diff(t)) || throw(ArgumentError("Knots must be strictly increasing"))
    return t
end

"""mgcv's cr default: type-7 quantiles of the unique predictor values."""
function default_knots(x,k::Int)
    v = sort(unique(Float64.(collect(x))))
    all(isfinite,v) || throw(ArgumentError("Predictors must be finite"))
    3 <= k <= length(v) || throw(ArgumentError("Require 3 <= k <= unique predictor count"))
    # Explicit interpolation fixes the quantile convention across languages.
    return [begin
        index = 1+(length(v)-1)*p
        lo = floor(Int,index)
        hi = min(lo+1,length(v))
        v[lo]+(index-lo)*(v[hi]-v[lo])
    end for p in range(0,1;length=k)]
end

"""Map values at knots to second derivatives of the natural cubic spline."""
function second_derivative_map(knots)
    t = checked_knots(knots)
    k = length(t)
    h = diff(t)
    A = SymTridiagonal(2 .* (h[1:end-1]+h[2:end]),h[2:end-1])
    D = zeros(k-2,k)
    for j in 1:k-2
        D[j,j] = 6/h[j]
        D[j,j+1] = -6*(1/h[j]+1/h[j+1])
        D[j,j+2] = 6/h[j+1]
    end
    F = zeros(k,k)
    F[2:k-1,:] = A\D
    return F
end

"""Evaluate cardinal spline shapes, or their first/second derivatives.

Natural cubic splines interpolate the knot values, have zero curvature at the
end knots, and extend linearly outside the knot interval.
"""
function evaluate_basis(knots,F::AbstractMatrix,x;derivative::Int=0)
    t = checked_knots(knots)
    k = length(t)
    size(F)==(k,k) || throw(DimensionMismatch("F must be k-by-k"))
    derivative in 0:2 || throw(ArgumentError("Derivative must be 0, 1 or 2"))
    xv = Float64.(collect(x))
    all(isfinite,xv) || throw(ArgumentError("Predictors must be finite"))
    B = zeros(length(xv),k)
    for (row,value) in enumerate(xv)
        if value < t[1] || value > t[end]
            derivative==2 && continue
            left = value < t[1]
            j = left ? 1 : k-1
            h = t[j+1]-t[j]
            slope = zeros(k)
            slope[j] = -1/h
            slope[j+1] = 1/h
            if left
                slope -= h .* F[2,:]/6
            else
                slope += h .* F[k-1,:]/6
            end
            if derivative==1
                B[row,:] = slope
            else
                endpoint = left ? 1 : k
                B[row,endpoint] = 1
                B[row,:] += (value-t[endpoint]) .* slope
            end
            continue
        end
        j = clamp(searchsortedlast(t,value),1,k-1)
        h = t[j+1]-t[j]
        a,b = (t[j+1]-value)/h,(value-t[j])/h
        if derivative==0
            B[row,j] = a
            B[row,j+1] = b
            B[row,:] += (h^2/6) .* ((a^3-a).*F[j,:]+(b^3-b).*F[j+1,:])
        elseif derivative==1
            B[row,j] = -1/h
            B[row,j+1] = 1/h
            B[row,:] += (h/6) .* ((1-3a^2).*F[j,:]+(3b^2-1).*F[j+1,:])
        else
            B[row,:] = a.*F[j,:]+b.*F[j+1,:]
        end
    end
    return B
end

"""Construct raw, scaled and centered B/S matrices independently of R."""
function build_spline(x,knots)
    t = checked_knots(knots)
    xv = Float64.(collect(x))
    length(unique(xv)) >= length(t) || throw(ArgumentError("Need at least k distinct predictor values"))
    F = second_derivative_map(t)
    B_raw = evaluate_basis(t,F,xv)
    h = diff(t)
    # Integral of the squared piecewise-linear second derivative.
    diag_w = vcat(h[1]/3,(h[1:end-1]+h[2:end])/3,h[end]/3)
    W = SymTridiagonal(diag_w,h/6)
    S_raw = Matrix(Symmetric(F'*W*F))
    # smoothCon's matrix norms: one-norm for S, infinity-norm for B.
    penalty_scale = opnorm(S_raw,1)/opnorm(B_raw,Inf)^2
    S_scaled = S_raw/penalty_scale
    constraint = vec(mean(B_raw;dims=1))
    # One-column QR uses a Householder reflection. Use the same sign convention
    # as R's default QR so that the centered entries can be compared directly.
    v = copy(constraint)
    v[1] += copysign(norm(constraint),constraint[1])
    Q = Matrix{Float64}(I,length(t),length(t))-2*(v*v')/dot(v,v)
    H = Q[:,2:end]
    B = B_raw*H
    S = Matrix(Symmetric(H'*S_scaled*H))
    return (;knots=t,F,B_raw,S_raw,S_scaled,penalty_scale,constraint,H,B,S)
end

end
