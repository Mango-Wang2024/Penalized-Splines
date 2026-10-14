module BSplineBasis

using LinearAlgebra, Statistics
export default_knots, extend_knots, evaluate_basis, build_spline

const DEGREE = 3

function checked_knots(knots)
    t = Float64.(collect(knots))
    length(t) >= 8 || throw(ArgumentError("A cubic B-spline needs at least eight full knot entries"))
    all(isfinite,t) || throw(ArgumentError("Knots must be finite"))
    all(>(0),diff(t)) || throw(ArgumentError("Full knots must be strictly increasing"))
    return t
end

"""Add three exterior knots on each side of the supplied interval breakpoints."""
function extend_knots(breaks)
    b = Float64.(collect(breaks))
    length(b)>=2 && all(isfinite,b) && all(>(0),diff(b)) ||
        throw(ArgumentError("Use at least two strictly increasing finite breakpoints"))
    left,right = b[2]-b[1],b[end]-b[end-1]
    return vcat([b[1]-j*left for j in 3:-1:1],b,[b[end]+j*right for j in 1:3])
end

"""mgcv bs default: uniform full knots, with a 0.1% data-range expansion."""
function default_knots(x,k::Int)
    xv = Float64.(collect(x))
    !isempty(xv) && all(isfinite,xv) || throw(ArgumentError("Predictors must be nonempty and finite"))
    k>=4 || throw(ArgumentError("Use at least four cubic B-spline basis functions"))
    xl,xu = extrema(xv)
    xr = xu-xl
    xr>0 || throw(ArgumentError("Predictor range must be positive"))
    xl -= 0.001*xr
    xu += 0.001*xr
    step = (xu-xl)/(k-DEGREE)
    return collect(range(xl-DEGREE*step,xu+DEGREE*step;length=k+DEGREE+1))
end

# Cox-de Boor recurrence, including first and second derivative recurrences.
function basis_row(t,value,derivative)
    values = Float64[t[i]<=value<t[i+1] for i in 1:length(t)-1]
    first = zeros(length(values))
    second = zeros(length(values))
    for p in 1:DEGREE
        count = length(values)-1
        next_values,next_first,next_second = zeros(count),zeros(count),zeros(count)
        for i in 1:count
            a = t[i+p]-t[i]
            b = t[i+p+1]-t[i+1]
            next_values[i] = (value-t[i])/a*values[i]+(t[i+p+1]-value)/b*values[i+1]
            next_first[i] = p/a*values[i]-p/b*values[i+1]
            next_second[i] = p/a*first[i]-p/b*first[i+1]
        end
        values,first,second = next_values,next_first,next_second
    end
    return derivative==0 ? values : derivative==1 ? first : second
end

"""Evaluate cubic B-spline shapes or derivatives, using mgcv's linear prediction tails."""
function evaluate_basis(knots,x;derivative::Int=0)
    t = checked_knots(knots)
    derivative in 0:2 || throw(ArgumentError("Derivative must be 0, 1 or 2"))
    xv = Float64.(collect(x))
    all(isfinite,xv) || throw(ArgumentError("Predictors must be finite"))
    lo,hi = t[4],t[end-3]
    B = zeros(length(xv),length(t)-4)
    for (i,value) in enumerate(xv)
        if lo<=value<=hi
            B[i,:] = basis_row(t,value,derivative)
        elseif derivative<2
            endpoint = value<lo ? lo : hi
            slope = basis_row(t,endpoint,1)
            B[i,:] = derivative==1 ? slope : basis_row(t,endpoint,0)+(value-endpoint)*slope
        end
    end
    return B
end

"""Construct cubic B-splines with the integrated squared-second-derivative penalty."""
function build_spline(x,knots)
    t = checked_knots(knots)
    xv = Float64.(collect(x))
    k = length(t)-4
    !isempty(xv) && all(isfinite,xv) || throw(ArgumentError("Predictors must be nonempty and finite"))
    length(unique(xv))>=k || throw(ArgumentError("Need at least k distinct predictor values"))
    all(v->t[4]<=v<=t[end-3],xv) || throw(ArgumentError("Training predictors must lie inside the middle knot interval"))
    B_raw = evaluate_basis(t,xv)
    D2 = evaluate_basis(t,xv;derivative=2)
    S_raw = zeros(k,k)
    # Cubic second derivatives are linear on each interval; two-point Gauss
    # quadrature integrates their products exactly, including uneven spacing.
    for j in 4:length(t)-4
        h = t[j+1]-t[j]
        mid = (t[j]+t[j+1])/2
        points = [mid-h/(2sqrt(3)),mid+h/(2sqrt(3))]
        D = evaluate_basis(t,points;derivative=2)
        S_raw += h/2*(D'*D)
    end
    S_raw = Matrix(Symmetric(S_raw))
    penalty_scale = opnorm(S_raw,1)/opnorm(B_raw,Inf)^2
    S_scaled = S_raw/penalty_scale
    constraint = vec(mean(B_raw;dims=1))
    # Preserve the existing mean-zero constraint and R QR sign convention.
    v = copy(constraint)
    v[1] += copysign(norm(constraint),constraint[1])
    Q = Matrix{Float64}(I,k,k)-2*(v*v')/dot(v,v)
    H = Q[:,2:end]
    B = B_raw*H
    S = Matrix(Symmetric(H'*S_scaled*H))
    return (;knots=t,D2,B_raw,S_raw,S_scaled,penalty_scale,constraint,H,B,S)
end

end
