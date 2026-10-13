# Construction equations

The coefficients represent the curve's values at ordered knots `t`. The gap between adjacent knots is `h[i] = t[i+1] - t[i]`.

## Raw basis

A natural cubic spline joins knot values smoothly and has zero second derivative at both end knots.

For interior knots, the second derivatives `m` satisfy:

```text
h[i-1] m[i-1] + 2(h[i-1]+h[i]) m[i] + h[i] m[i+1]
    = 6((z[i+1]-z[i])/h[i] - (z[i]-z[i-1])/h[i-1])
```

Solving this system gives `m = F z`. Evaluating the resulting cubic segments gives `f(x) = B_raw z`. Outside the end knots, the curve extends linearly.

## Raw penalty

The penalty adds up squared curvature over the knot interval:

```text
z' S_raw z = integral of (f''(x))²
S_raw = F' W F
```

`W` integrates the piecewise-linear second derivative. Each interval contributes `h/3` to its two diagonal entries and `h/6` to its off-diagonal entries.

## Scaling and centering

To reproduce `smoothCon` in the tested version:

```text
scale = one_norm(S_raw) / infinity_norm(B_raw)²
S_scaled = S_raw / scale
```

The one-norm is the largest absolute column sum. The infinity-norm is the largest absolute row sum.

A mean-zero constraint avoids duplicating the model intercept. An independently computed Householder transformation `H` uses R's one-column QR sign convention:

```text
B = B_raw H
S = H' S_scaled H
```

These are the prepared matrices passed to the unchanged spline-to-mixed-model conversion.

Sources: [natural cubic specification](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.cr.smooth.spec.html), [constraint and scaling interface](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smoothCon.html), and the locally inspected mgcv 1.9.4 functions.
