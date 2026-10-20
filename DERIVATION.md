# Why the spline and mixed-model representations agree

[Project overview](README.md) · [Validation method](VALIDATION.md)

This document explains basis construction and the fixed/random input preparation inside the combined Stage 1 workflow. A superscript $\top$ denotes matrix transpose. The current spline is cubic, with a penalty on its second derivative, matching mgcv `bs="bs", m=c(3,2)`.

## 1. Describe the curve and its bending

A spline represents the curve as a weighted sum of basis functions:

```math
\begin{aligned}
f(x)&=\sum_{j=1}^{k}c_j b_j(x),\\
(B_{\mathrm{raw}})_{ij}&=b_j(x_i).
\end{aligned}
```

Each column of $B_{\mathrm{raw}}$ is one shape evaluated at the observations. The unknown coefficients $c$ are weights, not values at individual knots. The full cubic knot vector has $k+4$ entries, including exterior knots. Our main example uses $k=8$.

The curvature penalty is:

```math
\begin{aligned}
(S_{\mathrm{raw}})_{ij}&=\int b_i''(x)b_j''(x)\,dx,\\
c^\top S_{\mathrm{raw}}c&=\int [f''(x)]^2\,dx.
\end{aligned}
```

Integration covers the middle knot interval. Straight lines have zero second derivative, so constant and linear curves receive no penalty. This gives two unpenalized directions before centering. Uneven knot spacing is accounted for through the actual interval widths.

Sources: [mgcv cubic B-splines and curvature penalties](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.bs.smooth.spec.html), [basis dimensions and derivatives](https://stat.ethz.ch/R-manual/R-devel/library/splines/html/splineDesign.html), [B-spline recurrence](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html#notes).

## 2. Scale and center the spline

A positive scale $a$ normalizes the penalty using mgcv's convention. A matrix $H$ restricts the spline to have mean zero over the training observations:

```math
\begin{aligned}
a&=\frac{\lVert S_{\mathrm{raw}}\rVert_1}{\lVert B_{\mathrm{raw}}\rVert_\infty^2},\\
B&=B_{\mathrm{raw}}H,\\
S&=H^\top(S_{\mathrm{raw}}/a)H.
\end{aligned}
```

These are induced matrix norms: the maximum absolute column sum and row sum, respectively. The columns of $H$ are orthonormal and satisfy the centering constraint. Centering separates the curve's average level from the model intercept.

With raw weights $c=Hz$, the prepared curve is $Bz$, and its normalized penalty is $z^\top Sz$. For the main example, $B$ has seven columns and $S$ has six positive-penalty directions and one zero-penalty direction.

Source: [mgcv scaling and constraints](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smoothCon.html).

## 3. Prepare fixed and random input columns

The eigenvectors of $S$ identify coefficient directions with separate penalty contributions. Split them into positive-penalty directions $U_+$, with positive eigenvalues in the diagonal matrix $D_+$, and zero-penalty directions $U_0$:

```math
\begin{aligned}
S&=U_+D_+U_+^\top,\\
SU_0&=0.
\end{aligned}
```

Rescale the positive directions by the inverse square roots of their eigenvalues:

```math
\begin{aligned}
C&=[U_+D_+^{-1/2}\;\;U_0],\\
z&=C\begin{bmatrix}r\\\alpha\end{bmatrix}.
\end{aligned}
```

```math
\begin{aligned}
Z&=BU_+D_+^{-1/2},\\
X_0&=BU_0.
\end{aligned}
```

The weights $r$ belong to penalized components; $\alpha$ belongs to unpenalized components. Two identities explain the conversion:

```math
\begin{aligned}
Bz&=Zr+X_0\alpha,\\
z^\top Sz&=r^\top r.
\end{aligned}
```

The curve is unchanged, and the penalty becomes the sum of squared penalized weights. Every coefficient direction is retained. Source: [mgcv conversion and coefficient mapping](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html).

## 4. Link smoothing to random-effect variance

Let $X$ contain ordinary terms such as an intercept. The same response model can be written as:

```math
\begin{aligned}
y&=X\beta+Bz+\varepsilon\\
 &=X_{\mathrm{fixed}}b+Zr+\varepsilon.
\end{aligned}
```

where $X_{\mathrm{fixed}}=[X,\;X_0]$ and $b$ combines $\beta$ and $\alpha$. The mixed model assumes independent Gaussian vectors:

```math
\begin{aligned}
r&\sim N(0,\tau^2I),\\
\varepsilon&\sim N(0,\sigma^2I).
\end{aligned}
```

$I$ is an identity matrix. All penalized coefficients share one variance. For fixed positive variances, their joint density leads to the penalized fitting objective:

```math
\begin{gathered}
\lVert y-X_{\mathrm{fixed}}b-Zr\rVert^2+\lambda r^\top r,\\
\lambda=\frac{\sigma^2}{\tau^2}.
\end{gathered}
```

Using the identities above, this equals the original objective:

```math
\lVert y-X\beta-Bz\rVert^2+\lambda z^\top Sz.
```

Thus variance-based smoothing and the spline penalty describe the same coefficient-fitting problem. Estimating the variances through ML or REML is a later fitting step. These random effects represent spline shapes; separate participant groups are not required.

Sources: [penalized fitting objective](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/magic.html), [smooths represented as mixed models](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gamm.html).

## Training settings and prediction

The inspected mgcv 1.9.4 default uses uniform knots over a predictor range expanded by 0.1% at each end, plus exterior knots. Explicit uneven knots are also supported. The current model uses strictly increasing full knots.

Future predictions will retain the training knots, centering, scaling and conversion. Outside the middle knot interval, the basis follows mgcv's linear extrapolation. Source: [mgcv construction and extrapolation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.bs.smooth.spec.html); the [development source](https://github.com/cran/mgcv/blob/master/R/smooth.r) may change, while the recorded results identify the reference versions used.
