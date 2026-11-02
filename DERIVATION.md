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

The implementation uses BSplineKit for basis values, derivatives and the integral defining $S_{\mathrm{raw}}$. Its cubic basis has order 4 (degree 3); the number of basis functions remains $k$. Supplying our full knot vector without automatic augmentation preserves the mgcv basis and integration interval.

Sources: [mgcv cubic B-splines and curvature penalties](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.bs.smooth.spec.html), [basis dimensions and derivatives](https://stat.ethz.ch/R-manual/R-devel/library/splines/html/splineDesign.html), [BSplineKit bases and full knots](https://jipolanco.github.io/BSplineKit.jl/stable/bsplines/), and [BSplineKit derivative integrals](https://jipolanco.github.io/BSplineKit.jl/stable/galerkin/).

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
z&=C\begin{bmatrix}\tau\\\alpha\end{bmatrix}.
\end{aligned}
```

```math
\begin{aligned}
Z&=BU_+D_+^{-1/2},\\
X_0&=BU_0.
\end{aligned}
```

The weights $\tau$ belong to penalized components; $\alpha$ belongs to unpenalized components. Two identities explain the conversion:

```math
\begin{aligned}
Bz&=Z\tau+X_0\alpha,\\
z^\top Sz&=\tau^\top \tau.
\end{aligned}
```

The curve is unchanged, and the penalty becomes the sum of squared penalized weights. Every coefficient direction is retained. Source: [mgcv conversion and coefficient mapping](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html).

## 4. Assemble any number of additive splines

Steps 1–3 apply separately to each predictor. For $m$ predictors, join the unpenalized parts to the ordinary fixed matrix, while retaining separate penalized blocks:

```math
X_{\mathrm{fixed}}=[X,X_{0,1},\ldots,X_{0,m}],\qquad
b=[\beta^\top,\alpha_1^\top,\ldots,\alpha_m^\top]^\top.
```

```math
y=X_{\mathrm{fixed}}b+\sum_{j=1}^{m}Z_j\tau_j+\varepsilon,
\qquad \tau_j\sim N(0,r_j^2I),\qquad
\varepsilon\sim N(0,\sigma^2I).
```

The random coefficient vectors and noise are mutually independent. Here $\tau_j$ is a vector of spline coefficients, $r_j$ is its shared standard deviation, and $\sigma$ is the noise standard deviation. For fixed positive variances, the joint density gives the coefficient-fitting objective:

```math
\left\|y-X_{\mathrm{fixed}}b-\sum_j Z_j\tau_j\right\|^2
+\sum_j\lambda_j\tau_j^\top\tau_j,
\qquad\lambda_j=\frac{\sigma^2}{r_j^2}.
```

Each spline receives its own smoothing strength. REML estimates the variances jointly; the ratios give the smoothing parameters. The implied response covariance is:

```math
\mathrm{Var}(y)=\sigma^2I+\sum_j r_j^2 Z_jZ_j^\top.
```

This identity also explains why Stage 1 compares each $Z_jZ_j^\top$: equivalent coefficient coordinates must represent the same contribution to covariance. Spline random effects represent curve shapes; participant groups are not required.

Sources: [penalized fitting](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/magic.html), [smooths as mixed models](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gamm.html).

## 5. Adapt the designs to fitting containers

The [Julia adapter](spline_fitting/src/SplineFit.jl) concatenates the penalized matrices only to find one common observation rotation. Write the combined design as $Z_f=[Z_1,\ldots,Z_m]=Q[R;0]$, and factor its small triangular matrix as $R=U D V^\top$. Define:

```math
T=\mathrm{blockdiag}(U^\top,I)Q^\top,
\qquad y_r=Ty,\qquad X_r=TX_{\mathrm{fixed}}.
```

The code applies $Q^\top$ implicitly, without forming a full observation matrix. Its line `Q' * hcat(y, X_fixed)` places the columns side by side and transforms each independently:

```math
T[y\;X_{\mathrm{fixed}}]=[Ty\;TX_{\mathrm{fixed}}].
```

No response column is added to a fixed-effect column. Next, each component receives its own orthogonal coefficient rotation $V_j$, obtained from an SVD of its columns of $R$:

```math
Z_{r,j}=TZ_jV_j,\qquad \tau_{r,j}=V_j^\top\tau_j.
```

A single spline can become diagonal under these transformations. With multiple splines, the upper part of $Z_{r,j}$ is generally dense; rows below the total penalized rank are zero. The row transformation is common to all model terms, and coefficient rotations stay **within** their component. Rotating coefficient directions between components could mix different variances and would change the model.

Orthogonality preserves squared residual lengths and the covariance $r_j^2I$ within each component. Each block therefore retains one variance parameter. Its native `MixedModels.ReMat` stores a dense multi-column term with one synthetic grouping level. [SharedVariance.jl](spline_fitting/src/SharedVariance.jl) ties that term's diagonal covariance entries to one parameter $\theta_j=r_j/\sigma$; off-diagonal entries remain zero. This wrapper uses the pinned MixedModels 5.9.0 low-level API.

The transformed model has the same covariance structure and REML objective. After fitting, coefficients return through $\tau_j=V_j\tau_{r,j}$, and `restore_observations` applies $T^\top$ to fitted responses. Thus each returned fitted value again corresponds to one original data row.

R uses its own internal QR operations in [nlme](https://github.com/cran/nlme/blob/3.1-169/src/matrix.c#L195-L223). Agreement is required for the model and fitted results, not for identical internal rotations.

Sources: pinned [MixedModels constructor](https://github.com/JuliaStats/MixedModels.jl/blob/v5.9.0/src/linearmixedmodel.jl), [random-term storage](https://github.com/JuliaStats/MixedModels.jl/blob/v5.9.0/src/remat.jl), and [estimation](https://juliastats.org/MixedModels.jl/stable/optimization/).

## 6. Evaluate the fitted model at new predictor values

Retain each component's training knots, $H_j$ and $C_j$. Evaluate its raw basis at the new predictor values and apply those same transformations to obtain $X_{0,j,\mathrm{new}}$ and $Z_{j,\mathrm{new}}$. Do not recompute centering from the new rows.

```math
\hat y_{\mathrm{new}}=X_{\mathrm{fixed,new}}b+\sum_j Z_{j,\mathrm{new}}\tau_j.
```

The current two-predictor grid comparison uses this formula. The model has additive components; interactions would require an additional model specification.
