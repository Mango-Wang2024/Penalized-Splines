# Matrices used for spline fitting

Julia's `uniform_k8` example has 200 observations and eight raw basis functions. The panels show model inputs before REML fitting, in construction order.

<a href="fitting_matrices.png"><img src="fitting_matrices.png" alt="Four matrix heatmaps in order: raw spline basis B_raw, centered basis B, penalized design Zf, and rotated penalized design Zr" width="800"></a>

- **(a) $B_{\mathrm{raw}}$ (200 × 8):** raw B-spline values. Each column is one basis function evaluated at the observations.
- **(b) $B$ (200 × 7):** the centered basis, $B=B_{\mathrm{raw}}H$. The centering transformation $H$ imposes a mean-zero constraint, separating the spline from the intercept.
- **(c) $Z_f$ (200 × 6):** penalized directions of $B$, rescaled using the penalty matrix $S$. The remaining unpenalized direction forms $X_0$.
- **(d) $Z_r$ (200 × 6):** the `ReMat` design after QR/SVD observation rotation $T$ and coefficient rotation $V$: $Z_r=T Z_f V$. Its six diagonal entries are singular values of $Z_f$. Only rows 1–10 are shown; rows 7–200 are zero. All observations remain in the model.

Colour indicates each matrix entry: red is positive, blue is negative, and white is zero or near zero. Each panel has its own linear colour scale.

**Sources:** [shared data and knots](../data/inputs/uniform_k8/), [basis construction and centering](../spline_basis/src/BSpline.jl#L69-L92), [penalized/unpenalized preparation](../model_preparation/src/ModelInputs.jl#L20-L49), and [input rotation and packaging](../spline_fitting/src/SplineFit.jl#L42-L69).
