# Stage 2 results

R and Julia independently prepare and fit the same additive model. R uses `mgcv::gamm` / `nlme::lme`; Julia uses `MixedModels.jl`. REML jointly estimates a separate spline variance for each predictor and one common noise variance.

**90/90 comparisons passed across 5 cases.** Spline-component counts: 2.

## Smoothing estimates

Each component has its own $\lambda_j=\sigma^2/r_j^2$. The table compares these smoothing strengths; the variance values and their errors are retained in the detailed estimates.

| Case | Component | R λ<sub>j</sub> | Julia λ<sub>j</sub> | Relative error |
| --- | --- | ---: | ---: | ---: |
| `uniform_k8` | $x_1$ | 0.56176456 | 0.56176455 | 2.99e-09 |
| `uniform_k8` | $x_2$ | 5.4194393 | 5.4194416 | 4.31e-07 |
| `default_k8` | $x_1$ | 0.56176456 | 0.56176455 | 2.96e-09 |
| `default_k8` | $x_2$ | 5.4194393 | 5.4194416 | 4.31e-07 |
| `uneven_k10` | $x_1$ | 4.1679192 | 4.1679195 | 8.24e-08 |
| `uneven_k10` | $x_2$ | 0.83912958 | 0.83912953 | 5.92e-08 |
| `duplicates_k6` | $x_1$ | 0.079999323 | 0.079999328 | 6.55e-08 |
| `duplicates_k6` | $x_2$ | 5.2855236 | 5.2855236 | 7.03e-09 |
| `small_k5` | $x_1$ | 0.0085183673 | 0.0085183663 | 1.20e-07 |
| `small_k5` | $x_2$ | 0.20653478 | 0.2065347 | 3.87e-07 |

Largest scaled difference in training fitted responses: **5.73e-09**; on the shared evaluation grid: **6.00e-09**.

## Joint response and individual components

<a href="spline_fitting/outputs/shared/comparison/fitted_surfaces.png"><img src="spline_fitting/outputs/shared/comparison/fitted_surfaces.png" alt="Six panels: known, R and Julia surfaces with training observations; and three surface difference maps" width="900"></a>

Panels (a–c) overlay the same **200 training observations** as grey points on translucent surfaces, using identical viewing angles and response limits. The points are noisy observed responses. Panels (d–f) compare surfaces on the prediction grid: (d) Julia minus R, with a maximum absolute difference of **1.19e-08**; (e–f) each fit minus the known surface, sharing a **−0.3 to 0.3** colour scale. Red is positive and blue is negative; all difference scales are centred at zero.

### Difference from the known surface

Across 961 grid points, RMSE is **0.0794 for R** and **0.0794 for Julia**. The largest absolute error across both fits is **0.2531**. These errors reflect recovery of the known curve from noisy training data; the R–Julia comparison above checks agreement between implementations.

[Numerical summary](spline_fitting/outputs/shared/comparison/truth_comparison.tsv). RMSE is the square root of the mean squared prediction error on this grid.

<a href="spline_fitting/outputs/shared/comparison/component_curves.png"><img src="spline_fitting/outputs/shared/comparison/component_curves.png" alt="Centered additive spline components" width="800"></a>

Each panel isolates one spline contribution. Known components are centered over the training observations to match the fitted smooths; their removed average levels belong to the overall intercept.

## Supporting checks

[Individual comparisons](spline_fitting/outputs/shared/comparison/checks.tsv) and [detailed estimates](spline_fitting/outputs/shared/comparison/estimates.tsv) record all cases. The acceptance criteria are defined in [VALIDATION.md](VALIDATION.md).

[Julia unit tests](spline_fitting/outputs/shared/unit_tests.txt) independently check the stored designs, covariance, REML objective, coefficients and basis invariance for 1, 2, 3 and 10 components, including unequal block sizes and the zero-variance boundary.

[R fitting status](spline_fitting/outputs/shared/r_reference/fit_summary.tsv); [Julia fitting status](spline_fitting/outputs/shared/julia/fit_summary.tsv); [source fingerprints](spline_fitting/outputs/shared/r_reference/source_checksums.tsv).

Both languages use the same saved observations, per-predictor settings and evaluation grid. Package versions are recorded in each branch's `versions.tsv`; Julia dependencies are pinned by the root `Manifest.toml`.
