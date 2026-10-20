# Stage 2 results

R and Julia independently construct the same spline model and estimate smoothing automatically with **REML**. R uses `mgcv::gamm` with `nlme::lme`; Julia uses prepared `FeTerm` and `ReMat` inputs with `MixedModels.LinearMixedModel`.

**50/50 reference comparisons passed across five cases.** No smoothing parameter is fixed before fitting.

## Estimated smoothing and fitted curves

The relative error compares each pair of variance or smoothing estimates. The fitted-curve error compares all 200 fitted responses in their original observation order, using the scaled norm defined in VALIDATION.md.

| Case | R lambda | Julia lambda | Lambda relative error | Fitted-curve scaled error |
| --- | ---: | ---: | ---: | ---: |
| `uniform_k8` | 0.60963055 | 0.609630639 | 1.45e-07 | 1.50e-09 |
| `default_k8` | 0.609630551 | 0.609630639 | 1.44e-07 | 1.50e-09 |
| `uneven_k10` | 4.42327209 | 4.42327175 | 7.49e-08 | 9.50e-10 |
| `duplicates_k6` | 0.0534448003 | 0.0534448028 | 4.61e-08 | 2.28e-10 |
| `small_k5` | 0.00962105738 | 0.00962105456 | 2.94e-07 | 2.06e-09 |

## Estimated variances

All penalized spline directions share one variance, tau squared. The remaining observation noise has variance sigma squared. Their ratio is the smoothing parameter: lambda = sigma squared / tau squared.

| Case | R noise variance | Julia noise variance | R spline variance | Julia spline variance |
| --- | ---: | ---: | ---: | ---: |
| `uniform_k8` | 0.140939245 | 0.140939245 | 0.231187963 | 0.23118793 |
| `default_k8` | 0.140939245 | 0.140939245 | 0.231187963 | 0.23118793 |
| `uneven_k10` | 0.141393509 | 0.141393508 | 0.0319658176 | 0.0319658199 |
| `duplicates_k6` | 0.104718631 | 0.104718631 | 1.95937922 | 1.95937913 |
| `small_k5` | 0.153465772 | 0.153465772 | 15.9510297 | 15.9510343 |

The largest relative variance error is **2.89e-07**. The acceptance thresholds are 1e-5 for smoothing and variances, and 1e-6 for fitted curves.

![R and Julia fitted curves](spline_fitting/outputs/shared/comparison/fitted_curves.png)

## Supporting checks

Independent Julia tests check that the input adapter preserves the covariance model, REML objective and fitted coefficients. They also check invariance to rotated spline directions, the zero-variance boundary and invalid inputs.

[Detailed comparisons](spline_fitting/outputs/shared/comparison/checks.tsv); [all estimates](spline_fitting/outputs/shared/comparison/estimates.tsv); [Julia test log](spline_fitting/outputs/shared/unit_tests.txt). [R fitting status](spline_fitting/outputs/shared/r_reference/fit_summary.tsv) and [Julia fitting status](spline_fitting/outputs/shared/julia/fit_summary.tsv) record the individual fits.

## Reproducibility

Versions: R version 4.6.1 (2026-06-24); mgcv 1.9.4; nlme 3.1.169; Julia 1.13.1; MixedModels 5.9.0; StatsModels 0.7.10.

The same saved observations and five knot settings are reused from Stage 1. [Input fingerprints](spline_fitting/outputs/shared/r_reference/source_checksums.tsv) confirm the R source files; both languages' saved predictor, response and knot values are compared directly. The Julia dependency versions are pinned in `spline_fitting/Manifest.toml`.
