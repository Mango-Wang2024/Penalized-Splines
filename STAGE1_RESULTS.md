# Stage 1 results

Saved comparisons for cubic B-splines with one curvature penalty: mgcv bs="bs", m=c(3,2). All comparisons use shared Julia-generated inputs.

## Overall results

| Check | Passed | Largest scaled error |
| --- | --- | --- |
| Basis construction only | 65/65 | 1.39e-15 |
| Conversion only | 9/9 | 2.08e-14 |
| Combined workflow | 35/35 | 1.73e-14 |

Mathematical tests completed before these comparisons. R input-capture checks passed: 7/7. Comparison threshold: 1e-8.

## Basis construction only: results by case

Five test cases, with 13 comparisons per case.

| Case | Raw basis size k | Passed | Largest scaled error |
| --- | --- | --- | --- |
| uniform_k8 | 8 | 13/13 | 3.33e-16 |
| default_k8 | 8 | 13/13 | 1.15e-15 |
| uneven_k10 | 10 | 13/13 | 9.78e-16 |
| duplicates_k6 | 6 | 13/13 | 1.39e-15 |
| small_k5 | 5 | 13/13 | 9.25e-16 |

## Combined workflow: results by case

Five test cases, with 7 comparisons per case.

| Case | Raw basis size k | Passed | Largest scaled error |
| --- | --- | --- | --- |
| uniform_k8 | 8 | 7/7 | 6.41e-15 |
| default_k8 | 8 | 7/7 | 9.52e-15 |
| uneven_k10 | 10 | 7/7 | 1.73e-14 |
| duplicates_k6 | 6 | 7/7 | 5.15e-15 |
| small_k5 | 5 | 7/7 | 1.92e-15 |

## Conversion only: individual checks

All nine comparisons use the main uniform_k8 case and the same R-generated B,S.

| Comparison | Passed | Scaled error |
| --- | --- | --- |
| Curve preserved | Yes | 0 |
| Penalty preserved | Yes | 2.01e-15 |
| All coefficient directions retained | Yes | 6.42e-15 |
| Fixed-effect space matches mgcv | Yes | 6.76e-15 |
| Random covariance factor matches mgcv | Yes | 2.08e-14 |
| Fixed-effect space matches captured lme inputs | Yes | 6.76e-15 |
| Random covariance factor matches captured lme inputs | Yes | 2.08e-14 |
| Curve preserved after coefficient mapping | Yes | 2.64e-15 |
| Penalty preserved after coefficient mapping | Yes | 1.69e-15 |

## Shared inputs

Data and knot fingerprints match the shared files for all five cases. The conversion check uses the same main observations and knots.
Main data MD5: `d7f3ce440112190e44cd54269ce48790`. Main knots MD5: `6665cb09426afb769ec771f896515c87`.
Cross-reference B and S scaled errors: 1.18e-15 and 1.2e-15.

## Detailed evidence

- [Basis construction checks](<spline_basis/outputs/shared/julia/checks.tsv>)
- [Conversion checks](<spline_conversion/outputs/shared/julia_conversion/validation.tsv>)
- [Combined workflow checks](<spline_workflow/outputs/shared/julia/checks.tsv>)
- [R input-capture checks](<spline_conversion/outputs/shared/r_reference/reference_checks.csv>)

## Prepared model

The [main Julia model inputs](<spline_workflow/outputs/shared/julia/uniform_k8>) contain y, X_fixed and Z, together with B, S, the coefficient map C, training settings and variance assumptions.

R reference versions: R version 4.6.1 (2026-06-24), mgcv 1.9.4, nlme 3.1.169.

## Sources

- [mgcv cubic B-splines and curvature penalty](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.bs.smooth.spec.html)
- [Basis dimensions](https://stat.ethz.ch/R-manual/R-devel/library/splines/html/splineDesign.html)
- [Spline conversion](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html)
