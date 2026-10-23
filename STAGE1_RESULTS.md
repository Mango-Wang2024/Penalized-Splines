# Stage 1 results

Saved comparisons for cubic B-splines with one curvature penalty: Julia uses BSplineKit; R uses mgcv bs="bs", m=c(3,2). All comparisons use shared Julia-generated inputs.

Stage 1 checks basis construction and the combined workflow. The required fixed/random preparation is part of that workflow.

## Overall results

| Check | Passed | Largest scaled error |
| --- | --- | --- |
| Basis construction only | 65/65 | 1.54e-15 |
| Combined workflow | 35/35 | 5.13e-14 |

Mathematical tests completed before these comparisons. R fitting-input checks passed: 7/7. Comparison threshold: 1e-8.

## Basis construction only: results by case

Five test cases, with 13 comparisons per case.

| Case | Raw basis size $k$ | Passed | Largest scaled error |
| --- | --- | --- | --- |
| uniform_k8 | 8 | 13/13 | 3.7e-16 |
| default_k8 | 8 | 13/13 | 9.26e-16 |
| uneven_k10 | 10 | 13/13 | 1.54e-15 |
| duplicates_k6 | 6 | 13/13 | 1.52e-15 |
| small_k5 | 5 | 13/13 | 9.06e-16 |

## Combined workflow: results by case

Five test cases, with 7 comparisons per case.

| Case | Raw basis size $k$ | Passed | Largest scaled error |
| --- | --- | --- | --- |
| uniform_k8 | 8 | 7/7 | 1.31e-14 |
| default_k8 | 8 | 7/7 | 1.84e-14 |
| uneven_k10 | 10 | 7/7 | 5.13e-14 |
| duplicates_k6 | 6 | 7/7 | 5.56e-15 |
| small_k5 | 5 | 7/7 | 2.76e-15 |

## Shared inputs

Data and knot fingerprints match the shared files for all five cases. The captured R fitting inputs use the same main observations and knots.
Main data MD5: `d7f3ce440112190e44cd54269ce48790`. Main knots MD5: `6665cb09426afb769ec771f896515c87`.
Cross-reference $B$ and $S$ scaled errors: 1.18e-15 and 1.2e-15.

## Detailed evidence

- [Basis construction checks](<spline_basis/outputs/shared/julia/checks.tsv>)
- [Combined workflow checks](<model_preparation/outputs/shared/julia/checks.tsv>)
- [R fitting-input checks](<model_preparation/outputs/shared/r_reference/lme_inputs/reference_checks.csv>)

## Prepared model

The [main Julia model inputs](<model_preparation/outputs/shared/julia/uniform_k8>) contain $y$, $X_{\mathrm{fixed}}$ and $Z$, together with $B$, $S$, the coefficient map $C$, training settings and variance assumptions.

R reference versions: R version 4.6.1 (2026-06-24), mgcv 1.9.4, nlme 3.1.169.

Julia 1.13.1; BSplineKit 0.19.3. Both stages use the pinned root Julia environment.

## Sources

- [mgcv cubic B-splines and curvature penalty](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.bs.smooth.spec.html)
- [Basis dimensions](https://stat.ethz.ch/R-manual/R-devel/library/splines/html/splineDesign.html)
- [Fixed/random preparation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html)
