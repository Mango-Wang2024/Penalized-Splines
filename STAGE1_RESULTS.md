# Stage 1 results

R and Julia independently construct each component and assemble the additive model from the same saved inputs.

## Overall results

| Check | Passed | Largest scaled error |
| --- | ---: | ---: |
| Basis construction only | 130/130 | 2.08e-15 |
| Combined additive preparation | 65/65 | 5.13e-14 |

Captured R fitting-input checks: **10/10 passed**. All preparation errors are below $10^{-8}$.
The basis references and captured fitting inputs agree; their largest scaled difference is 0.

## Results by case

| Case | Components | Basis checks | Preparation checks | Largest scaled error |
| --- | ---: | ---: | ---: | ---: |
| uniform_k8 | 2 | 26/26 | 13/13 | 1.31e-14 |
| default_k8 | 2 | 26/26 | 13/13 | 1.84e-14 |
| uneven_k10 | 2 | 26/26 | 13/13 | 5.13e-14 |
| duplicates_k6 | 2 | 26/26 | 13/13 | 6.62e-15 |
| small_k5 | 2 | 26/26 | 13/13 | 2.76e-15 |

## Prepared model

The main two-predictor case produces $X_{\mathrm{fixed}}$ with dimensions 200 × 3 and two penalized designs $Z_1,Z_2$, each 200 × 6. Each component retains its own basis, penalty, knots and coefficient transformation.

## Detailed evidence

- [Basis comparisons](spline_basis/outputs/shared/julia/checks.tsv)
- [Additive preparation comparisons](model_preparation/outputs/shared/julia/checks.tsv)
- [Captured R inputs](model_preparation/outputs/shared/r_reference/lme_inputs/INPUTS_EXPLAINED.txt)
- [Captured-input checks](model_preparation/outputs/shared/r_reference/lme_inputs/reference_checks.csv)

Sources: [mgcv construction](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smoothCon.html) and [mixed-model preparation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html).
