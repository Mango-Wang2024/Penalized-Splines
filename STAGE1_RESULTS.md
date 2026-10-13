# Stage 1 validation

Three separate checks passed for the tested natural cubic spline model. Inputs were generated once in Julia.

| Check | Passed | Largest scaled error |
| --- | --- | --- |
| Basis construction only | 65/65 | 6.14e-16 |
| Conversion only | 9/9 | 6.97e-14 |
| Combined workflow | 35/35 | 7.08e-14 |

The mathematical tests completed before these comparisons. R input-capture checks passed: 7/7. Comparison threshold: 1e-8.

## What each check establishes

1. **Basis:** R and Julia independently construct B and S from the same data and settings. No conversion is run.
2. **Conversion:** both converters use the exact same R-generated B and S. Julia is compared with mgcv and the actual lme inputs.
3. **Workflow:** each language constructs and converts its own spline from the shared raw inputs. Final model representations are compared.

The basis and workflow checks cover five cases. The captured lme conversion check uses the main eight-knot case.

## Shared inputs

All basis and workflow input fingerprints match the shared files. The conversion check uses the same main data and knots.
Main data MD5: `d7f3ce440112190e44cd54269ce48790`. Main knots MD5: `c411e2e123b06233043b39fdb5d783c4`.
Cross-reference B and S scaled errors: 7.93e-16 and 6.53e-16.

## Reports

- [Basis construction](<spline_basis/outputs/shared/RESULTS.md>)
- [Conversion](<spline_conversion/outputs/shared/julia_conversion/RESULTS_EXPLAINED.txt>)
- [Combined workflow](<spline_workflow/outputs/shared/RESULTS.md>)

## Stage 2 inputs

The [main Julia model inputs](<spline_workflow/outputs/shared/julia/uniform_k8>) contain y, X_fixed and Z. Matching B, S, C, knots, centering H, penalty scale and variance assumptions are retained.

Stage 1 prepares the model. Julia coefficient fitting, smoothing estimation and fitted predictions remain later stages.
The R input-capture example performs a reference fit. Its fitted parameters are not Julia estimates.
Reference versions: R version 4.6.1 (2026-06-24), mgcv 1.9.4, nlme 3.1.169. The historical source revision is not automatically reproduced.
