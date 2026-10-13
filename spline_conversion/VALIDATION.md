# Check 2: Conversion validation

Shared-data comparison completed on **13 October 2026** using R 4.6.1, mgcv 1.9.4, nlme 3.1.169 and Julia 1.13.1.

## Checks

| Check | Result |
| --- | --- |
| Julia examples with known mathematical answers | 44/44 passed |
| R input capture and reference conversion | 7/7 passed |
| Julia conversion compared with actual R outputs | 9/9 passed |

Recording R's inputs left its fitted predictions and smoothing estimate unchanged. The largest scaled difference in the live Julia comparison was **6.97 × 10⁻¹⁴**.

Both projects used the same Julia-generated data and knot files. The input hashes are recorded in the R metadata.

## Example dimensions

| Object | Dimensions | Meaning |
| --- | --- | --- |
| Spline basis | 200 × 7 | Curve shapes evaluated at 200 observations |
| Penalty | 7 × 7 | Roughness rule for those shapes |
| Full fixed matrix | 200 × 2 | Intercept and one unpenalized spline component |
| Random matrix | 200 × 6 | Six penalized spline components |

The requested basis size was 8. A centering constraint reduced it to 7 columns.

## What this establishes

Both languages convert the exact same R-generated B and S. Julia preserves the curve and penalty and agrees with mgcv's converted model. Different eigenvector signs do not affect that conclusion.

Independent construction is checked in [spline_basis](../spline_basis/outputs/shared/RESULTS.md). The complete construction/conversion comparison is checked in [spline_workflow](../spline_workflow/outputs/shared/RESULTS.md). Julia fitting, smoothing estimation and final fitted predictions remain untested. These results do not establish general robustness or historical-version compatibility.

Detailed reports:

- [R input explanation](outputs/shared/r_reference/INPUTS_EXPLAINED.txt)
- [Julia comparisons](outputs/shared/julia_conversion/RESULTS_EXPLAINED.txt)
