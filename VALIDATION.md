# How Stage 1 is validated

[Project overview](README.md) · [Mathematical explanation](DERIVATION.md) · [Recorded results](STAGE1_RESULTS.md)

We check basis construction on its own, then compare the complete preparation of the mixed-model inputs. Both languages use the same saved Julia-generated observations and matched spline settings.

## Two comparisons

| Check | Shared information | What we compare |
| --- | --- | --- |
| Basis construction | Predictor values and spline settings | Independently constructed R and Julia `B,S`, before conversion |
| Combined workflow | Raw observations and settings | Complete R preparation against complete Julia preparation, including the required fixed/random transformation |

Construction and combined-workflow checks cover five cases. The captured-input conversion check uses the main example.

Construction comparisons include basis values, derivatives, penalty, scaling, centering and evaluation beyond the training range. The combined workflow checks the final fixed/random representation and confirms that the curve, penalty and all coefficient directions are retained. Recording the R fitting inputs is also checked against a normal R fit, to ensure capture leaves its predictions and smoothing estimate unchanged.

## What counts as agreement

Equivalent converted models can have different column signs or rotations. We therefore compare the fixed-effect space, the random covariance factor $ZZ^\top$, and curves and penalties expressed in matching coefficient coordinates. The random variance multiplying this factor will be estimated in Stage 2.

The numerical criterion is:

```math
\frac{\lVert A_{\mathrm{Julia}}-A_R\rVert}
{\max(1,\lVert A_{\mathrm{Julia}}\rVert,\lVert A_R\rVert)} < 10^{-8}.
```

For matrices, the norm is the square root of the sum of squared entries. The denominator handles small or zero values. Basis and workflow reports also record the largest entry difference as a diagnostic; it is not a separate pass criterion.

Independent mathematical examples check polynomial shapes, derivatives, curvature integration, centering, known conversion outcomes and invalid inputs. File fingerprints confirm that the reference comparisons used the intended data and knots. These checks supplement agreement between the two languages.

Numerical results remain in [STAGE1_RESULTS.md](STAGE1_RESULTS.md). Supporting mathematics and primary sources are collected in [DERIVATION.md](DERIVATION.md).
