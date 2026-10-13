# Research plan

## Goal

Implement penalized spline fitting in Julia and connect it to Unfold. Begin with one continuous predictor, one spline penalty and independent, normally distributed errors.

The final implementation should run without R. During development, mgcv provides a reference for comparison.

## Evidence needed

| Claim | Evidence |
| --- | --- |
| Mathematical equivalence | A derivation showing that the curve and penalty are preserved |
| Correct implementation | Numerical checks and independently solved examples |
| Agreement with mgcv | Comparisons using matched basis, penalty and fitting settings |

A derivation does not prove that the code is correct. A valid conversion of supplied matrices does not prove that Julia independently constructs the same spline as mgcv.

## Stage 1: Construct and transform inputs

Build the spline basis `B` and penalty `S`. Separate unpenalized columns `X0` from penalized columns `Z`.

Check:

```text
B z = X0 alpha + Z u
z' S z = u' u
```

Here, `z` weights the original spline shapes. `alpha` and `u` weight the unpenalized and penalized shapes after conversion. No coefficient directions should be lost.

The complete fixed matrix is `[X X0]`, including ordinary predictors `X`. Compare fixed-effect spaces and `Z Z'` rather than requiring identical eigenvectors. Match knots, constraints and penalty scaling when checking mgcv compatibility.

Stage 1 has three separate checks:

1. Basis only: compare independently constructed R and Julia B,S from shared inputs.
2. Conversion only: give both converters the exact same R-generated B,S.
3. Combined workflow: compare R construction/conversion with Julia construction/conversion from the same raw inputs.

See the [Stage 1 summary](../STAGE1_RESULTS.md). The [combined Julia outputs](../spline_workflow/outputs/shared/julia/uniform_k8) provide y, the full fixed matrix and Z for Stage 2.

## Stage 2: Fit the model

### Checkpoint A: Fixed smoothing strength

Choose `lambda` and compare the transformed fit with an independent solution of:

```text
min ||y - X beta - B z||² + lambda z' S z
```

`y` contains observed responses and `beta` weights the ordinary predictors. Compare fitted values, objective values and identifiable coefficients in the same coordinates. This check can run entirely in Julia.

### Checkpoint B: Estimated smoothing strength

Estimate the variance parameters using ML or REML, established likelihood-based methods. For the normalized single-penalty model:

```text
Var(u) = sigma_u² I
lambda = sigma_error² / sigma_u²
```

`I` is the identity matrix. All penalized coefficients share one variance and have zero correlations. Compare with a matched R fit or an independently validated likelihood implementation.

## Stage 3: Predict new data

Reuse the training knots, constraints and transformation. Compare final predictions through both representations:

```text
yhat_new = X_new betahat + B_new zhat
         = X_new betahat + X0_new alphahat + Z_new uhat
```

The suffix `new` identifies matrices for new observations. The suffix `hat` identifies fitted coefficients. Include the fitted spline random coefficients `uhat`, which form part of the curve. R predictions provide an additional compatibility check.

## Integration

Use the exploratory notebook as an implementation reference. Its results remain unvalidated. Confirm that the selected Unfold fitting route supports the required shared-variance structure before integration.

## Sources

- [Conversion and prediction](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html)
- [Penalized fitting objective](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/magic.html)
- [Variance-based smoothing estimation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gamm.html)
- [Reference source revision](https://github.com/cran/mgcv/blob/1b6a4c8374612da27e36420b4459e93acb183f2d/R/gamm.r#L1398)
