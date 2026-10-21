# Penalized splines as mixed models in Julia

We construct cubic B-splines in Julia, prepare their mixed-model inputs, and estimate smoothing with **MixedModels.jl**. R's **mgcv** and **nlme** provide the validation reference. Future work will connect this approach to Unfold for EEG analysis.

The project follows **model preparation → fitting → prediction**. Stages 1 and 2 are implemented and checked against R; Stage 3 is planned.

[DATA.md](DATA.md) explains the observations and test cases; [DERIVATION.md](DERIVATION.md) explains the mathematics; [VALIDATION.md](VALIDATION.md) defines the comparison criteria. [STAGE1_RESULTS.md](STAGE1_RESULTS.md) records input-preparation checks, and [STAGE2_RESULTS.md](STAGE2_RESULTS.md) records fitting results.

The main dimensions below refer to our **200-observation example with $k=8$ raw basis functions**. Matrix dimensions are written as rows × columns; other test cases use different basis sizes.

## Stage 1: prepare the model inputs

**Purpose:** turn predictor values and spline settings into matrices that a mixed model can use.

| Role | Quantity and meaning | Main dimensions |
| --- | --- | --- |
| Input | $x$: training predictor values | 200 values |
| Input | Full knot vector: positions defining the spline pieces | 12 values |
| Constructed | $B,S$: centered spline basis and matching curvature penalty | 200 × 7; 7 × 7 |
| Output | $X_{\mathrm{fixed}}$: fixed-effect design | 200 × 2 |
| Output | $Z$: penalized spline design | 200 × 6 |

**Method and settings:** Julia uses our basis-construction and model-preparation functions; R uses mgcv as an independent reference. Both use cubic B-splines (degree 3), $k=8$ in the main example, a supplied or automatically calculated knot vector, and a curvature penalty that measures squared second derivatives.

Centering separates the spline from the intercept, reducing eight raw basis columns to seven. $X_{\mathrm{fixed}}=[X,X_0]$ combines the intercept with the unpenalized spline component; $Z$ contains the six penalized components. Training knots, centering, penalty scaling and coefficient transformations are retained for prediction. Observed responses $y$ are used for fitting in Stage 2.

We check preparation independently of fitting:

1. **Basis construction only:** R and Julia independently construct $B,S$ from the same predictor values and spline settings. We compare them before fixed/random preparation.
2. **Combined workflow:** each language constructs its spline and prepares $X_{\mathrm{fixed}},Z$ from the same raw inputs. We compare the final model representations.

### Comparison flow

![Stage 1 compares basis construction and the complete preparation of model inputs.](assets/stage1_two_comparisons.png)

## Stage 2: fit and estimate smoothing

**Purpose:** use the prepared model and observed responses to learn the curve and its smoothing strength. The Stage 2 test repeats Stage 1 preparation in each branch, starting from the same saved observations and settings.

| Role | Quantity and meaning | Main dimensions |
| --- | --- | --- |
| Input | $y$: observed training responses | 200 values |
| Input | $X_{\mathrm{fixed}},Z$: designs prepared in Stage 1 | 200 × 2; 200 × 6 |
| Output | $b$: fitted fixed coefficients | 2 values |
| Output | $r$: fitted penalized spline coefficients | 6 values |
| Output | $\sigma^2,\tau^2$: noise variance and shared spline variance | One value each |
| Output | $\lambda$: smoothing strength derived from the variances | One value |
| Output | $\hat y$: fitted responses at the training predictor values | 200 values |

**Tools:** R's `mgcv::gamm` prepares the spline and calls `nlme::lme`. Julia prepares its own matrices, packages them in `FeTerm` and `ReMat` containers, and fits with `MixedModels.LinearMixedModel` and `fit!`.

**Essential settings:** REML estimates variances while accounting for fitted fixed effects. All six penalized coefficients share one Gaussian variance $\tau^2$; noise is independent Gaussian with variance $\sigma^2$, and observations have equal weights. Smoothing is estimated automatically rather than fixed beforehand.

Using fitted coefficients and variances:

```math
\hat y=X_{\mathrm{fixed}}b+Zr,
\qquad
\lambda=\frac{\sigma^2}{\tau^2}.
```

Parameter symbols follow the derivation; here they denote fitted estimates. The hat on $\hat y$ distinguishes fitted responses from observed $y$. A larger $\lambda$ means stronger smoothing; it is derived from the variance fit, without an additional fitting step.

### Comparison flow

![Stage 2 repeats preparation in separate R and Julia branches, fits each model with REML, then compares smoothing strength and fitted responses.](assets/stage2_fitting_flow.svg)

The two main comparisons are **smoothing strength $\lambda$** and **fitted responses $\hat y$** across five cases. Variance estimates remain supporting checks. Equivalent representations can have different coefficient coordinates, so raw coefficient entries are not the cross-language comparison. The [results](STAGE2_RESULTS.md) include a summary table, an overlaid curve plot and links to detailed estimates.

## Stage 3: predict at new predictor values — planned

**Purpose:** apply the fitted curve to new predictor values. Let $m$ be the number of new observations.

| Role | Quantity and meaning | Main dimensions |
| --- | --- | --- |
| Input | $x_{\mathrm{new}}$: new predictor values | $m$ values |
| Input | Training knots and saved Stage 1 transformations | 12 knot values; transformations reused |
| Input | Fitted coefficients $b,r$ from Stage 2 | 2 values; 6 values |
| Constructed | $X_{\mathrm{fixed,new}},Z_{\mathrm{new}}$: designs evaluated at new values | $m$ × 2; $m$ × 6 |
| Output | $\hat y_{\mathrm{new}}$: predicted responses | $m$ values |

The planned Julia prediction step will reuse the training knots, centering, penalty scaling and coefficient transformations, then apply the fitted coefficients. The learned smoothing is already reflected in those coefficients. We will compare these predictions with R at the same new predictor values; observed responses at those points are not required to produce predictions.

After this, we plan to integrate with Unfold or UnfoldMixedModels where the required variance structure is supported with comparable efficiency, starting with controlled examples before EEG applications.

## Sources and package interfaces

Sources: [mgcv fitting](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gamm.html), [MixedModels constructors](https://juliastats.org/MixedModels.jl/stable/constructors/), and [MixedModels estimation](https://juliastats.org/MixedModels.jl/stable/optimization/). The [original R fitting call](https://github.com/cran/mgcv/blob/1b6a4c8374612da27e36420b4459e93acb183f2d/R/gamm.r#L1398) and [captured R inputs](spline_workflow/outputs/shared/r_reference/lme_inputs/INPUTS_EXPLAINED.txt) document how R packages the same fixed/random model. The Julia fitting environment pins MixedModels 5.9.0 because the adapter uses its low-level container API.
