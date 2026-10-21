# Penalized splines as mixed models in Julia

We construct cubic B-splines in Julia, prepare their mixed-model inputs, and estimate smoothing with **MixedModels.jl**. R's **mgcv** and **nlme** provide the validation reference. Future work will connect this approach to Unfold for EEG analysis.

The project follows **model preparation → fitting → prediction**. Stages 1 and 2 are implemented and checked against R; Stage 3 is next.

[DATA.md](DATA.md) explains the observations and test cases; [DERIVATION.md](DERIVATION.md) explains the mathematics; [VALIDATION.md](VALIDATION.md) defines the comparison criteria. [STAGE1_RESULTS.md](STAGE1_RESULTS.md) records input-preparation checks, and [STAGE2_RESULTS.md](STAGE2_RESULTS.md) records fitting results.

## Stage 1: prepare the model inputs

We use two comparisons to check preparation independently of fitting:

1. **Basis construction only:** R and Julia independently construct $B,S$ from the same predictor values and spline settings. We compare them before fixed/random preparation.
2. **Combined workflow:** each language constructs its spline and prepares $X_{\mathrm{fixed}},Z$ from the same raw inputs. We compare the final model representations.

### Comparison flow

![Stage 1 compares basis construction and the complete preparation of model inputs.](assets/stage1_two_comparisons.png)

## Stage 2: fit and estimate smoothing

The Stage 2 test runs the validated Stage 1 preparation again, then fits both models with **REML**, which estimates variances while accounting for the fitted fixed effects. Both branches start from the same saved observations and spline settings:

- **R:** `mgcv::gamm` prepares the spline and calls `nlme::lme`.
- **Julia:** our Stage 1 functions prepare $y,X_{\mathrm{fixed}},Z$. An adapter forms `FeTerm` and `ReMat` containers for `MixedModels.LinearMixedModel`, then `fit!` performs REML estimation.

Each fit estimates fixed coefficients $b$, penalized spline coefficients $r$, noise variance $\sigma^2$, and shared spline variance $\tau^2$. We keep the parameter symbols from the derivation; here they refer to fitted estimates. Using the fitted coefficients gives the training responses $\hat y$:

```math
\hat y=X_{\mathrm{fixed}}b+Zr,
\qquad
\lambda=\frac{\sigma^2}{\tau^2}.
```

$\lambda$ is derived from the same variance fit, rather than estimated in an additional fitting step. A larger ratio means stronger smoothing. The hat on $\hat y$ distinguishes fitted responses from observed $y$.

### Comparison flow

![Stage 2 repeats preparation in separate R and Julia branches, fits each model with REML, then compares smoothing strength and fitted responses.](assets/stage2_fitting_flow.svg)

The two main comparisons are **smoothing strength $\lambda$** and **fitted responses $\hat y$** across five cases. Variance estimates remain supporting checks. Coefficient coordinates can differ between equivalent representations, so raw coefficient entries are not the cross-language comparison. The [results](STAGE2_RESULTS.md) show one summary table and an overlaid curve plot, with detailed estimates linked there.

Sources: [mgcv fitting](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gamm.html), [MixedModels constructors](https://juliastats.org/MixedModels.jl/stable/constructors/), and [MixedModels estimation](https://juliastats.org/MixedModels.jl/stable/optimization/). The fitting environment pins MixedModels 5.9.0 because the adapter uses its low-level container API.

## Model and inputs

A spline combines simple curve shapes, and a curvature penalty controls their bending. Our main example uses 200 simulated observations, $k=8$ raw basis functions and 12 full knot entries. Centering separates the spline from the intercept, leaving seven basis columns.

| Input | Meaning | Main dimensions |
| --- | --- | --- |
| $y$ | Observed responses | 200 values |
| $B$, $S$ | Spline basis and matching penalty | 200 × 7; 7 × 7 |
| $X_{\mathrm{fixed}}$ | Intercept and unpenalized spline component | 200 × 2 |
| $Z$ | Penalized spline components sharing variance $\tau^2$ | 200 × 6 |

Training knots and transformations are retained for prediction. The repository includes the linked evidence tables, captured R inputs and main Julia model inputs.

## Connection to the R fitting inputs

The [referenced `gamm.r` call](https://github.com/cran/mgcv/blob/1b6a4c8374612da27e36420b4459e93acb183f2d/R/gamm.r#L1398) passes the prepared representation to `lme`:

| R argument | Meaning in our project |
| --- | --- |
| `fixed.formula` | Describes the response and components represented by $X_{\mathrm{fixed}}$ |
| `rand` | Describes the $Z$ components, grouping and shared variance |
| `data` | Contains the observations and numerical matrices |
| `correlation` | Residual relationships; our example assumes independent errors |
| `weights` | Residual variability; our example assumes equal variability |
| `method` | Fitting method; we use REML |
| `control` | Numerical stopping rules and limits |

The [captured R inputs](spline_workflow/outputs/shared/r_reference/lme_inputs/INPUTS_EXPLAINED.txt) provide a concrete example; the [official `lme` documentation](https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/lme.html) explains the arguments.

## Next steps

- **Stage 3 — prediction:** evaluate the fitted model at new predictor values using its training knots, transformations and learned coefficients, then compare predictions with R.
- **Unfold integration:** use Unfold or UnfoldMixedModels where the required variance structure is supported with comparable efficiency, starting with controlled examples before EEG applications.
