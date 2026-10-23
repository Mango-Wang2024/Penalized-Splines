# How the models are validated

[Project overview](README.md) · [Mathematical explanation](DERIVATION.md) · [Stage 1 results](STAGE1_RESULTS.md) · [Stage 2 results](STAGE2_RESULTS.md)

We check basis construction on its own, then compare the complete preparation of the mixed-model inputs. Julia uses BSplineKit and R uses mgcv, with the same saved Julia-generated observations and matched spline settings. The existing mathematical checks and cross-language criteria also validate the BSplineKit-based implementation.

## Stage 1: two comparisons

| Check | Shared information | What we compare |
| --- | --- | --- |
| Basis construction | Predictor values and spline settings | Independently constructed R and Julia $B,S$, before conversion |
| Combined workflow | Raw observations and settings | Complete R preparation against complete Julia preparation, including the required fixed/random transformation |

Construction and combined-workflow checks cover five cases. The captured R fitting-input check uses the main example.

Construction comparisons include basis values, derivatives, penalty, scaling, centering and evaluation beyond the training range. The combined workflow checks the final fixed/random representation and confirms that the curve, penalty and all coefficient directions are retained. Recording the R fitting inputs is also checked against a normal R fit, to ensure capture leaves its predictions and smoothing estimate unchanged.

## Stage 1 agreement

Equivalent converted models can have different column signs or rotations. We therefore compare the fixed-effect space, the random covariance factor $ZZ^\top$, and curves and penalties expressed in matching coefficient coordinates. The random variance multiplying this factor is estimated in Stage 2.

The numerical criterion is:

```math
\frac{\lVert A_{\mathrm{Julia}}-A_R\rVert}
{\max(1,\lVert A_{\mathrm{Julia}}\rVert,\lVert A_R\rVert)} < 10^{-8}.
```

For matrices, the norm is the square root of the sum of squared entries. The denominator handles small or zero values. Basis and workflow reports also record the largest entry difference as a diagnostic; it is not a separate pass criterion.

Independent mathematical examples check polynomial shapes, derivatives, curvature integration, centering, known conversion outcomes and invalid inputs. File fingerprints confirm that the reference comparisons used the intended data and knots. These checks supplement agreement between the two languages.

Numerical results remain in [STAGE1_RESULTS.md](STAGE1_RESULTS.md). Supporting mathematics and primary sources are collected in [DERIVATION.md](DERIVATION.md).

## Stage 2: REML fitting and smoothing estimation

R and Julia independently construct and fit all five cases from the same saved predictor and response values and spline settings. R uses `mgcv::gamm` / `nlme::lme`; Julia uses the low-level `MixedModels.LinearMixedModel` constructor with prepared `FeTerm` and `ReMat` containers, followed by `fit!(model; REML=true)`. The Stage 2 test repeats Stage 1 preparation in each branch and estimates smoothing through REML.

The model has one shared spline variance, $r^2$, and one residual variance, $\sigma^2$. Julia estimates $\theta=r/\sigma$ and obtains $\lambda=1/\theta^2$. R's spline variance is read independently from its fitted random-effect covariance, and its ratio to residual variance is checked against `gam$sp`. Both fits must converge.

### Error formulas and tolerances

Both errors measure agreement between R and Julia. For **smoothing strength**, we divide the difference by the larger estimate:

```math
E_{\lambda}=
\frac{|\lambda_{\mathrm{Julia}}-\lambda_R|}
{\max(|\lambda_{\mathrm{Julia}}|,|\lambda_R|,\epsilon_{\mathrm{mach}})}
\leq 10^{-5}.
```

This permits a relative difference of **0.001%**. The small guard $\epsilon_{\mathrm{mach}}$ is machine precision (about $2.22\times10^{-16}$), preventing division by zero. Variance comparisons use the same formula and threshold with $\sigma^2$ or $r^2$, checking the variance scale as well as its ratio.

For **fitted responses**, we combine differences at all 200 matching observations:

```math
E_{\hat y}=
\frac{\lVert\hat y_{\mathrm{Julia}}-\hat y_R\rVert_2}
{\max(1,\lVert\hat y_{\mathrm{Julia}}\rVert_2,\lVert\hat y_R\rVert_2)}
\leq 10^{-6}.
```

Here $\lVert v\rVert_2=\sqrt{\sum_i v_i^2}$. The denominator scales the difference by the size of the fitted-response vectors; the floor of 1 handles small values. This is one combined error across observations.

**Why two thresholds?** We chose these numerical acceptance thresholds for this project. The looser $10^{-5}$ parameter threshold accommodates differences between variance optimizers; the stricter $10^{-6}$ response threshold requires closer agreement in the resulting curve. These are comparison thresholds, separate from each fitter's stopping rules.

Supporting checks verify matching predictor, response and knot values, and exactly one shared spline variance in each model. Parameter symbols follow the derivation; in the results they denote fitted estimates.

Raw coefficient entries and package log-likelihood values are not used for cross-language agreement because coefficient coordinates and fixed-effect scaling can differ. Fitted responses are compared in the original observation coordinates. Results and individual checks are in [STAGE2_RESULTS.md](STAGE2_RESULTS.md).

### Checking the Julia input adapter

A dense spline matrix must be packaged consistently with `ReMat`'s native storage and shared variance. The adapter uses orthogonal changes of coordinates to make the random design diagonal while preserving the Gaussian covariance model and REML objective. It does not change the Stage 1 spline or add variance parameters. Predictions and random coefficients are returned to their original coordinates after fitting. The container's levels index transformed spline directions.

Independent Julia tests verify the adapter's covariance identity and round trip, compare its fitted solution and REML objective with a direct marginal Gaussian calculation, and check invariance to rotations of spline directions. They also cover zero estimated spline variance and invalid inputs. This calculation evaluates the REML-fitted solution; it is not a separate fitting stage with fixed $\lambda$.

Sources: [mgcv::gamm](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gamm.html), [nlme covariance extraction](https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/pdMatrix.html), [MixedModels estimation](https://juliastats.org/MixedModels.jl/stable/optimization/), and the pinned [constructor](https://github.com/JuliaStats/MixedModels.jl/blob/v5.9.0/src/linearmixedmodel.jl) and [ReMat implementation](https://github.com/JuliaStats/MixedModels.jl/blob/v5.9.0/src/remat.jl).
