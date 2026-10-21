# How the models are validated

[Project overview](README.md) · [Mathematical explanation](DERIVATION.md) · [Stage 1 results](STAGE1_RESULTS.md) · [Stage 2 results](STAGE2_RESULTS.md)

We check basis construction on its own, then compare the complete preparation of the mixed-model inputs. Both languages use the same saved Julia-generated observations and matched spline settings.

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

## Stage 2: automatic REML fitting

R and Julia independently construct and fit all five cases from the same saved predictor and response values and spline settings. R uses `mgcv::gamm` / `nlme::lme`; Julia uses the low-level `MixedModels.LinearMixedModel` constructor with prepared `FeTerm` and `ReMat` containers, followed by `fit!(model; REML=true)`. The Stage 2 test repeats Stage 1 preparation in each branch; $\lambda$ is estimated automatically.

The model has one shared spline variance, $\tau^2$, and one residual variance, $\sigma^2$. Julia estimates $\theta=\tau/\sigma$ and obtains $\lambda=1/\theta^2$. R's spline variance is read independently from its fitted random-effect covariance, and its ratio to residual variance is checked against `gam$sp`. Both fits must converge.

We compare:

- **$\hat\lambda$, $\hat\sigma^2$, and $\hat\tau^2$:** absolute difference divided by the larger magnitude, with tolerance $10^{-5}$.
- **All fitted responses $\hat y$:** the scaled norm above, with tolerance $10^{-6}$.
- **Inputs and model structure:** matching predictor, response and knot values, and exactly one shared spline variance in each model.

The fitting tolerances allow small differences between numerical optimizers. Raw coefficient entries and package log-likelihood values are not used for cross-language agreement because coefficient coordinates and fixed-effect scaling can differ. Fitted responses are compared in the original observation coordinates. Since $\hat\lambda=\hat\sigma^2/\hat\tau^2$, the smoothing and variance comparisons describe related outputs of the same fit. Results and individual checks are in [STAGE2_RESULTS.md](STAGE2_RESULTS.md).

### Checking the Julia input adapter

A dense spline matrix must be packaged consistently with `ReMat`'s native storage and shared variance. The adapter uses orthogonal changes of coordinates to make the random design diagonal while preserving the Gaussian covariance model and REML objective. It does not change the Stage 1 spline or add variance parameters. Predictions and random coefficients are returned to their original coordinates after fitting. The container's levels index transformed spline directions.

Independent Julia tests verify the adapter's covariance identity and round trip, compare its fitted solution and REML objective with a direct marginal Gaussian calculation, and check invariance to rotations of spline directions. They also cover zero estimated spline variance and invalid inputs. This calculation evaluates the automatically fitted solution; it is not a separate fitting stage with fixed $\lambda$.

Sources: [mgcv::gamm](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gamm.html), [nlme covariance extraction](https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/pdMatrix.html), [MixedModels estimation](https://juliastats.org/MixedModels.jl/stable/optimization/), and the pinned [constructor](https://github.com/JuliaStats/MixedModels.jl/blob/v5.9.0/src/linearmixedmodel.jl) and [ReMat implementation](https://github.com/JuliaStats/MixedModels.jl/blob/v5.9.0/src/remat.jl).
