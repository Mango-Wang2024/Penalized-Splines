# Data and test cases

[Project overview](README.md) · [Stage 1 results](STAGE1_RESULTS.md) · [Stage 2 results](STAGE2_RESULTS.md)

Julia generates the observations once; both languages read the same saved values. The example has **200 observations, two predictors and one response**:

```math
\begin{aligned}
f_1(x_1)&=\sin(1.7x_1)+0.15x_1,\\
f_2(x_2)&=0.5\cos(1.2x_2),\\
y&=f_1(x_1)+f_2(x_2)+\varepsilon,
\qquad\varepsilon\sim N(0,0.35^2).
\end{aligned}
```

The predictors are generated independently between −2 and 2. $x_1$ is sorted and includes both endpoints; $x_2$ retains its independent row order. Each row pairs its own $x_1,x_2,y$. Seeds are 42 for $x_1$, 43 for $x_2$, and 44 for noise.

`data.csv` contains `x1`, `x2`, `y`, the two known component values, and their sum `truth`. Only predictor values and noisy responses are used for fitting. Known curves are retained to interpret the simulation.

## Five cases

Each predictor has its own settings. $k_j$ is its raw basis size; a full cubic knot vector has $k_j+4$ entries.

| Case | $k_1$ | $k_2$ | Purpose |
| --- | ---: | ---: | --- |
| `uniform_k8` | 8 | 8 | Main example; explicitly supplied uniform knots |
| `default_k8` | 8 | 8 | Each language independently calculates the default knots |
| `uneven_k10` | 10 | 6 | Uneven knots for $x_1$; different basis sizes across components |
| `duplicates_k6` | 6 | 8 | Repeated $x_1$ values: 40 distinct values, each repeated five times |
| `small_k5` | 5 | 5 | Smaller bases for both components |

Four cases share the same observations. The repeated-value case changes $x_1$ and the corresponding response, while retaining $x_2$ and the same noise realization. Explicit and default uniform cases obtain the same knots through different paths.

## Saved information

[Each case](data/inputs/) contains its observations, a `terms.tsv` list of predictors and spline settings, and separate knot and evaluation files for each predictor. The term list lets preparation and fitting accept a variable number of predictors without changing their functions.

Each predictor also has 61 evaluation points from −2.7 to 2.7 for basis and extrapolation checks. A shared 31 × 31 grid within the training ranges supplies 961 predictor pairs for comparing fitted response surfaces. This grid is used after fitting and does not train the model.

Stage 1 has 13 basis comparisons per component and case: **130 in total**. Complete preparation has one joint fixed-space check and six checks per component: **65 in total**. The main case also captures the actual R fitting inputs. Stage 2 fits all five cases jointly and compares both smoothing strengths and the resulting responses.

Sources: [Julia data definition](data/generate.jl), [saved case list](data/inputs/cases.tsv), and [cubic B-spline dimensions](https://stat.ethz.ch/R-manual/R-devel/library/splines/html/splineDesign.html).
