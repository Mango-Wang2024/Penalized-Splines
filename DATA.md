# Data and test cases

[Project overview](README.md) · [Results by case](STAGE1_RESULTS.md)

We use simulated observations generated once in Julia. R and Julia read the same saved values, so the comparisons share their data. The random seed is 42.

## Observations

Every case has 200 observations. The main dataset contains 198 randomly distributed predictor values between −2 and 2, plus the two endpoints. Each observation has:

- **`x`:** the predictor value $x$.
- **`truth`:** the known curve, $f(x)=\sin(1.7x)+0.15x$.
- **`y`:** the observed response $y=f(x)+\varepsilon$, with normally distributed noise of standard deviation 0.35.

The basis depends on predictor values and spline settings. Stage 2 fits the observed responses $y$; the known curve helps interpret the simulation and is not supplied to the fitter.

## Five test cases

The case names identify different data or spline settings. $k$ is the number of raw basis functions; a full cubic knot vector has $k+4$ entries.

| Case | $k$ | Full knot entries | Purpose |
| --- | --- | --- | --- |
| `uniform_k8` | 8 | 12 | Main case with an explicitly supplied uniform knot vector |
| `default_k8` | 8 | 12 | Check that R and Julia independently obtain the same default knots |
| `uneven_k10` | 10 | 14 | Check uneven knot spacing and a larger basis |
| `duplicates_k6` | 6 | 10 | Check repeated predictor values: 40 distinct values, each repeated five times |
| `small_k5` | 5 | 9 | Check a smaller basis |

Four cases share exactly the same observations. `duplicates_k6` has its own predictor values and noisy responses. The explicit and default eight-function cases use the same knot vector but obtain it differently, allowing that step to be checked independently.

Each case also has 61 evaluation points from −2.7 to 2.7 for checking basis values within and beyond the training range.

## Connection to the results

Basis construction has **13 comparisons per case**, giving 65 across five cases. The combined workflow has **7 comparisons per case**, giving 35, and includes the required fixed/random preparation. The main `uniform_k8` case also provides the captured R fitting inputs.

Stage 2 reuses these same five cases and observed responses. Each fit estimates smoothing automatically; the [fitting results](STAGE2_RESULTS.md) compare R and Julia without changing the data or knots.

Saved observations, knots and evaluation points are under [data/inputs](data/inputs). Sources: the [Julia data definition](data/generate.jl), [recorded case settings](data/inputs/cases.tsv), and [cubic B-spline dimensions](https://stat.ethz.ch/R-manual/R-devel/library/splines/html/splineDesign.html).
