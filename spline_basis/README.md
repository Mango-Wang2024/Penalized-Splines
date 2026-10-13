# Check 1: Basis construction only

This check compares Julia's independently constructed natural cubic spline basis `B` and roughness penalty `S` with mgcv's construction. No conversion code is imported or run.

```text
Shared predictor values and spline settings
        -> R constructs B and S
        -> Julia constructs B and S
        -> Compare the matrices
```

The main inputs are `../data/inputs/uniform_k8/data.csv` and `knots.csv`. Julia reads the shared inputs, never R's matrices, to construct its own results. Response values do not affect this construction.

From this folder:

```sh
julia --startup-file=no test/runtests.jl
sh run.sh
```

If outputs already exist, use a new directory:

```sh
sh run.sh outputs/run2
```

| File | Purpose |
| --- | --- |
| `src/NaturalCubic.jl` | Construct and evaluate natural cubic splines |
| `scripts/reference.R` | Construct R reference matrices only |
| `scripts/compare.jl` | Compare basis, penalty, scaling and centering |
| `test/runtests.jl` | Check examples with known mathematical answers |
| `DERIVATION.md` | Explain the construction equations |

The comparison has 13 checks per case, 65 across five cases. There are also 28 mathematical tests. Cases cover regular knots, default quantile knots, uneven knots, repeated predictor values and three knots.

See the [results](outputs/shared/RESULTS.md), [detailed checks](outputs/shared/julia/checks.tsv) and [construction equations](DERIVATION.md).

Conversion-only validation is in [spline_conversion](../spline_conversion/README.md). The combined comparison is in [spline_workflow](../spline_workflow/README.md). Julia fitting and smoothing estimation are later stages.

Sources: [natural cubic splines](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.cr.smooth.spec.html) and [mgcv scaling and constraints](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smoothCon.html). Reference functions were inspected in mgcv 1.9.4. Each run records its versions.
