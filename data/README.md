# Shared Julia data

Both validation projects use **`inputs/uniform_k8/data.csv`** and its **`knots.csv`** for the main example.

Julia generates the data once. R reads the saved values and does not simulate another dataset.

| Column | Meaning |
| --- | --- |
| `x` | Predictor value |
| `y` | Simulated response with noise |
| `truth` | Known curve before noise |

The example has 200 observations and eight knots. Four knot cases share the same predictor and response values. A separate repeated-value case checks robustness.

Generate the data from the project root:

```sh
julia --startup-file=no data/generate.jl
```

Existing data are preserved. The project runners reuse them automatically. File checksums in the R reports identify the exact inputs used.
