# Stage 1: Prepare the spline model

Stage 1 has three separate validation checks. All use shared data generated once in Julia by `data/generate.jl`.

| Folder | Question |
| --- | --- |
| `data` | What observations and spline settings are shared? |
| `spline_basis` | Does Julia construct the same B and S as R? |
| `spline_conversion` | Do both languages correctly convert the same R-generated B and S? |
| `spline_workflow` | Do complete Julia and R construction/conversion workflows agree? |

Requirements: Julia, R, mgcv and nlme. No extra Julia packages or Unfold installation are required.

Run all checks from this folder:

```sh
sh run_stage1.sh
```

The runner checks the mathematical examples, runs all three comparisons and verifies the shared input fingerprints. It writes [STAGE1_RESULTS.md](STAGE1_RESULTS.md). Each folder stores its own outputs in `outputs/shared`.

Existing outputs are preserved. For another run:

```sh
sh run_stage1.sh outputs/run2
```

This writes `basis`, `conversion` and `workflow` results and a summary inside the chosen directory. Each folder also has a runner for its own comparison.

The reusable mathematical implementations are `spline_basis/src/NaturalCubic.jl` and `spline_conversion/src/SplineConversion.jl`. The combined check reuses them without changing their mathematics.

The main [Julia inputs for Stage 2](spline_workflow/outputs/shared/julia/uniform_k8) contain `y`, `X_fixed` and `Z`, with the matching spline matrices and training settings. They have not been fitted.

See the [data guide](data/README.md), [basis guide](spline_basis/README.md), [conversion guide](spline_conversion/README.md) and [combined-workflow guide](spline_workflow/README.md).
