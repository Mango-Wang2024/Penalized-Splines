# Check 3: Combined construction and conversion

This check starts from the same raw data and spline settings in both languages:

```text
R:      shared inputs -> R basis and penalty -> R conversion
Julia:  shared inputs -> Julia basis and penalty -> Julia conversion
                                  -> Compare final representations
```

Julia constructs its own `B,S` and passes them to the existing Julia converter. R results are read only for comparison. The R construction script is reused from `spline_basis`; both mathematical Julia modules are reused without modification.

From this folder:

```sh
sh run.sh
```

For another run:

```sh
sh run.sh outputs/run2
```

| File | Purpose |
| --- | --- |
| `run.sh` | Run R construction, R conversion and the Julia workflow |
| `scripts/reference.R` | Convert the prepared R spline objects |
| `scripts/compare.jl` | Independently construct and convert in Julia, compare and export |
| `scripts/report.R` | Summarize all three checks and verify shared input fingerprints |

The combined check has seven numerical comparisons per case, 35 across five cases. It compares the full fixed-effect space and random covariance, checks curve and penalty preservation, and maps all coefficient directions into common coordinates. Transformation matrices can differ in signs or rotations while representing the same model.

See the [results](outputs/shared/RESULTS.md) and [detailed checks](outputs/shared/julia/checks.tsv).

## Outputs for Stage 2

Each case in `outputs/shared/julia` contains:

| File | Meaning |
| --- | --- |
| `y.csv` | Responses to fit |
| `X_fixed.csv` | Intercept plus unpenalized spline components |
| `Z.csv` | Penalized spline components |
| `B.csv`, `S.csv` | Julia's basis and matching penalty in original coordinates |
| `C.csv` | Map converted coefficients back to original spline coefficients |
| `knots.csv`, `H.csv`, `penalty_scale.csv` | Training construction settings |
| `model_settings.tsv` | Dimensions, shared random variance and independent residual assumptions |

CSV matrices have no headers. The main example is `uniform_k8`. No Julia coefficients or smoothing strength have been estimated.

This validates preparation for the tested single-penalty natural cubic model. It does not validate fitting or new-data predictions.

Source: [mgcv spline conversion](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html).
