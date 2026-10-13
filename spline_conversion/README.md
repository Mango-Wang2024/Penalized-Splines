# Check 2: Conversion only

This check compares conversion in Julia with R's **mgcv**. Both converters receive the exact same R-generated basis `B` and penalty `S`. Julia basis construction is not run here.

The example reads 200 Julia-generated observations from `../data/inputs/uniform_k8/data.csv` and the matching `knots.csv`. R constructs the spline basis and penalty without generating another dataset. Julia converts them into fixed and random components and checks agreement with R.

## Current result

All 44 Julia mathematical tests, 7 R reference checks and 9 live R-Julia comparisons passed. This validates conversion for the supplied R matrices. Independent construction is checked in [spline_basis](../spline_basis/README.md), and the complete workflow in [spline_workflow](../spline_workflow/README.md). Julia fitting and final predictions remain later stages. See [validation results](VALIDATION.md).

## Requirements

- R with `mgcv` and `nlme`.
- Julia. The current Julia code uses standard libraries only.

Check the setup from the project directory:

```sh
Rscript scripts/00_check_environment.R
julia --startup-file=no test/runtests.jl
```

## Run the example

For the first run:

```sh
sh run_stage1.sh
```

Existing outputs are preserved. To repeat the comparison, choose new directories:

```sh
sh run_stage1.sh outputs/run2
```

The runner reuses the shared inputs, generating them in Julia only if missing. The R script accepts an output directory, data CSV and knot CSV. Constraints can reduce the number of basis columns.

## Code

| File | Purpose |
| --- | --- |
| `scripts/00_check_environment.R` | Check R and required packages |
| `scripts/01_inspect_gamm_inputs.R` | Record inputs to R's mixed-model fitter |
| `src/SplineConversion.jl` | Convert a basis and penalty into fixed/random components |
| `scripts/02_validate_conversion.jl` | Compare the conversion with R |
| `test/runtests.jl` | Check examples with known mathematical answers |

The R script also checks that recording the inputs leaves the reference fit unchanged. Each run records its package versions.

## Outputs

- `outputs/shared/r_reference/INPUTS_EXPLAINED.txt`: input meanings and dimensions.
- `outputs/shared/r_reference/matrix_dictionary.csv`: exported matrix descriptions.
- `outputs/shared/julia_conversion/RESULTS_EXPLAINED.txt`: numerical comparison results.

CSV files contain the matrices. RDS files preserve the R objects. The reference plot illustrates the fit.

## References

- [Research plan](RESEARCH_CONTEXT.md)
- [mgcv spline construction](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smoothCon.html)
- [mgcv conversion and prediction](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html)
- [Reference source revision](https://github.com/cran/mgcv/blob/1b6a4c8374612da27e36420b4459e93acb183f2d/R/gamm.r#L1398)

The comparison used the installed mgcv version. Compatibility with the linked historical revision has not been tested.
