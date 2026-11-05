# Direct fitting versus extra QR rotation

We compared fitting the prepared matrices directly with applying QR rotation before fitting. Both used the same inputs, variance structure and REML settings. The five cases each have 200 observations and two predictors, with 30 repetitions per method and case.

## Runtime

<a href="qr_runtime_comparison.png"><img src="qr_runtime_comparison.png" alt="Median fitting time for direct and QR fitting in five cases, with intervals showing the middle 50 percent of runs" width="720"></a>

Points show median time; lines cover the middle 50% of runs. **Further left means faster.** QR had 1–19% higher median runtime, although several ranges overlap. Timing includes model construction and fitting; the common spline-construction step is excluded.

## Memory allocation

<a href="qr_allocation_comparison.png"><img src="qr_allocation_comparison.png" alt="Median allocated memory per fit for direct and QR fitting in five cases" width="720"></a>

Points show median allocated memory per fit, in KiB (1 KiB = 1,024 bytes). **Further left means fewer allocated bytes.** QR allocated 18–24% more bytes. This measures total allocation during a fit, rather than peak RAM usage.

Both methods passed the R agreement checks and gave practically identical fitted results. The extra QR step provided no measured benefit in these cases, so we use direct fitting.

<details>
<summary>Supporting data and measurement settings</summary>

[Summary](../spline_fitting/outputs/shared/comparison/qr_benchmark_summary.tsv), [individual runs](../spline_fitting/outputs/shared/comparison/qr_benchmark_runs.tsv), [agreement checks](../spline_fitting/outputs/shared/comparison/qr_benchmark_agreement.tsv), and [settings](../spline_fitting/outputs/shared/comparison/qr_benchmark_settings.tsv) are saved in this repository. Measurements used an Apple M3, Julia 1.13.1 and MixedModels 5.9.0, with one Julia thread and one BLAS thread. Both fitting paths were warmed before timing. The [plot generator](plot_qr_comparison.R) reads these saved results; no model is refitted to draw the figures.

</details>
