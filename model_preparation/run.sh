#!/bin/sh
set -eu
project_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname -- "$project_dir")
shared="$root/data/inputs"
out=${1:-"$project_dir/outputs/shared"}
command -v Rscript >/dev/null 2>&1 || { printf '%s\n' 'Rscript is required.' >&2; exit 1; }
command -v julia >/dev/null 2>&1 || { printf '%s\n' 'Julia is required.' >&2; exit 1; }
if [ ! -f "$shared/cases.tsv" ]; then
    julia --startup-file=no "$root/data/generate.jl" "$shared"
fi
# Reuse the R construction script, then convert its prepared R objects.
Rscript "$root/spline_basis/scripts/reference.R" "$shared" "$out/r_reference"
Rscript "$project_dir/scripts/reference.R" "$out/r_reference"
Rscript "$project_dir/scripts/inspect_gamm_inputs.R" "$out/r_reference/lme_inputs" "$shared/uniform_k8/data.csv" "$shared/uniform_k8/knots.csv"
julia --startup-file=no "$project_dir/scripts/compare.jl" "$shared" "$out/r_reference" "$out/julia"
