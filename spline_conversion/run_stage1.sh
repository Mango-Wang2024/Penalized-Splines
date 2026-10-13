#!/bin/sh
set -eu
project_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
shared="$(dirname -- "$project_dir")/data/inputs"
out=${1:-"$project_dir/outputs/shared"}
if ! command -v Rscript >/dev/null 2>&1; then
    printf '%s\n' 'Rscript is unavailable. Install R, then open a new terminal.' >&2
    exit 1
fi
if ! command -v julia >/dev/null 2>&1; then
    printf '%s\n' 'Julia is unavailable on PATH.' >&2
    exit 1
fi
Rscript "$project_dir/scripts/00_check_environment.R"
if [ ! -f "$shared/cases.tsv" ]; then
    julia --startup-file=no "$(dirname -- "$project_dir")/data/generate.jl" "$shared"
fi
Rscript "$project_dir/scripts/01_inspect_gamm_inputs.R" "$out/r_reference" "$shared/uniform_k8/data.csv" "$shared/uniform_k8/knots.csv"
julia --startup-file=no "$project_dir/scripts/02_validate_conversion.jl" "$out/r_reference" "$out/julia_conversion"
