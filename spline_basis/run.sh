#!/bin/sh
set -eu
project_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
shared="$(dirname -- "$project_dir")/data/inputs"
out=${1:-"$project_dir/outputs/shared"}
if ! command -v Rscript >/dev/null 2>&1; then
    printf '%s\n' 'Rscript is required for the reference comparison.' >&2
    exit 1
fi
if ! command -v julia >/dev/null 2>&1; then
    printf '%s\n' 'Julia is required.' >&2
    exit 1
fi
if [ ! -f "$shared/cases.tsv" ]; then
    julia --startup-file=no "$(dirname -- "$project_dir")/data/generate.jl" "$shared"
fi
Rscript "$project_dir/scripts/reference.R" "$shared" "$out/r_reference"
julia --startup-file=no "$project_dir/scripts/compare.jl" "$shared" "$out/r_reference" "$out/julia"
