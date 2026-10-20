#!/bin/sh
set -eu
project_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname -- "$project_dir")
if [ "$#" -gt 2 ]; then
    printf '%s\n' 'Expected an optional new output directory and report path.' >&2
    exit 1
fi
out=${1:-"$project_dir/outputs/shared"}
report=${2:-"$out/STAGE2_RESULTS.md"}
command -v Rscript >/dev/null 2>&1 || { printf '%s\n' 'Rscript is required.' >&2; exit 1; }
command -v julia >/dev/null 2>&1 || { printf '%s\n' 'Julia is required.' >&2; exit 1; }
if [ -d "$out" ] && [ -n "$(ls -A "$out")" ]; then
    printf '%s\n' 'Outputs already exist. Supply a new output directory.' >&2
    exit 1
fi
if [ -e "$report" ]; then
    printf '%s\n' 'The summary already exists. Supply a new report path.' >&2
    exit 1
fi
if [ ! -f "$root/data/inputs/cases.tsv" ]; then
    printf '%s\n' 'The shared Stage 1 data are required.' >&2
    exit 1
fi
# Instantiate the pinned local environment, without changing a global environment.
julia --startup-file=no --project="$project_dir" -e 'using Pkg; Pkg.instantiate()'
mkdir -p "$out"
if ! julia --startup-file=no --project="$project_dir" "$project_dir/test/runtests.jl" > "$out/unit_tests.txt" 2>&1; then
    cat "$out/unit_tests.txt"
    exit 1
fi
cat "$out/unit_tests.txt"
Rscript "$project_dir/scripts/reference.R" "$root/data/inputs" "$out/r_reference"
julia --startup-file=no --project="$project_dir" "$project_dir/scripts/fit.jl" "$root/data/inputs" "$out/julia"
Rscript "$project_dir/scripts/report.R" "$root/data/inputs" "$out" "$report" "$out/unit_tests.txt"
