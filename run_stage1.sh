#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
if [ "$#" -gt 1 ]; then
    printf '%s\n' 'Usage: sh run_stage1.sh [new-output-directory]' >&2
    exit 1
fi
if [ "$#" -eq 1 ]; then
    basis="$1/basis"
    workflow="$1/workflow"
    report="$1/STAGE1_RESULTS.md"
else
    basis="$root/spline_basis/outputs/shared"
    workflow="$root/spline_workflow/outputs/shared"
    report="$root/STAGE1_RESULTS.md"
fi
command -v Rscript >/dev/null 2>&1 || { printf '%s\n' 'Rscript is required.' >&2; exit 1; }
command -v julia >/dev/null 2>&1 || { printf '%s\n' 'Julia is required.' >&2; exit 1; }
# Refuse to overwrite either comparison before starting a new run.
for out in "$basis" "$workflow"; do
    if [ -d "$out" ] && [ -n "$(ls -A "$out")" ]; then
        printf '%s\n' 'Outputs already exist. Supply a new output directory.' >&2
        exit 1
    fi
done
if [ -e "$report" ]; then
    printf '%s\n' 'The summary already exists. Supply a new output directory.' >&2
    exit 1
fi
julia --startup-file=no "$root/spline_basis/test/runtests.jl"
julia --startup-file=no "$root/spline_workflow/test/runtests.jl"
sh "$root/spline_basis/run.sh" "$basis"
sh "$root/spline_workflow/run.sh" "$workflow"
Rscript "$root/spline_workflow/scripts/report.R" "$root/data/inputs" "$basis" "$workflow" "$report"
