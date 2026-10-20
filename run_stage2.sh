#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
if [ "$#" -gt 1 ]; then
    printf '%s\n' 'Expected an optional new output directory.' >&2
    exit 1
fi
if [ "$#" -eq 1 ]; then
    sh "$root/spline_fitting/run.sh" "$1" "$1/STAGE2_RESULTS.md"
else
    sh "$root/spline_fitting/run.sh" "$root/spline_fitting/outputs/shared" "$root/STAGE2_RESULTS.md"
fi
