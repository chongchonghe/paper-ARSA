#!/usr/bin/env bash
# Run run.jl (the simulations; it writes the data files next to it) and save its printout in
# metadata.txt, after a header with the date, the commit of this repository ("-dirty": uncommitted
# changes to run.jl) and the version of Honghu.jl.
set -euo pipefail
cd "$(dirname "$0")"
source ../common.sh
{
    header "$(basename "$PWD")/run.jl" run.jl
    julia --project="$ROOT" -t 1 run.jl
} 2>&1 | tee metadata.txt
