#!/usr/bin/env bash
# Run mf_shocks.jl and keep its printout -- the derived parameters of every run (h, P, phi_*,
# tau_c, the eq:tau_one threshold, cells at phi = 1) and the errors -- in data/mf_shocks.txt,
# so that they can be read without a rerun. About 10 min.
set -euo pipefail
cd "$(dirname "$0")"
source common.sh
{
    header mf_shocks.jl mf_shocks.jl
    julia --project="$ROOT" -t 1 mf_shocks.jl
} 2>&1 | tee data/mf_shocks.txt
