#!/usr/bin/env bash
# Rebuild the data, the figures and Table 1 of He, Krumholz & Wibking, "An Adaptive Reduced Speed of Light
# Approximation". See README.md.
#
#   ./run_all.sh               everything: the simulations and calculations (in parallel, about 14 threads,
#                              25 min), then the figures and Table 1
#   ./run_all.sh figures       only the figures and Table 1, from the data in git (about 20 s)
#   ./run_all.sh <step> ...    only the steps named, e.g. ./run_all.sh mf_shocks (the steps are in step_cmd)
#
# Set up once: julia --project=. -e 'using Pkg; Pkg.instantiate()' and pip install -r requirements.txt.
# PY=<python> picks the Python (default python3). The figures need LaTeX (matplotlib's text.usetex).
# Each step writes its data to data/ and its printout to data/<step>.txt, after a header: the date, the
# commit of this repository ("-dirty": uncommitted changes in src/ or Manifest.toml) and the Honghu.jl
# version. A full run also keeps the printout of each step in logs/.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-python3}
mkdir -p data figures tables logs

ver() {   # commit of this repository, "-dirty" if src/ or Manifest.toml have uncommitted changes
    local c
    c=$(git rev-parse --short HEAD 2>/dev/null) || { echo no-git; return; }
    git diff --quiet HEAD -- src Manifest.toml || c+=-dirty
    echo "$c"
}
honghu() {   # the Honghu.jl release and tree hash pinned in Manifest.toml
    awk '/^\[\[deps\.Honghu\]\]/ {f = 1} f && /^git-tree-sha1/ {t = $3} f && /^repo-rev/ {r = $3} f && /^$/ {exit}
         END {gsub(/"/, "", t); gsub(/"/, "", r); print r " (tree " substr(t, 1, 7) ")"}' Manifest.toml
}
examples() {   # the examples folder of the installed Honghu.jl
    julia --project=. -e 'using Honghu; print(pkgdir(Honghu, "examples"))'
}

step_cmd() {   # step_cmd <step>: the script of the step in NAME, its command in CMD
    case $1 in
        checkerboard)        NAME="Honghu.jl examples/irsla_paper.jl, section checkerboard"             # Fig. 3
                             CMD=(julia --project=. -t 4 "$(examples)/irsla_paper.jl" --no-plot --outdir=data checkerboard) ;;
        linear_waves_rhd)    NAME="Honghu.jl examples/linear_waves_rhd.jl"                              # Fig. 5
                             CMD=(julia --project=. -t 4 "$(examples)/linear_waves_rhd.jl" --no-plot --outdir=data) ;;
        pulse_drift)         NAME="src/sim_pulse_drift.jl"; CMD=(julia --project=. src/sim_pulse_drift.jl) ;;
        mf_shocks | arsa_phi0 | diffusion_pulse | marshak_vaytet | floor_energy)
                             NAME="src/sim_$1.jl"; CMD=(julia --project=. -t 1 "src/sim_$1.jl") ;;
        rsla_accuracy | theta_scan | m1_stability | waves_heat_capacity)
                             NAME="src/calc_$1.py"; CMD=("$PY" "src/calc_$1.py") ;;
        m1_checks)           NAME="checks/m1/*.py"; CMD=(m1_checks) ;;          # App. C, by B. Wibking
        table_media | waves_rsla)                          # figure and table scripts whose numbers the paper quotes
                             NAME="src/make_$1.py"; CMD=("$PY" "src/make_$1.py") ;;
        *)                   echo "unknown step: $1 (the steps are in step_cmd in run_all.sh)" >&2; return 1 ;;
    esac
}
step() {   # step <step>: run it; header and printout in data/<step>.txt
    local NAME CMD
    step_cmd "$1"
    { echo "# $NAME, $(date '+%Y-%m-%d %H:%M'), paper-ARSA $(ver), Honghu.jl $(honghu)"; "${CMD[@]}"; } 2>&1 |
        tee "data/$1.txt"
}

m1_checks() {   # the checks of the M1 analysis; each writes its results to data/m1/
    local f
    for f in checks/m1/*.py; do
        echo "== $f"
        "$PY" "$f" || return 1
    done
}

job() {   # job <step>: run a step in the background, its printout also in logs/<step>.log
    step "$1" > "logs/$1.log" 2>&1 &
    pids+=($!); names+=("$1")
}
wait_jobs() {   # wait for the jobs started so far; stop if one failed
    local i bad=0
    for i in "${!pids[@]}"; do
        if wait "${pids[$i]}"; then echo "done: ${names[$i]}"
        else echo "FAILED: ${names[$i]} (see logs/${names[$i]}.log)"; bad=1; fi
    done
    pids=(); names=()
    [ $bad = 0 ]
}

figures() {
    step table_media > /dev/null                 # Table 1 (tables/tab_media.tex); numbers of Secs. 3.4 and 6.2
    step waves_rsla > /dev/null                  # Fig. 2; numbers of Sec. 2.4
    local f
    for f in rsla_accuracy checkerboard waves_arsa linear_waves diffusion_pulse pulse_drift mf_shocks marshak_vaytet theta_scan; do
        "$PY" "src/make_$f.py"                   # every figure of the paper
    done
    echo "figures in figures/, Table 1 in tables/"
}

case ${1:-all} in
    all)
        pids=(); names=()
        for s in checkerboard linear_waves_rhd mf_shocks arsa_phi0 diffusion_pulse marshak_vaytet \
                 rsla_accuracy theta_scan m1_stability waves_heat_capacity m1_checks pulse_drift; do
            job "$s"
        done
        wait_jobs
        job floor_energy                 # uses the references of the shocks, the pulse and the Marshak wave
        wait_jobs
        figures ;;
    figures)
        figures ;;
    *)
        for s in "$@"; do step "$s"; done ;;
esac
