#!/usr/bin/env bash
# Rebuild the data, the figures and Table 1 of He, Krumholz & Wibking, "An Adaptive Reduced Speed of Light
# Approximation". See README.md.
#
#   ./run_all.sh            the simulations and checks (in parallel, about 14 threads), then the figures
#   ./run_all.sh figures    only the figures, Table 1 and the printed numbers, from the data in git
#
# Set up once: julia --project=. -e 'using Pkg; Pkg.instantiate()' and pip install -r requirements.txt.
# PY=<python> picks the Python (default python3). The figures need LaTeX (matplotlib's text.usetex).
# Printouts of the simulations go to logs/; their data and metadata go to data/ and the fig_*/ folders.
set -euo pipefail
cd "$(dirname "$0")"
source common.sh
mkdir -p data figures tables logs

job() {   # job <name> <command ...>: start in the background, printout in logs/<name>.log
    "${@:2}" > "logs/$1.log" 2>&1 &
    pids+=($!); names+=("$1")
}
example() {   # example <driver> <name> [section ...]: a driver of Honghu.jl; header and printout in data/<name>.txt
    { header "Honghu.jl examples/$1"; julia --project=. -t 4 "$EXAMPLES/$1" --no-plot --outdir=data "${@:3}"; } \
        2>&1 | tee "data/$2.txt"
}
check() {   # wait for the jobs started so far; stop if one failed
    local i bad=0
    for i in "${!pids[@]}"; do
        if wait "${pids[$i]}"; then echo "done: ${names[$i]}"
        else echo "FAILED: ${names[$i]} (see logs/${names[$i]}.log)"; bad=1; fi
    done
    pids=(); names=()
    [ $bad = 0 ]
}

if [ "${1:-all}" != figures ]; then
    EXAMPLES=$(julia --project=. -e 'using Honghu; print(pkgdir(Honghu, "examples"))')
    pids=(); names=()
    job checkerboard example irsla_paper.jl checkerboard checkerboard        # Fig. 3
    job linear_waves example linear_waves_rhd.jl linear_waves_rhd            # Fig. 5
    job mf_shocks ./run_mf_shocks.sh                                         # Figs. 7, 8; Sec. 5.4
    job arsa_phi0 julia --project=. -t 1 arsa_phi0.jl                        # Sec. 5.4, the phi0 sweep
    job diffusion_pulse fig_diffusion_pulse/run.sh                           # Fig. 6; Sec. 5.3
    job marshak_vaytet fig_marshak_vaytet/run.sh                             # Fig. 9; Sec. 5.5
    job rsla_accuracy fig_rsla_accuracy/run.sh                               # Fig. 1
    job check_theta "$PY" check_theta.py                                     # Fig. B.10, Table B.3
    job check_m1 "$PY" check_m1.py                                           # App. B.2
    check
    job floor_energy floor_energy/run.sh     # Sec. 4; uses the references of the shock, pulse and Marshak runs
    check
fi

printout() {   # printout <script>: run a Python script; header and printout in data/<script>.txt
    { header "$1" "$1"; "$PY" "$1"; } 2>&1 | tee "data/${1%.py}.txt"
}
printout table_media.py                  # Table 1 (tables/tab_media.tex); numbers of Secs. 3.4 and 6.2
printout fig_waves_rsla.py               # Fig. 2; numbers of Sec. 2.4
printout check_waves_heat_capacity.py    # numbers of Sec. 2.4
"$PY" fig_waves.py                       # Fig. 4
"$PY" fig_tests.py checkerboard linear_waves mf_shocks   # Figs. 3, 5, 7, 8
"$PY" fig_rsla_accuracy/plot.py          # Fig. 1
"$PY" fig_diffusion_pulse/plot.py        # Fig. 6
"$PY" fig_marshak_vaytet/plot.py         # Fig. 9
"$PY" fig_theta.py                       # Fig. B.10
echo "figures in figures/, Table 1 in tables/"
