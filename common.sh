# Sourced by run_all.sh and the run.sh scripts: the Julia project, the Python, and the header line that each
# run writes above its printout (date, commit of this repository, Honghu.jl version).
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PY=${PY:-python3}

ver() {   # commit of this repository, "-dirty" if the files given (relative to the current folder) have
          # uncommitted changes; "no-git" outside git
    local c
    c=$(git rev-parse --short HEAD 2>/dev/null) || { echo no-git; return; }
    git diff --quiet HEAD -- "$@" || c+=-dirty
    echo "$c"
}

honghu() {   # the Honghu.jl release and tree hash pinned in Manifest.toml
    awk '/^\[\[deps\.Honghu\]\]/ {f = 1} f && /^git-tree-sha1/ {t = $3} f && /^repo-rev/ {r = $3} f && /^$/ {exit}
         END {gsub(/"/, "", t); gsub(/"/, "", r); print r " (tree " substr(t, 1, 7) ")"}' "$ROOT/Manifest.toml"
}

header() {   # header <script> [files whose changes make the run -dirty]
    echo "# $1, $(date '+%Y-%m-%d %H:%M'), paper-ARSA $(ver "${@:2}"), Honghu.jl $(honghu)"
}
