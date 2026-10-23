#!/usr/bin/env bash
# Run run.py (closed form; it writes systems.csv and curves.csv next to it) and save its printout in
# metadata.txt, after a header with the date and the commit of this repository ("-dirty": uncommitted
# changes to run.py or ../table_media.py).
set -euo pipefail
cd "$(dirname "$0")"
source ../common.sh
{
    header "$(basename "$PWD")/run.py" run.py ../table_media.py
    $PY run.py
} 2>&1 | tee metadata.txt
