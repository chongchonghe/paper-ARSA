#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
if command -v lake >/dev/null 2>&1; then
  task_lake=lake
elif [ -x .tools/lake ]; then
  task_lake=.tools/lake
else
  echo 'Install Lean via elan, then run this script again.' >&2
  exit 1
fi
"$task_lake" build
"$task_lake" env lean RSLA/AxiomAudit.lean
