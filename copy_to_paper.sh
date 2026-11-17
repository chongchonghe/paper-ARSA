#!/usr/bin/env bash
# Copy the figures and Table 1 into the paper's LaTeX source, ../overleaf (PAPER=<dir> for another).
set -euo pipefail
cd "$(dirname "$0")"
PAPER=${PAPER:-../overleaf}
cp figures/*.pdf "$PAPER/figures/"
cp tables/tab_media.tex "$PAPER/"
echo "copied figures/*.pdf and tables/tab_media.tex to $PAPER"
