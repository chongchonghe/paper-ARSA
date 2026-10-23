#!/usr/bin/env bash
# Copy the figures and Table 1 into the paper's LaTeX source, ../writing-latex-v2 (PAPER=<dir> for another).
set -euo pipefail
cd "$(dirname "$0")"
PAPER=${PAPER:-../writing-latex-v2}
cp figures/*.pdf "$PAPER/figures/"
cp tables/tab_media.tex "$PAPER/"
echo "copied figures/*.pdf and tables/tab_media.tex to $PAPER"
