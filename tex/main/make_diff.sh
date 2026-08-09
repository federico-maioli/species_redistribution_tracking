#!/usr/bin/env bash
# Build a coloured latexdiff PDF (main-d.pdf) between main_old.tex (previous
# submission) and main.tex (current revision).
#
# Usage (from tex/main/):
#   ./make_diff.sh
#
# Notes:
#   --allow-spaces lets latexdiff diff *inside* long text commands whose
#   arguments span paragraphs (matmethods).
#   --preamble=diff_preamble.tex injects a fix so PNAS's \matmethods does not
#   silently overwrite the added block with the deleted one when latexdiff
#   emits two \matmethods{} calls.
#   Only 'matmethods' is added to --append-textcmd: diffing inside \title{}
#   or \significancestatement{} breaks PNAS's own macros that re-expand those
#   arguments.

set -euo pipefail
cd "$(dirname "$0")"

latexdiff \
  --allow-spaces \
  --append-textcmd="matmethods" \
  --preamble=diff_preamble.tex \
  main_old.tex main.tex > main-d.tex

latexmk -pdf -interaction=nonstopmode main-d.tex
