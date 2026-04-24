#!/bin/bash
# render.sh
#
# Full publication-quality render of the current paper version.
# Produces PDF with:
#   - Full citation hyperlinking (biblatex \textcite)
#   - Section cross-references fully hyperlinked (cleveref \Cref)
#   - Word (.docx) output for Cees
#
# Usage (from the paper/ directory):
#   bash render.sh
#   bash render.sh ElasticBounds-r_v6b          # optional: specify version stem
#
# Why three steps?
#   Quarto resolves @sec-... cross-references before any filter runs,
#   producing Section~\ref{} in the .tex. The sed pass replaces these
#   with \Cref{} so cleveref can hyperlink the full "Section X" text.

set -e  # exit on error

STEM="${1:-ElasticBounds-r_v6b}"
QMD="${STEM}.qmd"
TEX="${STEM}.tex"

echo "==> Rendering ${QMD}..."
quarto render "${QMD}" --to pdf
quarto render "${QMD}" --to docx

echo "==> Applying cleveref section-reference fix..."
sed -i '' 's/Section~\\ref{\(sec-[^}]*\)}/\\Cref{\1}/g' "${TEX}"

echo "==> Recompiling LaTeX for full section hyperlinking..."
xelatex -interaction=nonstopmode "${TEX}" > /dev/null

echo ""
echo "Done. Output files:"
echo "  ${STEM}.pdf   (push to GitHub -> Overleaf pulls)"
echo "  ${STEM}.docx  (send to Cees)"
