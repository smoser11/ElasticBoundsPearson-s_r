#!/bin/bash
# cleveref_fix.sh — Quarto post-render script
#
# Quarto sets QUARTO_PROJECT_OUTPUT_FILES to the rendered output path(s).
# We only act when a PDF was rendered (keep-tex: true produces a .tex alongside it).
# The sed pass converts Quarto's Section~\ref{sec-...} → \Cref{sec-...} so
# that cleveref can hyperlink the full "Section X.X" text, then xelatex
# recompiles to apply the change.

case "${QUARTO_PROJECT_OUTPUT_FILES}" in
  *.pdf*) ;;
  *)      exit 0 ;;
esac

PDF=$(echo "${QUARTO_PROJECT_OUTPUT_FILES}" | tr '\n' ' ' | grep -o '[^ ]*\.pdf' | head -1)
STEM="${PDF%.pdf}"
TEX="${STEM}.tex"

if [ ! -f "$TEX" ]; then
  echo "cleveref_fix.sh: no .tex found for ${PDF}, skipping." >&2
  exit 0
fi

echo "==> cleveref fix: patching ${TEX} ..."
sed -i '' 's/Section~\\ref{\(sec-[^}]*\)}/\\Cref{\1}/g' "${TEX}"

echo "==> cleveref fix: recompiling with xelatex ..."
xelatex -interaction=nonstopmode "${TEX}" > /dev/null 2>&1

echo "==> cleveref fix: done."
