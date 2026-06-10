#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="/Library/TeX/texbin:$PATH"

TEXFILE="${1:-$SCRIPT_DIR/instructions/dephasing_diagnostics.tex}"
if [[ ! -f "$TEXFILE" && -f "$SCRIPT_DIR/$TEXFILE" ]]; then
  TEXFILE="$SCRIPT_DIR/$TEXFILE"
fi

TEXDIR="$(cd "$(dirname "$TEXFILE")" && pwd)"
TEXBASE="$(basename "$TEXFILE")"
OUTDIR="$TEXDIR/build"

mkdir -p "$OUTDIR"

/Library/TeX/texbin/latexmk \
  -pdf \
  -outdir="$OUTDIR" \
  -interaction=nonstopmode \
  -halt-on-error \
  "$TEXDIR/$TEXBASE"

/Library/TeX/texbin/latexmk \
  -c \
  -outdir="$OUTDIR" \
  "$TEXDIR/$TEXBASE"
