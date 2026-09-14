#!/bin/bash
# Build the eggNOG-mapper environment used by `make annotate`, as a plain venv
# (conda linking into home fails on Ceres; see viz-env.sh). DIAMOND is not installed
# here: the Makefile loads the cluster module (MOD_DIAMOND) when emapper runs.
#
# eggnog-mapper 2.1.13 pins biopython==1.76 and psutil==5.7.0, which don't install on
# current Python, so install it without its pins and add current versions of the three
# libraries it imports.
#
#   scripts/emapper-env.sh && make annotate
set -euo pipefail
VENV="${EMAPPER_VENV:-$HOME/.venvs/emapper}"

# --clear: rebuild from scratch, so an interrupted earlier build (no pip) can't stick
python3 -m venv --clear "$VENV"
"$VENV/bin/pip" -q install --upgrade pip wheel
"$VENV/bin/pip" -q install --no-deps eggnog-mapper==2.1.13
"$VENV/bin/pip" -q install biopython psutil xlsxwriter
"$VENV/bin/python" -c 'import eggnogmapper, Bio, psutil; print("emapper env ready:", eggnogmapper.__file__)'
