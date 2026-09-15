#!/bin/bash
# Build the Python environment the notebooks need, as a plain venv.
#
# Not conda: `mamba env create` solves fine but dies during file linking on this
# home filesystem (numpy failed to link, transaction rolled back, ~12 min lost).
# pip wheels into a venv do the same job in about a minute.
#
#   scripts/viz-env.sh && ~/.venvs/viz/bin/jupyter lab
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="${VIZ_VENV:-$HOME/.venvs/viz}"

python3 -m venv "$VENV"
"$VENV/bin/pip" -q install --upgrade pip wheel
"$VENV/bin/pip" install -q -r "$HERE/../notebooks/requirements.txt"
"$VENV/bin/python" -m ipykernel install --user --name viz --display-name "viz (gallrna)"
"$VENV/bin/python" -c 'import pandas, matplotlib, seaborn, scipy; print("viz env ready:", "$VENV")'
