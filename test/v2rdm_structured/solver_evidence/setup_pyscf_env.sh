#!/usr/bin/env bash
set -euo pipefail

UV="${UV:-uv}"
VENV="${VENV:-$HOME/.venvs/pyscf-h}"

"$UV" venv "$VENV"
"$UV" pip install --python "$VENV/bin/python" pyscf numpy

echo "PySCF venv ready at $VENV"
echo "Use: export PYSCF_PYTHON=$VENV/bin/python"
