#!/usr/bin/env bash
# Install CASPA on Linux / macOS via conda (environment.yml), then verify.
#   bash install_unix.sh
#
# On these platforms bioconda provides the full Python + R/Bioconductor stack,
# so a single conda environment is enough (unlike Windows — see install_windows.ps1).
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found. Install miniforge:" >&2
  echo "  https://github.com/conda-forge/miniforge/releases" >&2
  exit 1
fi

echo "=== [1/2] Creating/updating the 'caspa' conda environment ==="
if conda env list | grep -qE '^caspa\s'; then
  conda env update -n caspa -f "$repo/environment.yml" --prune
else
  conda env create -f "$repo/environment.yml"
fi

echo "=== [2/2] Verifying (caspa doctor) ==="
conda run -n caspa python "$repo/caspa/doctor.py"

echo
echo "Done. Activate and run with:"
echo "  conda activate caspa"
echo "  python caspa/run.py --workdir <workdir> --cores N"
