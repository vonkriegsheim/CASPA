#!/usr/bin/env python3
"""caspa doctor — verify a CASPA installation.

Checks the Python packages, the R/Bioconductor packages, and whether
``python`` / ``snakemake`` / ``Rscript`` resolve on PATH (the Snakefile shells
out to bare ``python`` and ``Rscript``), then prints exactly what is missing.

    python caspa/doctor.py
    python caspa/doctor.py --rscript "C:\\Program Files\\R\\R-4.6.0\\bin\\Rscript.exe"

Exit code 0 if everything REQUIRED is present, 1 otherwise. PATH issues are
warnings (the packages are installed; PATH is a per-shell concern), not failures.
"""

import argparse
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

OK = "OK"
MISSING = "MISSING"
OPT = "missing (optional)"

# (import_name, install_name, required)
PY_PACKAGES = [
    ("snakemake", "snakemake", True),
    ("numpy", "numpy", True),
    ("pandas", "pandas", True),
    ("scipy", "scipy", True),
    ("sklearn", "scikit-learn", True),
    ("umap", "umap-learn", True),
    ("scanpy", "scanpy", True),
    ("harmonypy", "harmonypy", True),
    ("leidenalg", "leidenalg", True),
    ("igraph", "igraph", True),
    ("matplotlib", "matplotlib", True),
    ("seaborn", "seaborn", True),
    ("directlfq", "directlfq", True),
    ("openai", "openai", True),
    ("adjustText", "adjustText", False),
    ("anthropic", "anthropic", False),
]

# (package, required). Includes transitive deps that are easy to miss because
# the dependent package still lists by name while failing to LOAD:
#   getopt -> optparse, GSEABase -> AUCell, GO.db -> clusterProfiler.
R_PACKAGES = [
    ("scp", True), ("SingleCellExperiment", True), ("AUCell", True),
    ("GSEABase", True), ("clusterProfiler", True), ("GO.db", True),
    ("AnnotationDbi", True), ("S4Vectors", True),
    ("org.Hs.eg.db", True), ("org.Mm.eg.db", True), ("fgsea", True),
    ("IHW", True), ("data.table", True), ("ggplot2", True), ("igraph", True),
    ("msigdbr", True), ("getopt", True), ("optparse", True), ("pheatmap", True),
    ("uwot", True), ("scales", True), ("matrixStats", True),
]


def find_rscript(override=None):
    """Locate the Rscript the pipeline will use — same precedence as run.py:
    explicit override -> CASPA_RSCRIPT -> bundled <CASPA>/R -> PATH -> system R."""
    if override:
        return override if Path(override).exists() else None
    env_r = os.environ.get("CASPA_RSCRIPT", "").strip()
    if env_r and Path(env_r).exists():
        return env_r
    exe = "Rscript.exe" if os.name == "nt" else "Rscript"
    caspa_dir = Path(__file__).resolve().parent.parent
    for cand in (caspa_dir / "R" / "bin" / exe, caspa_dir / "R" / "bin" / "x64" / exe):
        if cand.exists():
            return str(cand)
    found = shutil.which("Rscript")
    if found:
        return found
    candidates = []
    for base in (r"C:\Program Files\R", r"C:\Program Files\Microsoft\R Open"):
        p = Path(base)
        if p.is_dir():
            candidates += [d / "bin" / "Rscript.exe" for d in p.glob("R-*")
                           if (d / "bin" / "Rscript.exe").exists()]
    return str(sorted(candidates)[-1]) if candidates else None


def check_python():
    print(f"Python: {sys.version.split()[0]}  ({sys.executable})")
    missing = []
    for imp, name, required in PY_PACKAGES:
        present = importlib.util.find_spec(imp) is not None
        status = OK if present else (MISSING if required else OPT)
        if not present and required:
            missing.append(name)
        print(f"  {name:16} {status}")
    return missing


def check_tools(rscript):
    """PATH resolution of the bare commands the pipeline invokes."""
    print("\nTool resolution on PATH (needed at run time):")
    resolved = {}
    for tool in ("python", "snakemake", "Rscript"):
        path = shutil.which(tool)
        resolved[tool] = path
        shown = path or (f"{rscript}  (found, but NOT on PATH)" if tool == "Rscript"
                         and rscript else MISSING)
        print(f"  {tool:10} {shown}")
    return resolved


def check_r(rscript):
    if not rscript:
        print("\nR: Rscript not found - install R from CRAN. Cannot check R packages.")
        return None
    names = [n for n, _ in R_PACKAGES]
    # Verify by LOADING, not by name: a package can be in installed.packages()
    # yet fail to load because a dependency (GSEABase, GO.db, getopt, ...) is
    # missing — a false green the old name-only check could not see.
    expr = ("load1 <- function(p) isTRUE(tryCatch(suppressWarnings(suppressMessages("
            "requireNamespace(p, quietly=TRUE))), error=function(e) FALSE)); "
            "for (p in commandArgs(TRUE)) "
            "cat(p, if (load1(p)) as.character(packageVersion(p)) else 'MISSING', '\\n')")
    try:
        out = subprocess.run([rscript, "-e", expr] + names,
                             capture_output=True, text=True, timeout=180)
    except Exception as e:  # noqa: BLE001
        print(f"\nR: failed to run Rscript ({e})")
        return None
    versions = {}
    for line in out.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            versions[parts[0]] = parts[1]
    print(f"\nR: {rscript}")
    if os.environ.get("CASPA_RSCRIPT", "").strip():
        print(f"   (pinned via CASPA_RSCRIPT = {os.environ['CASPA_RSCRIPT']})")
    missing = []
    for name, required in R_PACKAGES:
        v = versions.get(name, "MISSING")
        if v == "MISSING" and required:
            missing.append(name)
        print(f"  {name:22} {v}")
    return missing


def main():
    ap = argparse.ArgumentParser(description="Verify a CASPA installation")
    ap.add_argument("--rscript", default=None,
                    help="Path to Rscript.exe (auto-detected if omitted)")
    args = ap.parse_args()

    print("=" * 64)
    print("CASPA doctor")
    print("=" * 64)

    rscript = find_rscript(args.rscript)
    py_missing = check_python()
    tools = check_tools(rscript)
    r_missing = check_r(rscript)

    problems, warnings = [], []
    if py_missing:
        problems.append(f"Python packages missing: {', '.join(py_missing)}")
        problems.append("    fix: pip install -r requirements-windows.txt"
                        "  (or conda env create -f environment.yml)")
    if rscript is None:
        problems.append("R not found — install R from https://cran.r-project.org/bin/windows/base/")
    elif r_missing:
        problems.append(f"R packages missing: {', '.join(r_missing)}")
        problems.append("    fix: Rscript pipeline/scripts/R/install_r_packages.R")

    path_python = tools.get("python")
    if path_python and Path(path_python).resolve() != Path(sys.executable).resolve():
        warnings.append(f"bare 'python' on PATH is {path_python}, NOT the interpreter "
                        f"with CASPA's packages ({sys.executable}). The Snakefile calls "
                        "bare 'python' - put the CASPA python's dir first on PATH.")
    if not tools.get("snakemake"):
        warnings.append("snakemake not on PATH - add <miniforge>\\Scripts to PATH "
                        "(or invoke via 'python -m snakemake')")
    if not tools.get("Rscript"):
        warnings.append("Rscript not on PATH - add <R>\\bin to PATH "
                        "(the Snakefile shells out to bare 'Rscript')")

    print("\n" + "=" * 64)
    for w in warnings:
        print("WARN: " + w)
    if problems:
        print("\nRESULT: problems found\n")
        for p in problems:
            print(" - " + p)
        sys.exit(1)
    print("\nRESULT: all required dependencies present - CASPA is ready.")
    print("Run:  python caspa/run.py --workdir <workdir> --cores N")
    sys.exit(0)


if __name__ == "__main__":
    main()
