#!/usr/bin/env python3
"""caspa run — invoke Snakemake for a CASPA workdir.

Usage:
    python caspa/run.py --workdir /path/to/MyExperiment --cores 30
    python caspa/run.py --workdir /path/to/MyExperiment --dry-run
    python caspa/run.py --workdir /path/to/MyExperiment --target scp_llm_annotation

The Snakefile shells out to bare `python` and `Rscript`, which trust PATH. PATH
order varies between shells (a conda env can put a packageless conda R first),
which is a common cause of "works for me / fails for them". To make this
deterministic, run.py arranges the subprocess PATH so:
  * bare `python`  -> the interpreter running run.py (has the CASPA python deps), and
  * bare `Rscript` -> a native CRAN R (where install_r_packages.R puts scp/AUCell),
                      unless the Rscript already on PATH is itself native.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys


CASPA_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAKEFILE  = os.path.join(CASPA_DIR, "Snakefile")


def _is_conda_r(path):
    """Heuristic: does this Rscript path look like a conda/bioconda R? Those lack
    CASPA's Bioconductor packages on Windows (no win-64 bioconda builds)."""
    if not path:
        return False
    p = path.replace("\\", "/").lower()
    return any(t in p for t in ("miniforge", "miniconda", "anaconda", "mambaforge",
                                "mamba", "/conda", "conda3", "/envs/", "/pkgs/"))


def _native_rscript_dir():
    """bin/ of a native (CRAN) R install — where install_r_packages.R installs
    scp/AUCell. Returns the highest-version one, or None."""
    bins = []
    if os.name == "nt":
        for base in (r"C:\Program Files\R", r"C:\Program Files\Microsoft\R Open"):
            for exe in glob.glob(os.path.join(base, "R-*", "bin", "Rscript.exe")):
                bins.append(os.path.dirname(exe))
    else:
        for c in ("/usr/local/bin", "/usr/bin", "/opt/homebrew/bin",
                  "/Library/Frameworks/R.framework/Resources/bin"):
            if os.path.exists(os.path.join(c, "Rscript")) and not _is_conda_r(c):
                bins.append(c)
    bins.sort(reverse=True)        # highest R-x.y.z dir first (Windows)
    return bins[0] if bins else None


def _build_env():
    """Return (env, python_used, rscript_used) with PATH arranged so the
    Snakefile's bare `python`/`Rscript` resolve to the right interpreters."""
    env = dict(os.environ)
    prepend = []

    # python -> this interpreter (it has the CASPA python deps)
    py_dir = os.path.dirname(os.path.abspath(sys.executable))
    prepend.append(py_dir)
    if os.name == "nt":
        prepend.append(os.path.join(py_dir, "Scripts"))
    python_used = sys.executable

    # Rscript -> native CRAN R. Only an issue on Windows: bioconda has no win-64
    # Bioconductor, so a conda R lacks scp/AUCell. On Linux/macOS the conda R is
    # correct (that's the supported conda path), so leave PATH as-is there.
    cur_r = shutil.which("Rscript", path=env.get("PATH"))
    rscript_used = cur_r
    if os.name == "nt":
        native = _native_rscript_dir()
        if native and (cur_r is None or _is_conda_r(cur_r)):
            prepend.insert(0, native)
            rscript_used = os.path.join(native, "Rscript.exe")

    env["PATH"] = os.pathsep.join(prepend) + os.pathsep + env.get("PATH", "")
    return env, python_used, rscript_used


def parse_args():
    p = argparse.ArgumentParser(description="Run the CASPA SCP pipeline")
    p.add_argument("--workdir", required=True,
                   help="Path to the experiment workdir (must contain config/caspa.json)")
    p.add_argument("--cores", type=int, default=8,
                   help="Number of CPU cores (default: 8)")
    p.add_argument("--dry-run", "-n", action="store_true",
                   help="Snakemake dry-run: print rules without executing")
    p.add_argument("--keep-going", "-k", action="store_true",
                   help="Continue running independent rules after a failure")
    p.add_argument("--rerun-incomplete", action="store_true", default=True,
                   help="Re-run rules with incomplete output (default: True)")
    p.add_argument("--target", metavar="RULE_OR_FILE", default=None,
                   help="Run up to a specific rule or output file instead of 'all'")
    p.add_argument("--snakemake-args", nargs=argparse.REMAINDER, default=[],
                   help="Additional Snakemake arguments (pass after --)")
    return p.parse_args()


def main():
    args = parse_args()
    workdir = os.path.abspath(args.workdir)

    # Sanity check
    caspa_json = os.path.join(workdir, "config", "caspa.json")
    if not os.path.isfile(caspa_json):
        print(f"[caspa run] ERROR: config/caspa.json not found in {workdir}")
        print(f"  Run: python {os.path.join(CASPA_DIR, 'caspa', 'init.py')} --workdir {workdir} ...")
        sys.exit(1)

    # Invoke snakemake via THIS interpreter so it uses the python that has the deps.
    cmd = [
        sys.executable, "-m", "snakemake",
        "--snakefile", SNAKEFILE,
        "--directory", workdir,
        "--cores", str(args.cores),
    ]

    if args.dry_run:
        cmd.append("--dry-run")
    if args.keep_going:
        cmd.append("--keep-going")
    if args.rerun_incomplete:
        cmd.append("--rerun-incomplete")
    if args.target:
        cmd.append(args.target)
    cmd.extend(args.snakemake_args)

    env, python_used, rscript_used = _build_env()

    print(f"[caspa run] Workdir : {workdir}")
    print(f"[caspa run] Python  : {python_used}")
    print(f"[caspa run] Rscript : {rscript_used or 'NOT FOUND — install R from CRAN and run install_r_packages.R'}")
    if os.name == "nt" and rscript_used and _is_conda_r(rscript_used):
        print("[caspa run] WARNING: Rscript resolves to a conda R, which usually lacks "
              "scp/AUCell (no win-64 bioconda builds). The R rules will fail. Install "
              "native CRAN R and run pipeline/scripts/R/install_r_packages.R.")
    print()

    result = subprocess.run(cmd, env=env)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
