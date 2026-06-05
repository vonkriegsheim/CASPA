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

    # Rscript -> the R that actually has CASPA's packages, chosen deterministically
    # so a stray system R can't shadow it. Priority:
    #   1. CASPA_RSCRIPT env override
    #   2. R bundled next to CASPA (the Windows installer puts it at <CASPA>\R)
    #   3. whatever is on PATH, unless it's a conda R on Windows (no win-64 scp)
    #   4. highest-version native CRAN R under Program Files
    exe = "Rscript.exe" if os.name == "nt" else "Rscript"
    rscript_used = None

    env_r = os.environ.get("CASPA_RSCRIPT", "").strip()
    if env_r and os.path.exists(env_r):
        rscript_used = env_r

    if not rscript_used:
        for cand in (os.path.join(CASPA_DIR, "R", "bin", exe),
                     os.path.join(CASPA_DIR, "R", "bin", "x64", exe)):
            if os.path.exists(cand):
                rscript_used = cand
                break

    cur_r = shutil.which("Rscript", path=env.get("PATH"))
    if not rscript_used and cur_r and not (os.name == "nt" and _is_conda_r(cur_r)):
        rscript_used = cur_r

    if not rscript_used and os.name == "nt":
        native = _native_rscript_dir()
        if native:
            rscript_used = os.path.join(native, exe)

    # Make the chosen R win for the Snakefile's bare `Rscript`.
    if rscript_used:
        prepend.insert(0, os.path.dirname(rscript_used))

    env["PATH"] = os.pathsep.join(prepend) + os.pathsep + env.get("PATH", "")
    return env, python_used, rscript_used


def _reexec_with_bundled_python_if_needed():
    """If the Python running run.py lacks snakemake (e.g. a stray pre-existing
    miniforge picked up from PATH), re-exec with the CASPA-bundled Python that has
    the deps — so `python caspa/run.py ...` works from a plain prompt, not just the
    CASPA Console. Does nothing if the current Python already has snakemake."""
    import importlib.util
    if importlib.util.find_spec("snakemake") is not None:
        return
    exe = "python.exe" if os.name == "nt" else "python"
    for cand in (os.path.join(CASPA_DIR, "miniforge3", exe),
                 os.path.join(CASPA_DIR, "miniforge3", "envs", "caspa", exe)):
        if os.path.exists(cand) and os.path.abspath(cand) != os.path.abspath(sys.executable):
            try:
                ok = subprocess.run([cand, "-c", "import snakemake"],
                                    capture_output=True).returncode == 0
            except Exception:                                  # noqa: BLE001
                ok = False
            if ok:
                print(f"[caspa run] this Python has no snakemake; switching to the "
                      f"bundled Python:\n  {cand}", file=sys.stderr)
                sys.stdout.flush(); sys.stderr.flush()
                os.execv(cand, [cand, os.path.abspath(__file__), *sys.argv[1:]])
    # Fall through: main() prints a clear "install the deps" message.


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
    _reexec_with_bundled_python_if_needed()   # may replace this process; returns if current Python is fine
    args = parse_args()

    import importlib.util
    if importlib.util.find_spec("snakemake") is None:
        print("[caspa run] ERROR: 'snakemake' is not installed in this Python:")
        print(f"  {sys.executable}")
        print("  Fix one of:")
        print("   - run from the 'CASPA Console' shortcut (Start Menu -> CASPA), or")
        print(f'   - "{sys.executable}" -m pip install -r '
              f'"{os.path.join(CASPA_DIR, "requirements-windows.txt")}"')
        sys.exit(1)

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
