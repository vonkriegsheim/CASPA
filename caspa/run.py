#!/usr/bin/env python3
"""caspa run — invoke Snakemake for a CASPA workdir.

Usage:
    python caspa/run.py --workdir /path/to/MyExperiment --cores 30
    python caspa/run.py --workdir /path/to/MyExperiment --dry-run
    python caspa/run.py --workdir /path/to/MyExperiment --target scp_llm_annotation
"""

import argparse
import os
import subprocess
import sys


CASPA_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAKEFILE  = os.path.join(CASPA_DIR, "Snakefile")


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

    cmd = [
        "snakemake",
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

    print(f"[caspa run] Workdir : {workdir}")
    print(f"[caspa run] Command : {' '.join(cmd)}")
    print()

    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
