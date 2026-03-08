#!/usr/bin/env python3
"""caspa init — scaffold a new CASPA workdir.

Usage:
    python caspa/init.py \\
        --workdir /path/to/MyExperiment \\
        --pg-matrix /path/to/report.pg_matrix.tsv \\
        --species human \\
        --name "My Experiment"

    # OR for Spectronaut input:
    python caspa/init.py \\
        --workdir /path/to/MyExperiment \\
        --spectronaut-tsv /path/to/export.tsv \\
        --species human \\
        --name "My Experiment"
"""

import argparse
import json
import os
import sys


CASPA_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE_PATH = os.path.join(CASPA_DIR, "config", "caspa.json.template")

SAMPLE_SHEET_HEADER = "sample_id\tsample_file\tbatch\n"


def parse_args():
    p = argparse.ArgumentParser(description="Scaffold a new CASPA workdir")
    p.add_argument("--workdir", required=True, help="Path for the new experiment workdir")
    p.add_argument("--name", default="", help="Experiment name (for caspa.json project.name)")
    p.add_argument("--species", default="human", choices=["human", "mouse"],
                   help="Species label (default: human)")
    p.add_argument("--description", default="",
                   help="One-sentence experiment description for LLM context")
    # Input — mutually exclusive
    inp = p.add_mutually_exclusive_group(required=True)
    inp.add_argument("--pg-matrix", metavar="PATH",
                     help="Path to pre-exported pg_matrix.tsv (DIA-NN or FragPipe)")
    inp.add_argument("--spectronaut-tsv", metavar="PATH",
                     help="Path to Spectronaut long-format TSV export")
    return p.parse_args()


def infer_sample_ids_from_pg_matrix(pg_matrix_path):
    """Read the header row and return all sample columns (after Protein.Group etc.)."""
    with open(pg_matrix_path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    # DIA-NN pg_matrix: first cols are Protein.Group, Protein.Ids, Protein.Names, Genes, First.Protein.Description
    skip_cols = {"Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                 "First.Protein.Description", "Proteotypic"}
    samples = [c for c in header if c not in skip_cols]
    return samples


def build_sample_sheet(samples):
    """Build ms_inputs.tsv content: sample_id | sample_file | batch."""
    lines = [SAMPLE_SHEET_HEADER]
    for s in samples:
        lines.append(f"{s}\t{s}\t1\n")
    return "".join(lines)


def main():
    args = parse_args()
    workdir = os.path.abspath(args.workdir)
    config_dir = os.path.join(workdir, "config")

    os.makedirs(config_dir, exist_ok=True)
    print(f"[caspa init] Scaffolding workdir: {workdir}")

    # Load template
    with open(TEMPLATE_PATH) as fh:
        cfg = json.load(fh)

    # Fill in project fields
    cfg["project"]["name"] = args.name or os.path.basename(workdir)
    cfg["project"]["species_label"] = args.species
    if args.description:
        cfg["project"]["description"] = args.description

    # Fill in input fields
    if args.pg_matrix:
        pg_abs = os.path.abspath(args.pg_matrix)
        cfg["input"]["pg_matrix"] = pg_abs
        cfg["input"]["spectronaut_tsv"] = None
        # Auto-populate sample sheet from pg_matrix headers
        try:
            samples = infer_sample_ids_from_pg_matrix(pg_abs)
            sheet_content = build_sample_sheet(samples)
            print(f"[caspa init] Detected {len(samples)} sample columns from pg_matrix header")
        except Exception as e:
            print(f"[caspa init] Warning: could not read pg_matrix header ({e}). "
                  "Writing empty sample sheet — fill it in manually.")
            sheet_content = SAMPLE_SHEET_HEADER
    else:
        cfg["input"]["pg_matrix"] = "TODO: /path/to/report.pg_matrix.tsv"
        cfg["input"]["spectronaut_tsv"] = os.path.abspath(args.spectronaut_tsv)
        sheet_content = SAMPLE_SHEET_HEADER
        print("[caspa init] Spectronaut input mode. "
              "Edit config/ms_inputs.tsv to add sample_id / batch mapping.")

    # Write caspa.json
    caspa_json_path = os.path.join(config_dir, "caspa.json")
    with open(caspa_json_path, "w") as fh:
        json.dump(cfg, fh, indent=2)
    print(f"[caspa init] Written: {caspa_json_path}")

    # Write ms_inputs.tsv
    sheet_path = os.path.join(config_dir, "ms_inputs.tsv")
    with open(sheet_path, "w") as fh:
        fh.write(sheet_content)
    print(f"[caspa init] Written: {sheet_path}")

    print()
    print("Next steps:")
    print(f"  1. Edit {caspa_json_path}")
    print(f"     - Set scp.llm.api_key (OpenAI key for cell type annotation)")
    print(f"     - Adjust scp parameters if needed")
    print(f"  2. Review {sheet_path}")
    print(f"     - Confirm sample_id values match pg_matrix column names exactly")
    print(f"     - Set batch column if samples span multiple acquisition runs")
    print(f"  3. Run CASPA:")
    print(f"     python {os.path.join(CASPA_DIR, 'caspa', 'run.py')} --workdir {workdir} --cores 30")


if __name__ == "__main__":
    main()
