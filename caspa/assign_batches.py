#!/usr/bin/env python3
"""caspa/assign_batches.py — Auto-assign batch numbers from Evosep file naming.

Reads config/ms_inputs.tsv, detects plate/slot identifiers in sample_id using
regex patterns defined in caspa/evosep_batch_patterns.json, and overwrites the
batch column with consecutive integers sorted by plate/slot number.

Patterns are tried in priority order (first match wins per sample):
  1. plate_keyword — 'plate1', 'plate2' etc. → physical plate = batch
  2. evosep_slot   — 's1-a1', 's2-h12' etc. → Evosep sample list = batch

Usage:
    python caspa/assign_batches.py --workdir /path/to/MyExperiment
    python caspa/assign_batches.py --workdir /path/to/MyExperiment --dry-run
    python caspa/assign_batches.py --workdir /path/to/MyExperiment \\
        --patterns /path/to/my_patterns.json
"""

import argparse
import json
import os
import re
import sys

CASPA_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_PATTERNS = os.path.join(CASPA_DIR, "caspa", "evosep_batch_patterns.json")


# ---------------------------------------------------------------------------
# Pattern loading
# ---------------------------------------------------------------------------

def load_patterns(patterns_path):
    with open(patterns_path, encoding="utf-8") as fh:
        data = json.load(fh)
    pats = [p for p in data.get("patterns", []) if not p.get("disabled")]
    if not pats:
        sys.exit(f"ERROR: no active patterns found in {patterns_path}")
    return pats


# ---------------------------------------------------------------------------
# Batch key extraction
# ---------------------------------------------------------------------------

def extract_batch_key(sample_id, patterns):
    """Return (pattern_name, raw_key, sort_key) for the first matching pattern,
    or (None, None, None) if no pattern matches.

    raw_key  — lowercased match string used as the unique batch identifier
    sort_key — integer extracted from trailing digits of raw_key, for ordering
    """
    for pat in patterns:
        m = re.search(pat["regex"], sample_id, re.IGNORECASE)
        if m:
            group = pat.get("group", 0)
            raw = m.group(group)
            digits = re.search(r"(\d+)$", raw)
            sort_key = int(digits.group(1)) if digits else 0
            return pat["name"], raw.lower(), sort_key
    return None, None, None


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------

def assign_batches(manifest_path, patterns, dry_run=False):
    with open(manifest_path, encoding="utf-8") as fh:
        lines = fh.readlines()

    if not lines:
        sys.exit(f"ERROR: manifest is empty: {manifest_path}")

    header = lines[0].rstrip("\n").split("\t")
    rows = [l.rstrip("\n").split("\t") for l in lines[1:] if l.strip()]

    if "sample_id" not in header:
        sys.exit(f"ERROR: 'sample_id' column not found in {manifest_path}")

    sid_idx = header.index("sample_id")

    # Detect or add batch column
    if "batch" in header:
        batch_idx = header.index("batch")
    else:
        header.append("batch")
        batch_idx = len(header) - 1

    # Extract batch key per row
    tagged = []
    unmatched = []
    for row in rows:
        sample_id = row[sid_idx] if sid_idx < len(row) else ""
        pat_name, raw_key, sort_key = extract_batch_key(sample_id, patterns)
        tagged.append((row, pat_name, raw_key, sort_key))
        if raw_key is None:
            unmatched.append(sample_id)

    # Build sorted unique-key → batch-number mapping
    seen = {}
    for _, _, raw_key, sort_key in tagged:
        if raw_key is not None and raw_key not in seen:
            seen[raw_key] = sort_key

    ordered_keys = sorted(seen.keys(), key=lambda k: seen[k])
    key_to_batch = {k: i + 1 for i, k in enumerate(ordered_keys)}

    # Report
    print(f"[assign_batches] Manifest : {manifest_path}")
    print(f"[assign_batches] Patterns : {DEFAULT_PATTERNS}")
    print(f"[assign_batches] Detected {len(key_to_batch)} batch(es):")
    for raw_key, batch_num in sorted(key_to_batch.items(), key=lambda x: x[1]):
        count = sum(1 for _, _, rk, _ in tagged if rk == raw_key)
        pat_name = next((pn for _, pn, rk, _ in tagged if rk == raw_key), "?")
        print(f"  batch {batch_num:>2}  '{raw_key}'  ({count} samples)  [{pat_name}]")

    if unmatched:
        print(f"\n[assign_batches] WARNING: {len(unmatched)} sample(s) did not match any pattern → assigned batch 0")
        for s in unmatched[:15]:
            print(f"    {s}")
        if len(unmatched) > 15:
            print(f"    ... and {len(unmatched) - 15} more")

    if dry_run:
        print("\n[assign_batches] Dry run — manifest not modified.")
        return

    # Write updated manifest
    out_lines = ["\t".join(header) + "\n"]
    for row, _, raw_key, _ in tagged:
        while len(row) < len(header):
            row.append("")
        row[batch_idx] = str(key_to_batch.get(raw_key, 0))
        out_lines.append("\t".join(row) + "\n")

    with open(manifest_path, "w", encoding="utf-8") as fh:
        fh.writelines(out_lines)

    print(f"\n[assign_batches] Written: {manifest_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description="Auto-assign Evosep batch numbers to config/ms_inputs.tsv"
    )
    ap.add_argument("--workdir", required=True,
                    help="CASPA workdir (must contain config/ms_inputs.tsv)")
    ap.add_argument("--manifest", default="config/ms_inputs.tsv",
                    help="Manifest path relative to workdir (default: config/ms_inputs.tsv)")
    ap.add_argument("--patterns", default=DEFAULT_PATTERNS,
                    help="JSON file with batch detection patterns "
                         f"(default: {DEFAULT_PATTERNS})")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print detected batches without modifying the manifest")
    return ap.parse_args()


def main():
    args = parse_args()
    workdir = os.path.abspath(args.workdir)
    manifest_path = os.path.join(workdir, args.manifest)

    if not os.path.isfile(manifest_path):
        sys.exit(f"ERROR: manifest not found: {manifest_path}")

    patterns = load_patterns(args.patterns)
    assign_batches(manifest_path, patterns, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
