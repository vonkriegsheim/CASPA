#!/usr/bin/env python3
"""caspa/assign_batches.py — Auto-assign batch numbers from Evosep file naming.

A batch is a single contiguous acquisition session on the Evosep. The same plate
(S{N}) can produce multiple batches if cells were loaded and run on different days,
or when consecutive 96-well plates are run back-to-back with almost no gap.

Two detectors (both configurable in evosep_batch_patterns.json):

  PRIMARY — run-number gap:
    Two consecutive runs from the same plate are a new batch when their injection
    counter gap exceeds `gap_threshold` (default 20). Within a session cells are
    sequential (gap 1–6); between sessions the gap is typically >50.

  SECONDARY — well-position drop:
    Within an already-contiguous session (gap < gap_threshold), a new batch is
    also triggered when the 96-well position drops by more than `well_drop_threshold`
    (default 80 positions). This catches back-to-back full-plate reloads where the
    well resets from H12→A1 (drop 95) with only a 1–2 run gap between plates.
    The threshold of 80 avoids false splits from small drops like A10→A1 (drop 9)
    caused by a few test injections before the main plate run.

Usage:
    python caspa/assign_batches.py --workdir /path/to/experiment
    python caspa/assign_batches.py --workdir /path/to/experiment --dry-run
    python caspa/assign_batches.py --workdir /path/to/experiment --gap-threshold 50
"""

import argparse
import json
import os
import re
import sys

CASPA_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_PATTERNS = os.path.join(CASPA_DIR, "caspa", "evosep_batch_patterns.json")

_WELL_ROW = {c: i for i, c in enumerate("abcdefgh")}


def well_to_int(well_str):
    """Row-major 96-well position: A1=1 … H12=96."""
    m = re.match(r"([a-h])(\d{1,2})$", well_str.strip().lower())
    if not m:
        return -1
    return _WELL_ROW[m.group(1)] * 12 + int(m.group(2))


def int_to_well(n):
    if n < 1:
        return "?"
    row_i, col_i = divmod(n - 1, 12)
    if row_i > 7:
        return "?"
    return f"{'ABCDEFGH'[row_i]}{col_i + 1}"


# ---------------------------------------------------------------------------
# Pattern loading and sample parsing
# ---------------------------------------------------------------------------

def load_config(patterns_path):
    with open(patterns_path, encoding="utf-8") as fh:
        data = json.load(fh)
    pats = [p for p in data.get("patterns", []) if not p.get("disabled")]
    if not pats:
        sys.exit(f"ERROR: no active patterns in {patterns_path}")
    gap_threshold       = data.get("gap_threshold", 20)
    well_drop_threshold = data.get("well_drop_threshold", 80)
    return pats, gap_threshold, well_drop_threshold


def parse_sample(sample_id, patterns):
    """Try each pattern; return (pat_name, plate_id, well_int, run_number) or Nones."""
    for pat in patterns:
        m = re.search(pat["regex"], sample_id, re.IGNORECASE)
        if m:
            plate = m.group(pat["group_plate"]).lower()
            well  = well_to_int(m.group(pat["group_well"]))
            run   = int(m.group(pat["group_run"]))
            return pat["name"], plate, well, run
    return None, None, None, None


# ---------------------------------------------------------------------------
# Batch assignment by run-number gap
# ---------------------------------------------------------------------------

def compute_batches(tagged, gap_threshold, well_drop_threshold):
    """
    tagged: list of (row, pat_name, plate, well_int, run_number)

    Sort all samples by run number. For each plate, detect a new acquisition
    session via two criteria (both configurable):
      1. Run-number gap > gap_threshold (primary)
      2. Well-position drop > well_drop_threshold within a contiguous block (secondary)
    Global batch numbers are assigned in first-appearance order.

    Returns {id(row): global_batch_int}
    """
    sorted_rows = sorted(tagged, key=lambda x: x[4])   # sort by run_number

    plate_session    = {}   # plate → current session index (1-based)
    plate_last_run   = {}   # plate → last run number seen
    plate_last_well  = {}   # plate → last well_int seen (for secondary detector)
    session_order    = []   # (plate, session) in first-appearance order
    session_seen     = set()
    row_to_key       = {}   # id(row) → (plate, session)

    for row, _pat, plate, well, run in sorted_rows:
        if plate not in plate_session:
            plate_session[plate]   = 1
            plate_last_run[plate]  = run
            plate_last_well[plate] = well
        else:
            run_gap   = run  - plate_last_run[plate]
            well_drop = plate_last_well[plate] - well   # positive when well went back

            if run_gap > gap_threshold:
                # Primary: large injection-counter gap → new session
                plate_session[plate] += 1
            elif well_drop > well_drop_threshold:
                # Secondary: well position reset within a contiguous block →
                # back-to-back plate reload (e.g. H12→A1, drop=95)
                plate_session[plate] += 1

            plate_last_run[plate]  = run
            plate_last_well[plate] = well

        key = (plate, plate_session[plate])
        row_to_key[id(row)] = key
        if key not in session_seen:
            session_seen.add(key)
            session_order.append(key)

    batch_map = {k: i + 1 for i, k in enumerate(session_order)}
    return {rid: batch_map[key] for rid, key in row_to_key.items()}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def assign_batches(manifest_path, patterns, gap_threshold, well_drop_threshold, dry_run=False):
    with open(manifest_path, encoding="utf-8") as fh:
        lines = fh.readlines()
    if not lines:
        sys.exit(f"ERROR: manifest is empty: {manifest_path}")

    header = lines[0].rstrip("\n").split("\t")
    rows   = [l.rstrip("\n").split("\t") for l in lines[1:] if l.strip()]

    if "sample_id" not in header:
        sys.exit(f"ERROR: 'sample_id' column not found in {manifest_path}")
    sid_idx = header.index("sample_id")

    if "batch" in header:
        batch_idx = header.index("batch")
    else:
        header.append("batch")
        batch_idx = len(header) - 1

    # Parse
    tagged    = []   # (row, pat_name, plate, well_int, run_number)
    unmatched = []

    for row in rows:
        sample_id = row[sid_idx] if sid_idx < len(row) else ""
        pat_name, plate, well, run = parse_sample(sample_id, patterns)
        if plate is None:
            unmatched.append((row, sample_id))
        else:
            tagged.append((row, pat_name, plate, well, run))

    if not tagged:
        sys.exit("ERROR: no samples matched any Evosep pattern. Check --patterns file.")

    # Compute batches
    row_to_batch = compute_batches(tagged, gap_threshold, well_drop_threshold)

    # Build per-batch summary for reporting
    batch_info = {}   # batch_num → {plate, pat, wells, runs, count}
    for row, pat_name, plate, well, run in sorted(tagged, key=lambda x: x[4]):
        b = row_to_batch[id(row)]
        if b not in batch_info:
            batch_info[b] = dict(plate=plate, pat=pat_name, wells=[], runs=[], count=0)
        batch_info[b]["wells"].append(well)
        batch_info[b]["runs"].append(run)
        batch_info[b]["count"] += 1

    # Report
    print(f"[assign_batches] Manifest           : {manifest_path}")
    print(f"[assign_batches] Gap threshold      : {gap_threshold}")
    print(f"[assign_batches] Well-drop threshold: {well_drop_threshold}")
    print(f"[assign_batches] Detected {len(batch_info)} batch(es):\n")
    print(f"  {'Batch':>5}  {'Plate':<6}  {'Pattern':<18}  "
          f"{'Run range':>14}  {'Well range':<13}  Samples")
    print(f"  {'-'*5}  {'-'*6}  {'-'*18}  {'-'*14}  {'-'*13}  {'-'*7}")
    for b in sorted(batch_info):
        info  = batch_info[b]
        runs  = info["runs"]
        wells = [w for w in info["wells"] if w > 0]
        run_range  = f"{min(runs)}\u2013{max(runs)}"
        well_range = (f"{int_to_well(min(wells))}\u2013{int_to_well(max(wells))}"
                      if wells else "?")
        print(f"  {b:>5}  {info['plate']:<6}  {info['pat']:<18}  "
              f"{run_range:>14}  {well_range:<13}  {info['count']}")

    if unmatched:
        print(f"\n[assign_batches] WARNING: {len(unmatched)} sample(s) did not match "
              f"any pattern \u2192 assigned batch 0")
        for _, sid in unmatched[:10]:
            print(f"    {sid}")
        if len(unmatched) > 10:
            print(f"    ... and {len(unmatched) - 10} more")

    if dry_run:
        print("\n[assign_batches] Dry run \u2014 manifest not modified.")
        return

    # Write updated manifest
    unmatched_ids = {id(row) for row, _ in unmatched}
    out_lines = ["\t".join(header) + "\n"]
    for row in rows:
        while len(row) < len(header):
            row.append("")
        row[batch_idx] = ("0" if id(row) in unmatched_ids
                          else str(row_to_batch[id(row)]))
        out_lines.append("\t".join(row) + "\n")

    with open(manifest_path, "w", encoding="utf-8") as fh:
        fh.writelines(out_lines)
    print(f"\n[assign_batches] Written: {manifest_path}")


def parse_args():
    ap = argparse.ArgumentParser(
        description="Auto-assign Evosep batch numbers to config/ms_inputs.tsv"
    )
    ap.add_argument("--workdir", required=True,
                    help="CASPA workdir (must contain config/ms_inputs.tsv)")
    ap.add_argument("--manifest", default="config/ms_inputs.tsv",
                    help="Manifest path relative to workdir (default: config/ms_inputs.tsv)")
    ap.add_argument("--patterns", default=DEFAULT_PATTERNS,
                    help=f"JSON patterns file (default: {DEFAULT_PATTERNS})")
    ap.add_argument("--gap-threshold", type=int, default=None,
                    help="Override gap_threshold from JSON (default 20). "
                         "New batch when run-number gap between consecutive same-plate "
                         "samples exceeds this value.")
    ap.add_argument("--well-drop-threshold", type=int, default=None,
                    help="Override well_drop_threshold from JSON (default 80). "
                         "Secondary detector: new batch when well position drops by "
                         "more than this within a contiguous run block "
                         "(catches H12->A1 back-to-back reloads, drop=95).")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print detected batches without modifying the manifest")
    return ap.parse_args()


def main():
    args = parse_args()
    workdir = os.path.abspath(args.workdir)
    manifest_path = os.path.join(workdir, args.manifest)
    if not os.path.isfile(manifest_path):
        sys.exit(f"ERROR: manifest not found: {manifest_path}")

    patterns, gap_threshold, well_drop_threshold = load_config(args.patterns)
    if args.gap_threshold is not None:
        gap_threshold = args.gap_threshold
    if args.well_drop_threshold is not None:
        well_drop_threshold = args.well_drop_threshold

    assign_batches(manifest_path, patterns, gap_threshold, well_drop_threshold,
                   dry_run=args.dry_run)


if __name__ == "__main__":
    main()
