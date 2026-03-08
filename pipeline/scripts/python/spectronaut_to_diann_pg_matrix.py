#!/usr/bin/env python

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

import pandas as pd


# ----------------------------
# Helpers
# ----------------------------
def write_status(status_json: str, payload: dict):
    Path(status_json).parent.mkdir(parents=True, exist_ok=True)
    with open(status_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def write_stub_pg_matrix(out_pg_matrix: str, sample_sheet: pd.DataFrame):
    annotation_cols = [
        "Protein.Group",
        "Protein.Names",
        "Genes",
        "First.Protein.Description",
        "N.Sequences",
        "N.Proteotypic.Sequences",
    ]
    run_cols = sample_sheet["sample_file"].astype(str).tolist()
    stub = pd.DataFrame(columns=annotation_cols + run_cols)
    Path(out_pg_matrix).parent.mkdir(parents=True, exist_ok=True)
    stub.to_csv(out_pg_matrix, sep="\t", index=False)


def normalize_run_from_rlabel(rlabel: str, run_prefix: str, run_ext: str) -> str:
    """
    Example:
      rlabel = "Acinar_Distal___1_s1-c7_1_9731.d"
      -> "input\\Acinar_Distal___1_s1-c7_1_9731.raw"
    """
    base = str(rlabel).strip()
    base = re.sub(r"^.*[\\/]", "", base)  # drop any path
    base = re.sub(r"\.d$", "", base, flags=re.IGNORECASE)
    return f"{run_prefix}{base}{run_ext}"


def replace_numeric_zeros_with_blank(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """
    For DIA-NN-like pg_matrix: convert 0 -> "" (blank), keep other numbers as-is,
    and also convert NaN -> "".

    We do this column-wise to avoid converting annotation columns.
    """
    out = df.copy()
    for c in cols:
        if c not in out.columns:
            continue
        # Work with numeric where possible; if not numeric, just replace literal "0" and 0.0 too.
        s = out[c]

        # If it's numeric dtype, easy:
        if pd.api.types.is_numeric_dtype(s):
            out[c] = s.mask((s == 0) | (s.isna()), "")
        else:
            # try coerce numeric; if many are numeric-looking, treat them as numeric for masking
            num = pd.to_numeric(s, errors="coerce")
            # If coercion succeeds for at least some rows, apply numeric mask where numeric exists
            # and blank out zeros and NaNs.
            if num.notna().any():
                masked = num.mask((num == 0) | (num.isna()), "")
                # Keep non-numeric original values (shouldn't happen for intensity cols, but safe)
                out[c] = masked.where(num.notna(), s)
            else:
                # Fallback: string replace
                out[c] = s.replace({0: "", "0": "", "0.0": ""}).fillna("")
    return out


# ----------------------------
# Main
# ----------------------------
def main():
    p = argparse.ArgumentParser(
        description="Convert Spectronaut (fragment long TSV) -> directLFQ -> DIA-NN-like wide protein matrix"
    )
    p.add_argument("--spectronaut-tsv", required=True, help="Spectronaut long TSV export (fragment quant scheme).")
    p.add_argument("--sample-sheet", required=True, help="CSV with at least sample_file column (DIA-NN style).")
    p.add_argument("--out-pg-matrix", required=True, help="Output DIA-NN-like wide protein matrix TSV.")
    p.add_argument("--work-prefix", required=True, help="Prefix for intermediate directLFQ files (no extension).")
    p.add_argument("--status-json", required=True, help="Status JSON output (OK/FAILED).")

    # Run naming
    p.add_argument("--run-source", default="R.Label", choices=["R.Label", "R.FileName"], help="Which column to map runs from.")
    p.add_argument("--run-prefix", default="input\\", help="Prefix for run columns, e.g. input\\")
    p.add_argument("--run-ext", default=".d", help="Extension for run columns. Default: .d (Bruker). Use .raw for Thermo.")

    # Filtering
    p.add_argument("--remove-decoys", action="store_true", help="Drop EG.IsDecoy == TRUE rows.")
    p.add_argument("--fdr-fragment", type=float, default=0.01, help="Filter FG.Qvalue <= threshold. Set to -1 to skip.")
    p.add_argument("--fdr-ep", type=float, default=-1, help="Optional EG.Qvalue <= threshold. Default: skip.")

    # directLFQ options
    p.add_argument("--min-nonan", type=int, default=1, help="directLFQ min_nonan.")
    p.add_argument("--num-cores", type=int, default=0, help="directLFQ num_cores (0 lets tool decide).")
    p.add_argument("--deactivate-normalization", action="store_true", help="Disable between-sample normalization.")
    p.add_argument("--filename-suffix", default="", help="Suffix to append to directLFQ output files (if supported).")

    # Debug
    p.add_argument("--debug", action="store_true", help="Print debug row counts and distributions.")
    args = p.parse_args()

    t0 = time.time()

    try:
        # --------------------------
        # Load sample sheet
        # --------------------------
        ss = pd.read_csv(args.sample_sheet)
        if "sample_file" not in ss.columns:
            raise ValueError("sample_sheet.csv must contain a 'sample_file' column.")
        sample_files = ss["sample_file"].astype(str).tolist()
        sample_file_set = set(sample_files)

        # pre-create stub output (non-fatal pipeline contract)
        write_stub_pg_matrix(args.out_pg_matrix, ss)

        # --------------------------
        # Load Spectronaut TSV
        # --------------------------
        df = pd.read_csv(args.spectronaut_tsv, sep="\t", low_memory=False)

        # Required columns for fragment export scheme
        required_cols = [
            args.run_source,
            "PG.ProteinGroups",
            "PG.ProteinAccessions",
            "PG.ProteinDescriptions",
            "PG.Genes",
            "PEP.StrippedSequence",
            "PEP.IsProteotypic",
            "EG.IsDecoy",
            "EG.ModifiedSequence",
            "FG.Charge",
            "FG.IntMID",
            "FG.Qvalue",
            "FG.Quantity",
        ]
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns in Spectronaut TSV: {missing}")

        if args.debug:
            print("[DEBUG] running script:", __file__)
            print(f"[DEBUG] loaded rows: {df.shape[0]:,}")

        # NOTE: DO NOT filter on E.Errors.
        if args.debug and "E.Errors" in df.columns:
            vc = df["E.Errors"].value_counts(dropna=False).head(5)
            print("[DEBUG] E.Errors distribution (top 5):")
            print(vc.to_string())

        # --------------------------
        # Map runs to DIA-NN-like run column names
        # --------------------------
        df["Run"] = df[args.run_source].map(lambda x: normalize_run_from_rlabel(x, args.run_prefix, args.run_ext))

        # Keep only runs in sample sheet (ensures exact columns expected downstream)
        df = df[df["Run"].isin(sample_file_set)].copy()
        if args.debug:
            print(f"[DEBUG] after restricting to sample_sheet runs: {df.shape[0]:,} rows")
            print(f"[DEBUG] unique runs after restriction: {df['Run'].nunique()}")

        if df.shape[0] == 0:
            raise ValueError("After mapping runs and restricting to sample_sheet runs, 0 rows remain. Check run mapping / sample_sheet.")

        # --------------------------
        # Coerce numeric columns (robustness)
        # --------------------------
        df["FG.Qvalue"] = pd.to_numeric(df["FG.Qvalue"], errors="coerce")
        df["FG.Quantity"] = pd.to_numeric(df["FG.Quantity"], errors="coerce")
        df["FG.Charge"] = pd.to_numeric(df["FG.Charge"], errors="coerce")

        if args.fdr_ep is not None and args.fdr_ep >= 0:
            if "EG.Qvalue" not in df.columns:
                raise ValueError("EG.Qvalue filter requested but EG.Qvalue column not present.")
            df["EG.Qvalue"] = pd.to_numeric(df["EG.Qvalue"], errors="coerce")

        # --------------------------
        # Filtering
        # --------------------------
        if args.remove_decoys:
            df = df[df["EG.IsDecoy"] == False]
            if args.debug:
                print(f"[DEBUG] after remove-decoys: {df.shape[0]:,} rows")

        if args.fdr_fragment is not None and args.fdr_fragment >= 0:
            df = df[df["FG.Qvalue"].notna() & (df["FG.Qvalue"] <= args.fdr_fragment)]
            if args.debug:
                print(f"[DEBUG] after FG.Qvalue <= {args.fdr_fragment}: {df.shape[0]:,} rows")

        if args.fdr_ep is not None and args.fdr_ep >= 0:
            df = df[df["EG.Qvalue"].notna() & (df["EG.Qvalue"] <= args.fdr_ep)]
            if args.debug:
                print(f"[DEBUG] after EG.Qvalue <= {args.fdr_ep}: {df.shape[0]:,} rows")

        # Drop rows missing required fields for ion building / intensity
        df = df.dropna(subset=[
            "Run",
            "PG.ProteinGroups",
            "EG.ModifiedSequence",
            "FG.Charge",
            "FG.IntMID",
            "FG.Quantity",
        ])
        if args.debug:
            print(f"[DEBUG] after dropna required: {df.shape[0]:,} rows")

        df = df[df["FG.Quantity"] > 0].copy()
        if args.debug:
            print(f"[DEBUG] after FG.Quantity > 0: {df.shape[0]:,} rows")

        if df.shape[0] == 0:
            raise ValueError("0 rows remain after filtering. Check filters and required columns.")

        # --------------------------
        # Build directLFQ generic matrix (.aq_reformat.tsv)
        # --------------------------
        df["protein"] = df["PG.ProteinGroups"].astype(str)

        # Fragment-level ion id
        # Keep FG.IntMID as-is (string) because it contains sequence-like identifiers in your file
        df["ion"] = (
            df["EG.ModifiedSequence"].astype(str)
            + "_z" + df["FG.Charge"].astype(int).astype(str)
            + "_fg" + df["FG.IntMID"].astype(str)
        )

        df["intensity"] = df["FG.Quantity"].astype(float)

        df2 = (df.groupby(["protein", "ion", "Run"], as_index=False)
                 .agg(intensity=("intensity", "sum")))

        aq = df2.pivot_table(
            index=["protein", "ion"],
            columns="Run",
            values="intensity"
            # IMPORTANT: no fill_value here; missing stays NaN for directLFQ
        ).reset_index()

        aq_path = args.work_prefix + ".aq_reformat.tsv"
        Path(aq_path).parent.mkdir(parents=True, exist_ok=True)
        aq.to_csv(aq_path, sep="\t", index=False)

        if args.debug:
            print(f"[DEBUG] wrote aq_reformat: {aq_path}")
            print(f"[DEBUG] aq rows (ions): {aq.shape[0]:,}, runs: {aq.shape[1]-2:,}")

        if aq.shape[0] == 0:
            raise ValueError("aq_reformat matrix has 0 rows (no ions). Cannot run directLFQ.")

        # --------------------------
        # Run directLFQ
        # --------------------------
        import directlfq.lfq_manager as lfq_manager

        run_kwargs = dict(
            input_file=aq_path,
            min_nonan=args.min_nonan,
        )
        if args.num_cores and args.num_cores > 0:
            run_kwargs["num_cores"] = args.num_cores
        if args.filename_suffix:
            run_kwargs["filename_suffix"] = args.filename_suffix
        if args.deactivate_normalization:
            run_kwargs["deactivate_normalization"] = True

        lfq_manager.run_lfq(**run_kwargs)

        # Locate protein intensities output
        prot_path = aq_path.replace(".aq_reformat.tsv", ".protein_intensities.tsv")
        if not os.path.exists(prot_path):
            candidates = sorted(Path(Path(aq_path).parent).glob("*.protein_intensities.tsv"))
            if len(candidates) == 1:
                prot_path = str(candidates[0])
            else:
                raise FileNotFoundError(f"Could not find directLFQ protein intensities output. Looked for {prot_path}")

        prot = pd.read_csv(prot_path, sep="\t")

        # directLFQ output can be long or wide; handle both
        if set(["protein", "run", "protein_intensity"]).issubset(prot.columns):
            wide_int = prot.pivot_table(index="protein", columns="run", values="protein_intensity").reset_index()
            wide_int = wide_int.rename(columns={"protein": "Protein.Group"})
        else:
            wide_int = prot.rename(columns={prot.columns[0]: "Protein.Group"})

        # Ensure all sample sheet runs exist as columns (even if missing entirely)
        for run in sample_files:
            if run not in wide_int.columns:
                wide_int[run] = pd.NA

        # --------------------------
        # Build annotations (preserve exact ';' strings and repeats)
        # --------------------------
        meta_base = df[[
            "PG.ProteinGroups",
            "PG.ProteinAccessions",
            "PG.ProteinDescriptions",
            "PG.Genes",
            "PEP.StrippedSequence",
            "PEP.IsProteotypic",
        ]].copy()

        meta_base["Protein.Group"] = meta_base["PG.ProteinGroups"].astype(str)
        meta_base["PEP.IsProteotypic"] = meta_base["PEP.IsProteotypic"].astype(str).str.upper().isin(["TRUE", "1"])

        n_seq = (meta_base.groupby("Protein.Group")["PEP.StrippedSequence"]
                 .nunique(dropna=True)
                 .rename("N.Sequences")
                 .reset_index())

        n_prot = (meta_base[meta_base["PEP.IsProteotypic"]]
                  .groupby("Protein.Group")["PEP.StrippedSequence"]
                  .nunique(dropna=True)
                  .rename("N.Proteotypic.Sequences")
                  .reset_index())

        rep = (meta_base.groupby("Protein.Group", as_index=False)
               .agg({
                   "PG.ProteinAccessions": "first",
                   "PG.ProteinDescriptions": "first",
                   "PG.Genes": "first",
               }))

        rep["Protein.Names"] = rep["PG.ProteinAccessions"].fillna("").astype(str)
        rep["Genes"] = rep["PG.Genes"].fillna("").astype(str)
        rep["First.Protein.Description"] = (
            rep["PG.ProteinDescriptions"].fillna("").astype(str).str.split(";").str[0].str.strip()
        )

        meta = rep[["Protein.Group", "Protein.Names", "Genes", "First.Protein.Description"]]
        meta = meta.merge(n_seq, on="Protein.Group", how="left").merge(n_prot, on="Protein.Group", how="left")
        meta[["N.Sequences", "N.Proteotypic.Sequences"]] = meta[["N.Sequences", "N.Proteotypic.Sequences"]].fillna(0).astype(int)

        # --------------------------
        # Merge + write DIA-NN-like wide table
        # --------------------------
        out = meta.merge(wide_int, on="Protein.Group", how="right")

        annotation_cols = [
            "Protein.Group",
            "Protein.Names",
            "Genes",
            "First.Protein.Description",
            "N.Sequences",
            "N.Proteotypic.Sequences",
        ]

        # Keep only expected run columns (in sample_sheet order)
        out = out[annotation_cols + sample_files]

        # Convert NaN -> blank, and 0 -> blank (to match DIA-NN pg_matrix style)
        out = out.fillna("")
        out = replace_numeric_zeros_with_blank(out, sample_files)

        Path(args.out_pg_matrix).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out_pg_matrix, sep="\t", index=False)

        write_status(args.status_json, {
            "tool": "spectronaut_to_diann_pg_matrix",
            "status": "OK",
            "spectronaut_tsv": args.spectronaut_tsv,
            "sample_sheet": args.sample_sheet,
            "out_pg_matrix": args.out_pg_matrix,
            "aq_reformat_tsv": aq_path,
            "directlfq_protein_intensities_tsv": prot_path,
            "n_rows_spectronaut_after_filters": int(df.shape[0]),
            "n_proteins_out": int(out.shape[0]),
            "n_runs_out": int(len(sample_files)),
            "runtime_seconds": round(time.time() - t0, 3),
        })

    except Exception as e:
        # Non-fatal: keep stub output + status FAILED
        try:
            ss = pd.read_csv(args.sample_sheet)
            write_stub_pg_matrix(args.out_pg_matrix, ss)
        except Exception:
            pass

        write_status(args.status_json, {
            "tool": "spectronaut_to_diann_pg_matrix",
            "status": "FAILED",
            "error": str(e),
            "runtime_seconds": round(time.time() - t0, 3),
        })
        return 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
