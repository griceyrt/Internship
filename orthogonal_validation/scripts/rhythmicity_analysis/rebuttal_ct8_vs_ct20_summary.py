#!/usr/bin/env python3
"""
Rebuttal point (d) -- final summary numbers, WT CT8 vs CT20. Combines
rebuttal_ct8_vs_ct20_de.py/.R and rebuttal_ct8_vs_ct20_splicing.sh outputs
with Fig 1G / Fig 3B source data: how many of Fig 1G's 501 rhythmic genes are
significant DE CT8 vs CT20, how many of the 185 splicing-only genes have
significant differential splicing, and the Fig 3B-style events-vs-transcripts
comparison (PerKO vs CT8-vs-CT20).
Author: Gricey

Prefers the R DESeq2 output if present, falls back to pydeseq2.

Prerequisites:
    python3 scripts/rebuttal_ct8_vs_ct20_de.py   (or Rscript rebuttal_ct8_vs_ct20_de.R, preferred)
    bash    scripts/rebuttal_ct8_vs_ct20_splicing.sh

Run from: orthogonal_validation/
    python3 scripts/rebuttal_ct8_vs_ct20_summary.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" / "rhythmicity_analysis" / "rebuttal_ct8_vs_ct20"

DPSI_CUTOFF = 0.1
PVAL_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.5
FDR_CUTOFF = 0.05

EVENT_TYPES = ["SE", "A3", "A5", "MX", "RI", "AF", "AL"]


def load_de_results() -> tuple[pd.DataFrame, str]:
    r_path = RESULTS_DIR / "deseq2_R_tx_level_ct8_vs_ct20.csv"
    py_path = RESULTS_DIR / "deseq2_tx_level_ct8_vs_ct20.csv"
    if r_path.exists():
        df = pd.read_csv(r_path)
        return df, "R DESeq2 (deseq2_R_tx_level_ct8_vs_ct20.csv)"
    if py_path.exists():
        df = pd.read_csv(py_path, index_col=0).reset_index(names="transcript_id")
        return df, "pydeseq2 (deseq2_tx_level_ct8_vs_ct20.csv)"
    raise FileNotFoundError(
        "No DE results found -- run rebuttal_ct8_vs_ct20_de.py or rebuttal_ct8_vs_ct20_de.R first"
    )


def load_diffsplice_events(diff_dir: Path, pattern: str) -> pd.DataFrame:
    frames = []
    for et in EVENT_TYPES:
        path = diff_dir / pattern.format(event=et)
        df = pd.read_csv(path, sep="\t", header=None, names=["event_id", "dpsi", "pval"], skiprows=1)
        df["event_type"] = et
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["gene_id"] = out["event_id"].str.split(";").str[0]
    out["dpsi"] = pd.to_numeric(out["dpsi"], errors="coerce")
    out["pval"] = pd.to_numeric(out["pval"], errors="coerce")
    return out.dropna(subset=["dpsi", "pval"])


def main():
    t2 = pd.read_excel(DATA_DIR / "tables" / "Table2-significantly_cycling_transcripts.xlsx")
    t5 = pd.read_excel(DATA_DIR / "tables" / "Table5-rhythmic_AS_events.xlsx")
    rhythmic_genes = set(t2["gene_id.x"].dropna())
    as_genes = set(t5["geneID"].dropna())
    overlap = rhythmic_genes & as_genes
    splicing_only = as_genes - rhythmic_genes
    print(f"Fig 1G check: rhythmic-only={len(rhythmic_genes - as_genes)} (expect 457), "
          f"overlap={len(overlap)} (expect 44), splicing-only={len(splicing_only)} (expect 184-185)")

    de, de_source = load_de_results()
    print(f"Using DE results from: {de_source}")
    de_valid = de.dropna(subset=["padj"])
    de_sig = de_valid[(de_valid["log2FoldChange"].abs() > LOG2FC_CUTOFF) & (de_valid["padj"] < FDR_CUTOFF)]

    tested_rhythmic = rhythmic_genes & set(de_valid["gene_id"].dropna())
    sig_rhythmic = de_sig[de_sig["gene_id"].isin(rhythmic_genes)]["gene_id"].nunique()
    print(f"\n[Placeholder 1] Of {len(rhythmic_genes)} rhythmic genes, {len(tested_rhythmic)} testable, "
          f"{sig_rhythmic} significant DE CT8 vs CT20 (padj<{FDR_CUTOFF} & |log2FC|>{LOG2FC_CUTOFF})")

    ct_events = load_diffsplice_events(RESULTS_DIR / "diff", "diff_{event}.dpsi.temp.0")
    ct_sig = ct_events[(ct_events["dpsi"].abs() > DPSI_CUTOFF) & (ct_events["pval"] < PVAL_CUTOFF)]

    tested_splicing_only = splicing_only & set(ct_events["gene_id"].unique())
    sig_splicing_only = ct_sig[ct_sig["gene_id"].isin(splicing_only)]["gene_id"].nunique()
    print(f"[Placeholder 2] Of {len(splicing_only)} splicing-only genes, {len(tested_splicing_only)} testable, "
          f"{sig_splicing_only} significant differential splicing CT8 vs CT20 "
          f"(|dPSI|>{DPSI_CUTOFF} & p<{PVAL_CUTOFF})")

    perko_events = load_diffsplice_events(DATA_DIR / "nanopore_suppa", "res_{event}.dpsi.temp.0")
    perko_sig_events = perko_events[(perko_events["dpsi"].abs() > DPSI_CUTOFF) & (perko_events["pval"] < PVAL_CUTOFF)]
    perko_sig_tx = pd.read_excel(DATA_DIR / "Table3-differential_expressiong_WTvsPerKO.xlsx", sheet_name="perKO_sig")

    ct_sig_tx = de_sig

    print("\n[Fig 3B comparison -- significant splicing EVENTS vs significant DE TRANSCRIPTS]")
    print(f"  PerKO (WT-CT20 vs KO-CT20):  {len(perko_sig_events)} events / {len(perko_sig_tx)} transcripts "
          f"= {len(perko_sig_events)/len(perko_sig_tx):.2f}x")
    print(f"  CT8 vs CT20 (this analysis): {len(ct_sig)} events / {len(ct_sig_tx)} transcripts "
          f"= {len(ct_sig)/len(ct_sig_tx):.2f}x")


if __name__ == "__main__":
    main()
